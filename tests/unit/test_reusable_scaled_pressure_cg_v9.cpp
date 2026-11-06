#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::core;

static void demand(bool value, const char* reason) {
    if (!value) throw std::runtime_error(reason);
}

int main() {
    try {
        // SPD, non-uniform diagonal scale comparable to a stretched mesh.
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 0.005);
        A.push_back(0, 1, -0.010);
        A.push_back(1, 0, -0.010);
        A.push_back(1, 1, 2.0);
        A.push_back(1, 2, -3.0);
        A.push_back(2, 1, -3.0);
        A.push_back(2, 2, 1000.0);
        A.finalize();
        Vector reference(3, 0.0);
        reference(0) = 2.0;
        reference(1) = -1.0;
        reference(2) = 0.5;
        const auto right = A.matvec(reference);
        Vector b(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) b(i) = right[i];
        JacobiPreconditioner preconditioner;
        ReusableCgContext context(preconditioner);
        Vector x(3, 0.0);
        auto first = solve_reusable_scaled_cg(
            A, b, x, context, LinearScalingModel::SymmetricDiagonal,
            100, 1e-10);
        demand(first.status == SolverStatus::CONVERGED, "scaled CG failed");
        demand(first.residual_relative <= 1e-10,
               "scaled residual incorrectly certified physical operator");
        for (std::size_t i = 0; i < 3; ++i)
            demand(std::abs(x(i) - reference(i)) < 1e-7,
                   "scaled pressure unknown was not unscaled");
        x = Vector(3, 0.0);
        auto second = solve_reusable_scaled_cg(
            A, b, x, context, LinearScalingModel::SymmetricDiagonal,
            100, 1e-10);
        demand(second.status == SolverStatus::CONVERGED,
               "second scaled CG failed");
        const auto stats = context.stats();
        demand(stats.solves == 2 && stats.full_setups == 1 &&
               stats.unchanged_reuses == 1,
               "scaled CG did not reuse its numeric preconditioner");

        // A strict tolerance below the achieved FP64 residual must never be
        // marked converged simply because the scaled system converged.
        Vector strict_x(3, 0.0);
        const auto strict = solve_reusable_scaled_cg(
            A, b, strict_x, context, LinearScalingModel::SymmetricDiagonal,
            30, 1e-30);
        demand(strict.status != SolverStatus::CONVERGED,
               "false scaled pressure certification");
        demand(std::isfinite(strict.residual_relative),
               "physical pressure failure did not report its true residual");

        // Homogeneous pressure correction: b=0 does NOT mean an infinite
        // relative residual. Normalize with the initial original-system
        // residual and report the actual nonzero absolute residual.
        SparseMatrix homogeneous(3, 3);
        homogeneous.push_back(0, 0, 2.7);
        homogeneous.push_back(0, 1, -0.33);
        homogeneous.push_back(1, 0, -0.33);
        homogeneous.push_back(1, 1, 4.1);
        homogeneous.push_back(1, 2, -0.72);
        homogeneous.push_back(2, 1, -0.72);
        homogeneous.push_back(2, 2, 3.6);
        homogeneous.finalize();
        Vector hb(3, 0.0), hx(3, 0.0);
        hx(0) = 0.27; hx(1) = -0.41; hx(2) = 1.2;
        const double initial_residual = physical_residual_norm(homogeneous, hx, hb);
        JacobiPreconditioner homogeneous_pc;
        ReusableCgContext homogeneous_context(homogeneous_pc);
        const auto homogeneous_result = solve_reusable_scaled_cg(
            homogeneous, hb, hx, homogeneous_context,
            LinearScalingModel::SymmetricDiagonal, 100, 1e-10);
        const double final_residual = physical_residual_norm(homogeneous, hx, hb);
        demand(homogeneous_result.status == SolverStatus::CONVERGED,
               "zero-RHS scaled CG falsely rejected a converged pressure correction");
        demand(final_residual <= 1e-10 * initial_residual,
               "zero-RHS scaled CG failed its original-system physical gate");
        demand(std::abs(homogeneous_result.residual - final_residual) < 1e-24,
               "homogeneous scaled CG lost the nonzero absolute residual");
        demand(std::abs(homogeneous_result.residual_relative -
                        final_residual / initial_residual) < 1e-15,
               "homogeneous scaled CG returned the wrong relative reference");

        // Independently rescaling rows/columns of an SPD matrix need not
        // preserve symmetry, so CG must reject such a request.
        bool rejected_nonsymmetric_scaling = false;
        try {
            Vector rejected(3, 0.0);
            (void)solve_reusable_scaled_cg(
                A, b, rejected, context, LinearScalingModel::RowColumn,
                20, 1e-10);
        } catch (const std::invalid_argument&) { rejected_nonsymmetric_scaling = true; }
        demand(rejected_nonsymmetric_scaling,
               "reusable CG accepted nonsymmetric algebraic scaling");
        LinearSolverRequest invalid_dispatch;
        invalid_dispatch.krylov = KrylovModel::CG;
        invalid_dispatch.preconditioner = PreconditionerModel::Jacobi;
        invalid_dispatch.scaling = LinearScalingModel::RowColumn;
        bool rejected_dispatch = false;
        try {
            Vector rejected(3, 0.0);
            (void)solve_linear_system(
                A, b, rejected, LinearProblemKind::PressurePoisson,
                invalid_dispatch, 20, 1e-10);
        } catch (const std::invalid_argument&) { rejected_dispatch = true; }
        demand(rejected_dispatch,
               "linear CG dispatch accepted nonsymmetric algebraic scaling");

        bool rejected_nullspace = false;
        try {
            NullSpaceProjector nullspace = NullSpaceProjector::constant(3);
            Vector rejected(3, 0.0);
            (void)solve_reusable_scaled_cg(
                A, b, rejected, context, LinearScalingModel::SymmetricDiagonal,
                10, 1e-10, &nullspace);
        } catch (const std::invalid_argument&) { rejected_nullspace = true; }
        demand(rejected_nullspace,
               "scaled CG silently combined inconsistent null-space projection");

        SparseMatrix singular(2, 2);
        singular.push_back(0, 0, 1.0);
        singular.finalize();
        Vector zeros(2, 0.0), rhs(2, 1.0);
        const auto unusable = solve_reusable_scaled_cg(
            singular, rhs, zeros, context, LinearScalingModel::SymmetricDiagonal,
            10, 1e-9);
        demand(unusable.status == SolverStatus::NOT_APPLICABLE,
               "zero diagonal accepted as pressure scaling");
        std::cout << "PASS scaled reusable CG / original residual / AMG context contract\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
