#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_solver_context.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/native_fieldsplit.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_poisson(std::size_t n) {
    SparseMatrix matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        matrix.push_back(i, i, 2.0);
        if (i > 0) matrix.push_back(i, i - 1, -1.0);
        if (i + 1 < n) matrix.push_back(i, i + 1, -1.0);
    }
    matrix.finalize();
    return matrix;
}

SparseMatrix make_scaled_poisson(std::size_t n, double scale) {
    SparseMatrix matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        matrix.push_back(i, i, 2.0 * scale);
        if (i > 0) matrix.push_back(i, i - 1, -scale);
        if (i + 1 < n) matrix.push_back(i, i + 1, -scale);
    }
    matrix.finalize();
    return matrix;
}

SparseMatrix make_two_field_system() {
    SparseMatrix matrix(4, 4);
    matrix.push_back(0, 0, 4.0); matrix.push_back(0, 1, 1.0);
    matrix.push_back(0, 2, 1.0);
    matrix.push_back(1, 0, 1.0); matrix.push_back(1, 1, 3.0);
    matrix.push_back(1, 3, 1.0);
    matrix.push_back(2, 0, 2.0); matrix.push_back(2, 2, 5.0);
    matrix.push_back(2, 3, 1.0);
    matrix.push_back(3, 1, 1.0); matrix.push_back(3, 2, 1.0);
    matrix.push_back(3, 3, 4.0);
    matrix.finalize();
    return matrix;
}

Vector multiply(const SparseMatrix& matrix, const Vector& x) {
    const auto values = matrix.matvec(x);
    Vector result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) result(i) = values[i];
    return result;
}

double relative_true_residual(const SparseMatrix& matrix,
                              const Vector& x,
                              const Vector& rhs) {
    const auto ax = matrix.matvec(x);
    double residual2 = 0.0;
    double rhs2 = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double residual = rhs(i) - ax[i];
        residual2 += residual * residual;
        rhs2 += rhs(i) * rhs(i);
    }
    return std::sqrt(residual2 / std::max(rhs2, 1e-300));
}

} // namespace

int main() {
    run_case("boomer_style_amg_reduces_true_residual", [] {
        const auto matrix = make_poisson(64);
        Vector rhs(64);
        for (std::size_t i = 0; i < rhs.size(); ++i)
            rhs(i) = (i % 3 == 0) ? 1.0 : -0.25;

        NativeBoomerAMGPreconditioner amg;
        Vector correction;
        EXPECT_TRUE(!amg.apply(rhs, correction));
        EXPECT_TRUE(!amg.last_error().empty());
        EXPECT_TRUE(amg.setup(matrix));
        EXPECT_TRUE(amg.is_ready());
        EXPECT_TRUE(amg.coarse_size() > 0);
        EXPECT_TRUE(amg.coarse_size() < matrix.n_rows());
        EXPECT_TRUE(amg.apply(rhs, correction));
        EXPECT_TRUE(correction.is_valid());
        EXPECT_TRUE(relative_true_residual(matrix, correction, rhs) < 1.0);
    });

    run_case("cg_accepts_native_boomer_style_amg", [] {
        const auto matrix = make_poisson(96);
        Vector exact(96);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.05 * static_cast<double>(i + 1));
        const Vector rhs = multiply(matrix, exact);

        NativeBoomerAMGPreconditioner amg;
        amg.configure(AMGMemoryPolicy::Balanced);
        Vector solution(96, 0.0);
        const auto result = solve_cg(matrix, rhs, solution, amg, 500, 1e-10);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        const double amg_true_residual =
            relative_true_residual(matrix, solution, rhs);
        EXPECT_TRUE(amg_true_residual < 1e-8);
        EXPECT_NEAR(result.residual_relative, amg_true_residual, 1e-13);

        Vector native_solution(96, 0.0);
        const auto native = solve_cg(matrix, rhs, native_solution, 500, 1e-10);
        EXPECT_TRUE(native.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, native_solution, rhs) < 1e-8);
        EXPECT_TRUE((solution - native_solution).norm_inf() < 1e-7);
    });

    run_case("cg_accepts_native_smoothed_aggregation_amg", [] {
        const auto matrix = make_poisson(96);
        const auto updated = make_scaled_poisson(96, 1.25);
        Vector exact(96);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.05 * static_cast<double>(i + 1));

        NativeSmoothedAggregationAMGPreconditioner amg;
        EXPECT_TRUE(std::string(amg.name()) ==
                    "NativeSmoothedAggregationAMG");
        EXPECT_TRUE(amg.method() == NativeAMGMethod::SmoothedAggregation);
        ReusableCgContext context(amg);

        Vector solution(96, 0.0);
        const Vector rhs = multiply(matrix, exact);
        const auto first = context.solve(
            matrix, rhs, solution, 300, 1e-10);
        EXPECT_TRUE(first.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-9);
        EXPECT_TRUE(context.stats().full_setups == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);

        solution.fill(0.0);
        const Vector updated_rhs = multiply(updated, exact);
        const auto refreshed = context.solve(
            updated, updated_rhs, solution, 300, 1e-10);
        EXPECT_TRUE(refreshed.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(
            updated, solution, updated_rhs) < 1e-9);
        EXPECT_TRUE(context.stats().numeric_updates == 1);
        EXPECT_TRUE(amg.numeric_updates() == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
    });

    run_case("reusable_gmres_updates_amg_without_rebuilding_hierarchy", [] {
        const auto matrix = make_poisson(96);
        const auto updated = make_scaled_poisson(96, 1.25);
        Vector exact(96);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::cos(0.03 * static_cast<double>(i + 1));

        NativeBoomerAMGPreconditioner amg;
        ReusableGmresContext context(amg);
        KrylovControls controls;
        controls.adaptive_restart = false;
        controls.restart_min = 24;
        controls.restart_max = 24;

        Vector solution(96, 0.0);
        const Vector rhs = multiply(matrix, exact);
        const auto first = context.solve(
            matrix, rhs, solution, 24, 300, 1e-10, controls);
        EXPECT_TRUE(first.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-9);
        EXPECT_TRUE(context.stats().full_setups == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
        const auto allocations = context.workspace_reallocations();

        solution.fill(0.0);
        const auto repeated = context.solve(
            matrix, rhs, solution, 24, 300, 1e-10, controls);
        EXPECT_TRUE(repeated.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-9);
        EXPECT_TRUE(context.stats().unchanged_reuses == 1);
        EXPECT_TRUE(context.workspace_reallocations() == allocations);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);

        solution.fill(0.0);
        const Vector updated_rhs = multiply(updated, exact);
        const auto refreshed = context.solve(
            updated, updated_rhs, solution, 24, 300, 1e-10, controls);
        EXPECT_TRUE(refreshed.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(updated, solution, updated_rhs) < 1e-9);
        EXPECT_TRUE(context.stats().numeric_updates == 1);
        EXPECT_TRUE(amg.numeric_updates() == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
    });

    run_case("reusable_cg_updates_amg_without_rebuilding_hierarchy", [] {
        const auto matrix = make_poisson(96);
        const auto updated = make_scaled_poisson(96, 1.25);
        Vector exact(96);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.02 * static_cast<double>(i + 1));

        NativeBoomerAMGPreconditioner amg;
        ReusableCgContext context(amg);

        Vector solution(96, 0.0);
        const Vector rhs = multiply(matrix, exact);
        const auto first = context.solve(matrix, rhs, solution, 300, 1e-10);
        EXPECT_TRUE(first.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-9);
        EXPECT_TRUE(context.stats().full_setups == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);

        solution.fill(0.0);
        const auto repeated = context.solve(matrix, rhs, solution, 300, 1e-10);
        EXPECT_TRUE(repeated.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-9);
        EXPECT_TRUE(context.stats().unchanged_reuses == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);

        solution.fill(0.0);
        const Vector updated_rhs = multiply(updated, exact);
        const auto refreshed = context.solve(
            updated, updated_rhs, solution, 300, 1e-10);
        EXPECT_TRUE(refreshed.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(updated, solution, updated_rhs) < 1e-9);
        EXPECT_TRUE(context.stats().numeric_updates == 1);
        EXPECT_TRUE(context.stats().solves == 3);
        EXPECT_TRUE(amg.numeric_updates() == 1);
        EXPECT_TRUE(amg.hierarchy_builds() == 1);

        const auto resized = make_poisson(64);
        Vector resized_exact(64);
        for (std::size_t i = 0; i < resized_exact.size(); ++i)
            resized_exact(i) = std::cos(0.04 * static_cast<double>(i + 1));
        const Vector resized_rhs = multiply(resized, resized_exact);
        Vector resized_solution(64, 0.0);
        const auto rebuilt = context.solve(
            resized, resized_rhs, resized_solution, 300, 1e-10);
        EXPECT_TRUE(rebuilt.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(
            resized, resized_solution, resized_rhs) < 1e-9);
        EXPECT_TRUE(context.stats().full_setups == 2);
        EXPECT_TRUE(context.stats().solves == 4);
        EXPECT_TRUE(amg.hierarchy_builds() == 2);
    });

    run_case("solver_models_are_selected_by_problem_and_can_be_overridden", [] {
        const auto pressure = select_linear_solver(
            LinearProblemKind::PressurePoisson, 96);
        EXPECT_TRUE(pressure.krylov == KrylovModel::CG);
        EXPECT_TRUE(pressure.preconditioner == PreconditionerModel::NativeAMG);

        const auto momentum = select_linear_solver(
            LinearProblemKind::Momentum, 96);
        EXPECT_TRUE(momentum.krylov == KrylovModel::BiCGStab);
        EXPECT_TRUE(momentum.preconditioner == PreconditionerModel::ILU0);

        const auto coupled = select_linear_solver(
            LinearProblemKind::CoupledPressureVelocity, 512);
        EXPECT_TRUE(coupled.krylov == KrylovModel::FGMRES);
        EXPECT_TRUE(coupled.preconditioner ==
                    PreconditionerModel::CoupledBlockSchur);

        LinearSolverRequest request;
        request.krylov = KrylovModel::FGMRES;
        request.preconditioner = PreconditionerModel::GaussSeidel;
        const auto forced = select_linear_solver(
            LinearProblemKind::Momentum, 96, request);
        EXPECT_TRUE(forced.krylov == KrylovModel::FGMRES);
        EXPECT_TRUE(forced.preconditioner == PreconditionerModel::GaussSeidel);
        EXPECT_TRUE(!forced.automatic_krylov);
        EXPECT_TRUE(!forced.automatic_preconditioner);

        request.krylov = KrylovModel::LGMRES;
        EXPECT_THROW(select_linear_solver(
            LinearProblemKind::Momentum, 96, request), std::invalid_argument);

        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::ILU0;
        EXPECT_THROW(select_linear_solver(
            LinearProblemKind::PressurePoisson, 96, request),
            std::invalid_argument);

        const auto identity = make_scalar_preconditioner(
            PreconditionerModel::None);
        EXPECT_TRUE(identity != nullptr);
        EXPECT_TRUE(std::string(identity->name()) == "Identity");

        LinearSolverRequest smoothed_request;
        smoothed_request.krylov = KrylovModel::CG;
        smoothed_request.preconditioner =
            PreconditionerModel::SmoothedAggregationAMG;
        const auto smoothed = select_linear_solver(
            LinearProblemKind::PressurePoisson, 96, smoothed_request);
        EXPECT_TRUE(smoothed.preconditioner ==
                    PreconditionerModel::SmoothedAggregationAMG);
        const auto smoothed_pc = make_scalar_preconditioner(
            PreconditionerModel::SmoothedAggregationAMG);
        EXPECT_TRUE(std::string(smoothed_pc->name()) ==
                    "NativeSmoothedAggregationAMG");

        EXPECT_THROW(select_linear_solver(
            LinearProblemKind::Momentum, 96, smoothed_request),
            std::invalid_argument);
    });

    run_case("selected_pressure_and_momentum_models_solve_their_systems", [] {
        const auto pressure_matrix = make_poisson(64);
        Vector exact(64);
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.04 * static_cast<double>(i + 1));
        const Vector pressure_rhs = multiply(pressure_matrix, exact);
        Vector pressure_solution(64, 0.0);
        const auto pressure = solve_linear_system(
            pressure_matrix, pressure_rhs, pressure_solution,
            LinearProblemKind::PressurePoisson, {}, 400, 1e-10);
        EXPECT_TRUE(pressure.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(
            pressure_matrix, pressure_solution, pressure_rhs) < 1e-9);

        SparseMatrix momentum_matrix(64, 64);
        for (std::size_t i = 0; i < 64; ++i) {
            momentum_matrix.push_back(i, i, 3.0);
            if (i > 0) momentum_matrix.push_back(i, i - 1, -1.25);
            if (i + 1 < 64) momentum_matrix.push_back(i, i + 1, -0.75);
        }
        momentum_matrix.finalize();
        const Vector momentum_rhs = multiply(momentum_matrix, exact);
        Vector momentum_solution(64, 0.0);
        const auto momentum = solve_linear_system(
            momentum_matrix, momentum_rhs, momentum_solution,
            LinearProblemKind::Momentum, {}, 400, 1e-10);
        EXPECT_TRUE(momentum.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(
            momentum_matrix, momentum_solution, momentum_rhs) < 1e-9);
    });

    run_case("fieldsplit_schur_variants_are_well_defined", [] {
        const auto matrix = make_two_field_system();
        Vector rhs(4);
        rhs(0) = 1.0; rhs(1) = -2.0; rhs(2) = 3.0; rhs(3) = 0.5;

        const std::vector<FieldSplitStrategy> strategies{
            FieldSplitStrategy::Additive,
            FieldSplitStrategy::DiagonalSchur,
            FieldSplitStrategy::LowerSchur,
            FieldSplitStrategy::UpperSchur,
            FieldSplitStrategy::FullSchur};
        for (const auto strategy : strategies) {
            NativeFieldSplitPreconditioner split({0, 1}, {2, 3}, strategy);
            Vector correction(4, 0.0);
            EXPECT_TRUE(split.setup(matrix));
            EXPECT_TRUE(split.apply(rhs, correction));
            EXPECT_TRUE(correction.is_valid());
            EXPECT_TRUE(std::string(split.name()).find("NativeFieldSplit") == 0);
            if (strategy == FieldSplitStrategy::FullSchur)
                EXPECT_TRUE(relative_true_residual(matrix, correction, rhs) < 1e-12);
        }
    });

    run_case("fieldsplit_diagonal_uses_petsc_schur_sign", [] {
        SparseMatrix matrix(2, 2);
        matrix.push_back(0, 0, 2.0);
        matrix.push_back(1, 1, 4.0);
        matrix.finalize();
        Vector rhs(2);
        rhs(0) = 6.0;
        rhs(1) = 8.0;

        NativeFieldSplitPreconditioner diagonal(
            {0}, {1}, FieldSplitStrategy::DiagonalSchur);
        Vector correction(2, 0.0);
        EXPECT_TRUE(diagonal.setup(matrix));
        EXPECT_TRUE(diagonal.apply(rhs, correction));
        EXPECT_NEAR(correction(0), 3.0, 1e-14);
        EXPECT_NEAR(correction(1), -2.0, 1e-14);

        NativeFieldSplitPreconditioner full(
            {0}, {1}, FieldSplitStrategy::FullSchur);
        EXPECT_TRUE(full.setup(matrix));
        EXPECT_TRUE(full.apply(rhs, correction));
        EXPECT_NEAR(correction(0), 3.0, 1e-14);
        EXPECT_NEAR(correction(1), 2.0, 1e-14);
    });

    run_case("gmres_uses_full_schur_fieldsplit", [] {
        const auto matrix = make_two_field_system();
        Vector exact(4);
        exact(0) = 1.0; exact(1) = -2.0; exact(2) = 0.5; exact(3) = 3.0;
        const Vector rhs = multiply(matrix, exact);

        NativeFieldSplitPreconditioner split(
            {0, 1}, {2, 3}, FieldSplitStrategy::FullSchur);
        Vector solution(4, 0.0);
        const auto result = solve_gmres(
            matrix, rhs, solution, 4, 20, 1e-12, &split);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(matrix, solution, rhs) < 1e-11);
        EXPECT_TRUE((solution - exact).norm_inf() < 1e-10);
    });

    run_case("fieldsplit_rejects_overlap_without_stale_success", [] {
        const auto matrix = make_two_field_system();
        NativeFieldSplitPreconditioner split(
            {0, 1}, {2, 3}, FieldSplitStrategy::LowerSchur);
        EXPECT_TRUE(split.setup(matrix));

        SparseMatrix wrong_size(3, 3);
        wrong_size.push_back(0, 0, 1.0);
        wrong_size.push_back(1, 1, 1.0);
        wrong_size.push_back(2, 2, 1.0);
        wrong_size.finalize();
        EXPECT_TRUE(!split.setup(wrong_size));
        EXPECT_TRUE(!split.last_error().empty());
        Vector rhs(4, 1.0), correction(4, 0.0);
        EXPECT_TRUE(!split.apply(rhs, correction));
    });

    run_case("schur_core_clears_state_after_failed_setup", [] {
        const auto matrix = make_two_field_system();
        SchurComplementPreconditioner schur(
            {0, 1}, {2, 3}, SchurFactorization::Full);
        EXPECT_TRUE(schur.setup(matrix));

        SparseMatrix nonsquare(4, 3);
        nonsquare.push_back(0, 0, 1.0);
        nonsquare.finalize();
        EXPECT_TRUE(!schur.setup(nonsquare));
        Vector rhs(4, 1.0), correction(4, 0.0);
        EXPECT_TRUE(!schur.apply(rhs, correction));
    });

    return run_all();
}
