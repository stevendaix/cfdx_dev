#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/linear_solver_context.h"
#include "cfdx/core/linalg/null_space.h"
#include "common/test_harness.h"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_neumann_laplacian(std::size_t n, double scale = 1.0) {
    SparseMatrix matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        matrix.push_back(
            i, i, scale * ((i == 0 || i + 1 == n) ? 1.0 : 2.0));
        if (i > 0) matrix.push_back(i, i - 1, -scale);
        if (i + 1 < n) matrix.push_back(i, i + 1, -scale);
    }
    matrix.finalize();
    return matrix;
}

Vector multiply(const SparseMatrix& matrix, const Vector& x) {
    const auto values = matrix.matvec(x);
    Vector result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) result(i) = values[i];
    return result;
}

double true_residual(const SparseMatrix& matrix,
                     const Vector& x,
                     const Vector& rhs) {
    const auto ax = matrix.matvec(x);
    double squared = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double value = rhs(i) - ax[i];
        squared += value * value;
    }
    return std::sqrt(squared);
}

} // namespace

int main() {
    run_case("constant_projector_removes_mean", [] {
        const auto null_space = NullSpaceProjector::constant(4);
        Vector values(4);
        values(0) = 1.0;
        values(1) = 2.0;
        values(2) = 7.0;
        values(3) = 10.0;
        null_space.remove(values);
        EXPECT_NEAR(values(0) + values(1) + values(2) + values(3), 0.0, 1e-13);
        EXPECT_TRUE(null_space.component_norm(values) < 1e-13);
    });

    run_case("projector_checks_declared_modes_against_operator", [] {
        const auto singular = make_neumann_laplacian(4);
        const auto constant = NullSpaceProjector::constant(4);
        EXPECT_TRUE(constant.is_null_space(singular));
        EXPECT_TRUE(constant.operator_residual(singular) < 1e-14);

        SparseMatrix identity(4, 4);
        for (std::size_t i = 0; i < 4; ++i) identity.push_back(i, i, 1.0);
        identity.finalize();
        EXPECT_TRUE(!constant.is_null_space(identity));
        Vector rhs(4, 0.0);
        Vector solution(4, 0.0);
        EXPECT_TRUE(solve_cg(identity, rhs, solution, constant).status ==
                    SolverStatus::NOT_APPLICABLE);
    });

    run_case("basis_is_orthonormalized_and_validated", [] {
        Vector constant(4, 1.0);
        Vector alternating(4);
        alternating(0) = 1.0;
        alternating(1) = -1.0;
        alternating(2) = 1.0;
        alternating(3) = -1.0;
        const NullSpaceProjector null_space({constant, alternating});
        EXPECT_TRUE(null_space.basis_size() == 2);
        EXPECT_NEAR(null_space.basis()[0].dot(null_space.basis()[1]), 0.0, 1e-14);
        EXPECT_NEAR(null_space.basis()[0].norm2(), 1.0, 1e-14);
        EXPECT_NEAR(null_space.basis()[1].norm2(), 1.0, 1e-14);
        EXPECT_THROW(NullSpaceProjector({constant, constant}), std::invalid_argument);
        Vector wrong_dimension(3);
        EXPECT_THROW(null_space.remove(wrong_dimension), std::invalid_argument);
    });

    run_case("projected_cg_solves_singular_compatible_system", [] {
        const auto matrix = make_neumann_laplacian(5);
        Vector exact(5);
        exact(0) = -2.0;
        exact(1) = -1.0;
        exact(2) = 0.0;
        exact(3) = 1.0;
        exact(4) = 2.0;
        const Vector rhs = multiply(matrix, exact);
        const auto null_space = NullSpaceProjector::constant(5);

        Vector solution(5, 8.0);
        const auto result = solve_cg(
            matrix, rhs, solution, null_space, 100, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(matrix, solution, rhs) < 1e-11);
        EXPECT_TRUE(null_space.component_norm(solution) < 1e-12);
        EXPECT_TRUE((solution - exact).norm_inf() < 1e-10);
    });

    run_case("projected_cg_rejects_incompatible_rhs", [] {
        const auto matrix = make_neumann_laplacian(4);
        const auto null_space = NullSpaceProjector::constant(4);
        Vector rhs(4, 0.0);
        rhs(0) = 1.0;
        Vector solution(4, 0.0);
        const auto result = solve_cg(matrix, rhs, solution, null_space);
        EXPECT_TRUE(!null_space.is_compatible(rhs));
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_NEAR(solution.norm2(), 0.0, 0.0);
    });

    run_case("dispatcher_exposes_constant_pressure_null_space", [] {
        const auto matrix = make_neumann_laplacian(5);
        Vector exact(5);
        exact(0) = -1.0;
        exact(1) = -0.5;
        exact(2) = 0.0;
        exact(3) = 0.5;
        exact(4) = 1.0;
        const Vector rhs = multiply(matrix, exact);
        Vector solution(5, 3.0);

        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::Jacobi;
        request.null_space = NullSpaceModel::Constant;
        const auto report = solve_linear_system(
            matrix, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 100, 1e-12);
        EXPECT_TRUE(report.plan.null_space == NullSpaceModel::Constant);
        EXPECT_TRUE(report.plan.preconditioner == PreconditionerModel::Jacobi);
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(matrix, solution, rhs) < 1e-11);
        EXPECT_TRUE(std::string(to_string(report.plan.null_space)) == "constant");

        EXPECT_THROW(select_linear_solver(
            LinearProblemKind::Momentum, matrix.n_rows(), request),
            std::invalid_argument);
    });

    run_case("projected_cg_automatically_uses_native_amg", [] {
        constexpr std::size_t n = 65;
        const auto matrix = make_neumann_laplacian(n);
        Vector exact(n);
        for (std::size_t i = 0; i < n; ++i)
            exact(i) = static_cast<double>(i) -
                       0.5 * static_cast<double>(n - 1);
        const Vector rhs = multiply(matrix, exact);
        Vector solution(n, 4.0);

        LinearSolverRequest request;
        request.null_space = NullSpaceModel::Constant;
        const auto report = solve_linear_system(
            matrix, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 300, 1e-10);
        EXPECT_TRUE(report.plan.krylov == KrylovModel::CG);
        EXPECT_TRUE(report.plan.preconditioner == PreconditionerModel::NativeAMG);
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(matrix, solution, rhs) < 1e-8);
        EXPECT_TRUE(NullSpaceProjector::constant(n).component_norm(solution) < 1e-10);
    });

    run_case("projected_cg_supports_smoothed_aggregation_amg", [] {
        constexpr std::size_t n = 65;
        const auto matrix = make_neumann_laplacian(n);
        Vector exact(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double coordinate = static_cast<double>(i) /
                                      static_cast<double>(n - 1);
            exact(i) = std::cos(2.0 * 3.14159265358979323846 * coordinate);
        }
        NullSpaceProjector::constant(n).remove(exact);
        const Vector rhs = multiply(matrix, exact);
        Vector solution(n, -3.0);

        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::SmoothedAggregationAMG;
        request.null_space = NullSpaceModel::Constant;
        const auto report = solve_linear_system(
            matrix, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 300, 1e-10);
        EXPECT_TRUE(report.plan.preconditioner ==
                    PreconditionerModel::SmoothedAggregationAMG);
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(matrix, solution, rhs) < 1e-8);
        EXPECT_TRUE(NullSpaceProjector::constant(n).component_norm(solution) < 1e-10);
    });

    run_case("singular_amg_reuses_hierarchy_for_numeric_updates", [] {
        constexpr std::size_t n = 65;
        const auto matrix = make_neumann_laplacian(n);
        const auto scaled_matrix = make_neumann_laplacian(n, 1.25);
        const auto null_space = NullSpaceProjector::constant(n);
        Vector exact(n);
        for (std::size_t i = 0; i < n; ++i)
            exact(i) = static_cast<double>(i) -
                       0.5 * static_cast<double>(n - 1);

        NativeSmoothedAggregationAMGPreconditioner preconditioner(true);
        EXPECT_TRUE(preconditioner.uses_constant_null_space());
        ReusableCgContext context(preconditioner);
        Vector solution(n, 0.0);
        const Vector rhs = multiply(matrix, exact);
        const auto first = context.solve(
            matrix, rhs, solution, 300, 1e-10, {}, {}, &null_space);
        EXPECT_TRUE(first.status == SolverStatus::CONVERGED);

        solution.fill(2.0);
        const Vector scaled_rhs = multiply(scaled_matrix, exact);
        const auto second = context.solve(
            scaled_matrix, scaled_rhs, solution, 300, 1e-10,
            {}, {}, &null_space);
        EXPECT_TRUE(second.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(scaled_matrix, solution, scaled_rhs) < 1e-8);
        EXPECT_TRUE(context.stats().full_setups == 1);
        EXPECT_TRUE(context.stats().numeric_updates == 1);
        EXPECT_TRUE(context.stats().solves == 2);
        EXPECT_TRUE(preconditioner.hierarchy_builds() == 1);
        EXPECT_TRUE(preconditioner.numeric_updates() == 1);
    });

    run_case("singular_amg_application_is_symmetric_positive_on_projected_space", [] {
        constexpr std::size_t n = 65;
        const auto matrix = make_neumann_laplacian(n);
        const auto null_space = NullSpaceProjector::constant(n);
        Vector left(n);
        Vector right(n);
        for (std::size_t i = 0; i < n; ++i) {
            left(i) = std::sin(0.17 * static_cast<double>(i));
            right(i) = std::cos(0.11 * static_cast<double>(i));
        }
        null_space.remove(left);
        null_space.remove(right);

        for (const auto method : {NativeAMGMethod::BoomerStyle,
                                  NativeAMGMethod::SmoothedAggregation}) {
            NativeAMGPreconditioner preconditioner(method, true);
            EXPECT_TRUE(preconditioner.setup(matrix));
            Vector applied_left(n);
            Vector applied_right(n);
            EXPECT_TRUE(preconditioner.apply(left, applied_left));
            EXPECT_TRUE(preconditioner.apply(right, applied_right));
            const double lhs = left.dot(applied_right);
            const double rhs = right.dot(applied_left);
            const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
            EXPECT_NEAR(lhs, rhs, 1e-11 * scale);
            EXPECT_TRUE(left.dot(applied_left) > 0.0);
            EXPECT_TRUE(right.dot(applied_right) > 0.0);
            EXPECT_TRUE(null_space.component_norm(applied_left) < 1e-12);
            EXPECT_TRUE(null_space.component_norm(applied_right) < 1e-12);
        }
    });

    return run_all();
}
