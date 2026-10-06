#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_reference_matrix() {
    SparseMatrix A(4, 4);
    A.push_back(0, 0, 4.0);
    A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -2.0);
    A.push_back(1, 1, 5.0);
    A.push_back(1, 2, -1.0);
    A.push_back(2, 1, -3.0);
    A.push_back(2, 2, 7.0);
    A.push_back(2, 3, -2.0);
    A.push_back(3, 2, -1.0);
    A.push_back(3, 3, 3.0);
    A.finalize();
    return A;
}

static Vector make_exact_solution() {
    Vector x(4);
    x(0) = 0.75;
    x(1) = -1.25;
    x(2) = 2.0;
    x(3) = -0.5;
    return x;
}

static Vector rhs_from(const SparseMatrix& A, const Vector& x) {
    const auto values = A.matvec(x);
    Vector b(x.size());
    for (std::size_t i = 0; i < values.size(); ++i) b(i) = values[i];
    return b;
}

static SparseMatrix make_perturbed_matrix() {
    SparseMatrix A(4, 4);
    A.push_back(0, 0, 4.0 + 4.0e-10);
    A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -2.0);
    A.push_back(1, 1, 5.0);
    A.push_back(1, 2, -1.0 - 2.0e-10);
    A.push_back(2, 1, -3.0);
    A.push_back(2, 2, 7.0);
    A.push_back(2, 3, -2.0);
    A.push_back(3, 2, -1.0);
    A.push_back(3, 3, 3.0 + 3.0e-10);
    A.finalize();
    return A;
}

static SparseMatrix change_units(
    const SparseMatrix& A,
    const double row_scale,
    const double column_scale) {
    SparseMatrix transformed(A.n_rows(), A.n_cols());
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k) {
            const auto j = A.columns_data()[k];
            transformed.push_back(i, j, row_scale * A.values_data()[k] * column_scale);
        }
    }
    transformed.finalize();
    return transformed;
}

int main() {
    run_case("solution_is_invariant_under_row_and_variable_unit_changes", [] {
        const auto A = make_reference_matrix();
        const auto exact = make_exact_solution();
        const auto b = rhs_from(A, exact);

        LinearSolverRequest request;
        request.krylov = KrylovModel::GMRES;
        request.preconditioner = PreconditionerModel::Jacobi;
        request.scaling = LinearScalingModel::RowColumn;

        Vector reference_solution(4, 0.0);
        const auto reference = solve_linear_system(
            A, b, reference_solution, LinearProblemKind::General,
            request, 200, 1.0e-12);
        EXPECT_TRUE(reference.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(reference.physical_residual_relative < 1.0e-12);

        const double row_scale = 1.0e12;
        const double column_scale = 1.0e-12;
        const auto transformed = change_units(A, row_scale, column_scale);
        Vector transformed_rhs(4);
        for (std::size_t i = 0; i < 4; ++i) transformed_rhs(i) = row_scale * b(i);

        Vector transformed_solution(4, 0.0);
        const auto changed = solve_linear_system(
            transformed, transformed_rhs, transformed_solution,
            LinearProblemKind::General, request, 200, 1.0e-12);

        EXPECT_TRUE(changed.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(changed.physical_residual_relative < 1.0e-12);
        for (std::size_t i = 0; i < exact.size(); ++i)
            EXPECT_NEAR(transformed_solution(i), exact(i) / column_scale, 1.0e-8);

        for (std::size_t i = 0; i < exact.size(); ++i)
            EXPECT_NEAR(reference_solution(i), transformed_solution(i) * column_scale, 1.0e-9);
    });

    run_case("solver_is_independent_of_initial_guess", [] {
        const auto A = make_reference_matrix();
        const auto exact = make_exact_solution();
        const auto b = rhs_from(A, exact);

        LinearSolverRequest request;
        request.krylov = KrylovModel::GMRES;
        request.preconditioner = PreconditionerModel::Jacobi;
        request.scaling = LinearScalingModel::SymmetricDiagonal;

        Vector zero_guess(4, 0.0);
        const auto zero = solve_linear_system(
            A, b, zero_guess, LinearProblemKind::PressurePoisson,
            request, 200, 1.0e-12);

        Vector perturbed_guess(4);
        perturbed_guess(0) = 1.0e6;
        perturbed_guess(1) = -1.0e6;
        perturbed_guess(2) = 3.0e5;
        perturbed_guess(3) = -7.0e5;
        const auto perturbed = solve_linear_system(
            A, b, perturbed_guess, LinearProblemKind::PressurePoisson,
            request, 200, 1.0e-12);

        EXPECT_TRUE(zero.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(perturbed.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(zero.physical_residual_relative < 1.0e-12);
        EXPECT_TRUE(perturbed.physical_residual_relative < 1.0e-12);
        for (std::size_t i = 0; i < exact.size(); ++i) {
            EXPECT_NEAR(zero_guess(i), exact(i), 1.0e-9);
            EXPECT_NEAR(perturbed_guess(i), exact(i), 1.0e-9);
        }
    });

    run_case("small_coefficient_perturbation_preserves_a_bounded_physical_residual", [] {
        const auto A = make_perturbed_matrix();
        const auto exact = make_exact_solution();
        const auto b = rhs_from(A, exact);

        // Deliberately perturb several coefficients by O(1e-10). This is a
        // sensitivity probe, not a relaxed accuracy requirement.
        LinearSolverRequest request;
        request.krylov = KrylovModel::GMRES;
        request.preconditioner = PreconditionerModel::Jacobi;
        request.scaling = LinearScalingModel::RowColumn;

        Vector solution(4, 0.0);
        const auto report = solve_linear_system(
            A, b, solution, LinearProblemKind::General,
            request, 200, 1.0e-12);

        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(report.physical_residual_relative));
        EXPECT_TRUE(report.physical_residual_relative < 1.0e-12);
        for (std::size_t i = 0; i < exact.size(); ++i)
            EXPECT_NEAR(solution(i), exact(i), 1.0e-9);
    });

    return run_all();
}
