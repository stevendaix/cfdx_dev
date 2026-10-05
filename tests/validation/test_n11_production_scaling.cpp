#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_system() {
    SparseMatrix A(4, 4);
    A.push_back(0, 0, 2.0e-8);
    A.push_back(0, 1, -1.0e-8);
    A.push_back(1, 0, -3.0e3);
    A.push_back(1, 1, 4.0e3);
    A.push_back(1, 2, -1.0e3);
    A.push_back(2, 1, -2.0e-2);
    A.push_back(2, 2, 5.0e-2);
    A.push_back(2, 3, -1.0e-2);
    A.push_back(3, 2, -7.0e4);
    A.push_back(3, 3, 9.0e4);
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

int main() {
    run_case("production_scaling_improves_physical_solution_and_residual", [] {
        const auto A = make_system();
        const auto exact = make_exact_solution();
        const auto b = rhs_from(A, exact);

        Vector unscaled_solution(4, 0.0);
        LinearSolverRequest unscaled_request;
        unscaled_request.krylov = KrylovModel::GMRES;
        unscaled_request.preconditioner = PreconditionerModel::Jacobi;
        const auto unscaled = solve_linear_system(
            A, b, unscaled_solution, LinearProblemKind::General,
            unscaled_request, 200, 1.0e-10);

        Vector scaled_solution(4, 0.0);
        LinearSolverRequest scaled_request = unscaled_request;
        scaled_request.scaling = LinearScalingModel::RowColumn;
        const auto scaled = solve_linear_system(
            A, b, scaled_solution, LinearProblemKind::General,
            scaled_request, 200, 1.0e-10);

        // The unscaled solver's Krylov convergence criterion is not a
        // physical-error guarantee for this deliberately ill-conditioned
        // system. The N11 contract is that explicit scaling restores the
        // requested physical residual without changing the solution itself.
        EXPECT_TRUE(unscaled.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(unscaled.physical_residual_relative));
        EXPECT_TRUE(scaled.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(scaled.applied_scaling == LinearScalingModel::RowColumn);
        EXPECT_TRUE(std::isfinite(scaled.physical_residual_relative));
        EXPECT_TRUE(scaled.physical_residual_relative < 1.0e-10);
        EXPECT_TRUE(scaled.physical_residual_relative < unscaled.physical_residual_relative);

        for (std::size_t i = 0; i < exact.size(); ++i) {
            EXPECT_NEAR(scaled_solution(i), exact(i), 1.0e-8);
        }
    });

    run_case("default_solver_path_is_unchanged", [] {
        const auto A = make_system();
        const auto exact = make_exact_solution();
        const auto b = rhs_from(A, exact);

        Vector solution(4, 0.0);
        const auto report = solve_linear_system(
            A, b, solution, LinearProblemKind::General,
            LinearSolverRequest{}, 200, 1.0e-10);

        EXPECT_TRUE(report.applied_scaling == LinearScalingModel::None);
        EXPECT_TRUE(std::isfinite(report.physical_residual_relative));
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(report.physical_residual_relative < 1.0e-10);
    });

    return run_all();
}
