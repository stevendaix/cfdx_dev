#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <cmath>

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

static Vector reference_solution() {
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
    for (std::size_t i = 0; i < values.size(); ++i)
        b(i) = values[i];
    return b;
}

int main() {
    run_case("production_scaling_preserves_physical_solution", [] {
        const auto A = make_system();
        const auto exact = reference_solution();
        const auto b = rhs_from(A, exact);

        LinearSolverRequest unscaled_request;
        unscaled_request.krylov = KrylovModel::GMRES;
        unscaled_request.preconditioner = PreconditionerModel::Jacobi;

        Vector unscaled_solution(4, 0.0);
        const auto unscaled = solve_linear_system(
            A, b, unscaled_solution, LinearProblemKind::General,
            unscaled_request, 200, 1.0e-11);

        EXPECT_TRUE(unscaled.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(unscaled.scaling == MatrixScaling::None);
        EXPECT_TRUE(unscaled.physical_residual_relative < 1.0e-11);

        LinearSolverRequest scaled_request = unscaled_request;
        scaled_request.scaling = MatrixScaling::RowColumn;

        Vector scaled_solution(4, 0.0);
        const auto scaled = solve_linear_system(
            A, b, scaled_solution, LinearProblemKind::General,
            scaled_request, 200, 1.0e-11);

        EXPECT_TRUE(scaled.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(scaled.scaling == MatrixScaling::RowColumn);
        EXPECT_TRUE(scaled.physical_residual_relative < 1.0e-11);
        for (std::size_t i = 0; i < 4; ++i)
            EXPECT_NEAR(scaled_solution(i), exact(i), 1.0e-9);
    });

    run_case("production_scaling_does_not_repair_singular_rows", [] {
        SparseMatrix A(2, 2);
        A.push_back(1, 1, 3.0);
        A.finalize();

        Vector b(2, 0.0);
        Vector x(2, 0.0);
        LinearSolverRequest request;
        request.krylov = KrylovModel::GMRES;
        request.preconditioner = PreconditionerModel::Jacobi;
        request.scaling = MatrixScaling::RowColumn;

        const auto report = solve_linear_system(
            A, b, x, LinearProblemKind::General, request, 20, 1.0e-11);

        EXPECT_TRUE(report.result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(std::isinf(report.physical_residual_relative));
        EXPECT_NEAR(x(0), 0.0, 0.0);
    });

    return run_all();
}
