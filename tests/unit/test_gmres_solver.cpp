// M0.8-T06 — Tests for GMRES solver

#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

SparseMatrix make_poisson_1d(std::size_t n) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

int main() {
    run_case("gmres_solve_diagonal", []() {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 2.0);
        A.push_back(1, 1, 3.0);
        A.push_back(2, 2, 4.0);
        A.finalize();

        Vector b(3);
        b(0) = 2.0; b(1) = 6.0; b(2) = 12.0;

        Vector x(3, 0.0);
        auto result = solve_gmres(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        EXPECT_TRUE(result.residual_relative < 1e-10);
        EXPECT_NEAR(x(0), 1.0, 1e-10);
        EXPECT_NEAR(x(1), 2.0, 1e-10);
        EXPECT_NEAR(x(2), 3.0, 1e-10);
    });

    run_case("gmres_solve_poisson_1d", []() {
        const std::size_t n = 5;
        auto A = make_poisson_1d(n);

        Vector b(n);
        b(0) = 0.0; b(1) = 1.0; b(2) = 1.0; b(3) = 1.0; b(4) = 0.0;

        Vector x(n, 0.0);
        auto result = solve_gmres(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.residual_relative < 1e-10);
        auto Ax = A.matvec(x);
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(Ax[i], b(i), 1e-10);
        }
    });

    run_case("gmres_already_converged", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(2, 0.0);
        Vector x(2, 0.0);
        auto result = solve_gmres(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations == 0);
    });

    run_case("gmres_dimension_mismatch", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(3);
        Vector x(2, 0.0);
        auto result = solve_gmres(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_with_restart", []() {
        const std::size_t n = 10;
        auto A = make_poisson_1d(n);

        Vector b(n);
        for (std::size_t i = 0; i < n; ++i) b(i) = 1.0;

        Vector x(n, 0.0);
        auto result = solve_gmres(A, b, x, 5);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.residual_relative < 1e-10);

        auto Ax = A.matvec(x);
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(Ax[i], b(i), 1e-10);
        }
    });

    return run_all();
}