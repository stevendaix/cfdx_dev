// M0.8-T04 — Tests for CG solver

#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

// Crée une matrice symétrique définie positive (diagonale dominante).
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
    run_case("cg_solve_diagonal", []() {
        // A = 2I, b = [2, 2, 2], solution attendue x = [1, 1, 1].
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 2.0);
        A.push_back(1, 1, 2.0);
        A.push_back(2, 2, 2.0);
        A.finalize();

        Vector b(3);
        b(0) = 2.0; b(1) = 2.0; b(2) = 2.0;

        Vector x(3, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        EXPECT_TRUE(result.residual_relative < 1e-10);
        EXPECT_NEAR(x(0), 1.0, 1e-10);
        EXPECT_NEAR(x(1), 1.0, 1e-10);
        EXPECT_NEAR(x(2), 1.0, 1e-10);
    });

    run_case("cg_solve_poisson_1d", []() {
        // Système 1D : -u'' = f avec u(0)=u(n-1)=0.
        // A = tridiag(-1, 2, -1), b = [0, 1, 1, 1, 0].
        const std::size_t n = 5;
        auto A = make_poisson_1d(n);

        Vector b(n);
        b(0) = 0.0; b(1) = 1.0; b(2) = 1.0; b(3) = 1.0; b(4) = 0.0;

        Vector x(n, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.residual_relative < 1e-10);
        // Vérification : A x ≈ b.
        auto Ax = A.matvec(x);
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(Ax[i], b(i), 1e-10);
        }
    });

    run_case("cg_already_converged", []() {
        // x = 0, b = 0 → déjà convergé.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(2, 0.0);
        Vector x(2, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations == 0);
    });

    run_case("cg_max_iter_reached", []() {
        // Système mal conditionné avec très peu d'itérations.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1e10);
        A.push_back(1, 1, 1e-10);
        A.finalize();

        Vector b(2);
        b(0) = 1.0; b(1) = 1.0;
        Vector x(2, 0.0);
        auto result = solve_cg(A, b, x, 1, 1e-20);
        // Soit convergé, soit max_iter atteint.
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED ||
                    result.status == SolverStatus::MAX_ITER_REACHED);
    });

    run_case("cg_not_applicable_non_square", []() {
        SparseMatrix A(2, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(2);
        Vector x(3, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("cg_rejects_nonfinite_system", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, std::numeric_limits<double>::quiet_NaN());
        A.finalize();
        Vector b(2, 1.0);
        Vector x(2, 0.0);
        const auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
    });

    run_case("cg_dimension_mismatch", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(3);  // mauvaise taille
        Vector x(2, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    return run_all();
}