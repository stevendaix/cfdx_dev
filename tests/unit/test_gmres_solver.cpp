// M0.8-T06 — Tests for GMRES solver

#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"
#include <cmath>

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

    run_case("gmres_exact_nonsymmetric_3x3_reference_solution", []() {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 4.0); A.push_back(0, 1, 1.0);
        A.push_back(1, 0, -2.0); A.push_back(1, 1, 5.0); A.push_back(1, 2, 1.0);
        A.push_back(2, 0, 1.0); A.push_back(2, 1, -1.0); A.push_back(2, 2, 3.0);
        A.finalize();

        Vector x_exact(3);
        x_exact(0) = 1.0; x_exact(1) = -2.0; x_exact(2) = 0.5;
        const auto b_values = A.matvec(x_exact);
        Vector b(3);
        for (std::size_t i = 0; i < 3; ++i) b(i) = b_values[i];
        Vector x(3, 0.0);

        const auto result = solve_gmres(A, b, x, 3, 20, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.residual_relative < 1e-11);
        const auto Ax = A.matvec(x);
        double true_residual_sq = 0.0;
        double b_sq = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            const double r = Ax[i] - b(i);
            true_residual_sq += r * r;
            b_sq += b(i) * b(i);
            EXPECT_TRUE(std::isfinite(x(i)));
            EXPECT_NEAR(x(i), x_exact(i), 1e-10);
        }
        EXPECT_TRUE(std::sqrt(true_residual_sq / b_sq) < 1e-9);
    });

    run_case("gmres_rejects_invalid_controls", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 2.0);
        A.push_back(1, 1, 3.0);
        A.finalize();
        Vector b(2, 1.0);
        Vector x(2, 0.0);

        EXPECT_TRUE(solve_gmres(A, b, x, 0, 20, 1e-12).status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(solve_gmres(A, b, x, 2, 0, 1e-12).status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(solve_gmres(A, b, x, 2, 20, 0.0).status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_honors_fixed_restart_controls", []() {
        const std::size_t n = 12;
        auto A = make_poisson_1d(n);
        Vector b(n, 1.0);
        Vector x(n, 0.0);
        KrylovControls controls;
        controls.restart_min = static_cast<int>(n);
        controls.restart_max = static_cast<int>(n);
        controls.adaptive_restart = false;
        const auto result = solve_gmres(A, b, x, static_cast<int>(n), 50, 1e-12,
                                        nullptr, controls);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.residual_relative < 1e-10);
    });

    run_case("gmres_saddle_point_breakdown_reports_true_residual", []() {
        // Small singular saddle-point system without a pressure gauge:
        //
        // [ 1  0  1 ] [u]   [1]
        // [ 0  1  1 ] [v] = [1]
        // [-1 -1  0 ] [p]   [3]
        //
        // The first two rows imply u=v=1, while continuity then requires
        // -u-v=3, so the RHS is inconsistent with the null space.
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 1.0); A.push_back(0, 2, 1.0);
        A.push_back(1, 1, 1.0); A.push_back(1, 2, 1.0);
        A.push_back(2, 0, -1.0); A.push_back(2, 1, -1.0);
        A.finalize();

        Vector b(3);
        b(0) = 1.0; b(1) = 1.0; b(2) = 3.0;
        Vector x(3, 0.0);

        const auto result = solve_gmres(A, b, x, 3, 3, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
        EXPECT_TRUE(result.residual > 1e-8);
        EXPECT_TRUE(result.residual_relative > 1e-8);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_TRUE(std::isfinite(result.residual_relative));
    });

    run_case("gmres_nonhappy_breakdown_reports_true_residual", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.finalize();

        Vector b(2);
        b(0) = 1.0;
        b(1) = 1.0;
        Vector x(2, 0.0);

        const auto result = solve_gmres(A, b, x, 2, 2, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_TRUE(result.residual > 0.5);
        EXPECT_TRUE(std::isfinite(result.residual_relative));
        EXPECT_TRUE(result.residual_relative > 0.5);
        EXPECT_NEAR(x(0), 1.0, 1e-12);
        EXPECT_NEAR(x(1), 0.0, 1e-12);
    });

    run_case("gmres_lucky_breakdown_converges_without_nan", []() {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.push_back(2, 2, 1.0);
        A.finalize();

        Vector b(3);
        b(0) = 2.0; b(1) = -3.0; b(2) = 4.0;
        Vector x(3, 0.0);

        const auto result = solve_gmres(A, b, x, 3, 3, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations <= 3);
        EXPECT_TRUE(std::isfinite(result.residual_relative));
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(std::isfinite(x(i)));
            EXPECT_NEAR(x(i), b(i), 1e-12);
        }
    });

    return run_all();
}
