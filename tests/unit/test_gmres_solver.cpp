// M0.8-T06 — Tests for GMRES solver

#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

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

    run_case("gmres_tiny_operator_scale_is_not_breakdown", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1e-15);
        A.finalize();
        Vector b(1);
        b(0) = 1e-15;
        Vector x(1, 0.0);
        const auto result = solve_gmres(A, b, x, 1, 2, 1e-14);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0), 1.0, 1e-12);
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
        // [-1 -1 -2 ] [p]   [0]
        //
        // The third row is exactly the negative sum of the first two rows,
        // so the matrix is singular.  With b=(1,1,0), the left-null-space
        // compatibility condition is violated: b_3 != -(b_1+b_2).
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 1.0); A.push_back(0, 2, 1.0);
        A.push_back(1, 1, 1.0); A.push_back(1, 2, 1.0);
        A.push_back(2, 0, -1.0); A.push_back(2, 1, -1.0);
        A.push_back(2, 2, -2.0);
        A.finalize();

        Vector b(3);
        b(0) = 1.0; b(1) = 1.0; b(2) = 0.0;
        Vector x(3, 0.0);

        const auto result = solve_gmres(A, b, x, 3, 3, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
        EXPECT_TRUE(result.iterations == 2);
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
        // GMRES returns the minimum-residual iterate in the reachable Krylov
        // space; for diag(1,0) this is x=(1,1), not the arbitrary null-space
        // representative (1,0).
        EXPECT_NEAR(x(1), 1.0, 1e-12);
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

    run_case("gmres_zero_rhs_converges_from_nonzero_initial_guess", []() {
        LinearOperator op;
        op.size = 1;
        op.apply = [](const Vector& in, Vector& out) {
            if (out.size() != 1) out.resize(1);
            out(0) = in(0);
        };
        Vector b(1, 0.0);
        Vector x(1, 1.0);
        const auto result = solve_gmres(op, b, x, 1, 4, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        EXPECT_NEAR(x(0), 0.0, 1e-12);
        EXPECT_TRUE(result.residual < 1e-12);
        EXPECT_TRUE(std::isfinite(result.residual_relative));
        EXPECT_TRUE(result.residual_relative < 1e-12);
    });

    run_case("gmres_arnoldi_uses_owned_range_for_projections", []() {
        LinearOperator op;
        op.size = 2;
        op.apply = [](const Vector& in, Vector& out) {
            if (out.size() != 2) out.resize(2);
            out(0) = in(0);
            out(1) = in(1);
        };
        Vector b(2); b(0) = 1.0; b(1) = 1.0;
        Vector x(2, 0.0);
        KrylovControls controls;
        controls.reduction.local_begin = 0;
        controls.reduction.local_end = 1; // second entry models a ghost copy
        const auto result = solve_gmres(op, b, x, 2, 4, 1e-12,
                                        nullptr, controls);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0), 1.0, 1e-12);
        EXPECT_TRUE(result.residual_relative < 1e-12);
    });

    run_case("gmres_sparse_mixed_precision_rejects_false_convergence", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0 + 1e-8);
        A.finalize();
        Vector b(1, 1.0);
        PrecisionPolicy precision;
        precision.enabled = true;
        precision.operator_precision = SolverPrecision::FP32;
        precision.reduction_precision = SolverPrecision::FP64;

        Vector x(1, 0.0);
        const auto result = solve_gmres(
            A, b, x, 1, 4, 1e-12, nullptr, KrylovControls{}, precision);
        const auto Ax64 = A.matvec(x);
        const double independent_residual = std::abs(b(0) - Ax64[0]);
        EXPECT_TRUE(result.status != SolverStatus::CONVERGED);
        EXPECT_TRUE(independent_residual > 1e-12);
        EXPECT_NEAR(result.residual, independent_residual, 1e-14);
    });

    run_case("fgmres_sparse_mixed_precision_rejects_false_convergence", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0 + 1e-8);
        A.finalize();
        Vector b(1, 1.0);
        PrecisionPolicy precision;
        precision.enabled = true;
        precision.operator_precision = SolverPrecision::FP32;
        precision.reduction_precision = SolverPrecision::FP64;

        Vector x(1, 0.0);
        const auto result = solve_fgmres(
            A, b, x, 1, 4, 1e-12, nullptr, KrylovControls{}, precision);
        const auto Ax64 = A.matvec(x);
        const double independent_residual = std::abs(b(0) - Ax64[0]);
        EXPECT_TRUE(result.status != SolverStatus::CONVERGED);
        EXPECT_TRUE(independent_residual > 1e-12);
        EXPECT_NEAR(result.residual, independent_residual, 1e-14);
    });

    run_case("gmres_nonzero_rhs_uses_rhs_relative_tolerance", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.finalize();
        Vector b(2, 1.0);
        Vector x(2, 1e6);
        const auto result = solve_gmres(A, b, x, 1, 1, 1e-8);
        const auto Ax = A.matvec(x);
        double r2 = 0.0, b2 = 0.0;
        for (std::size_t i = 0; i < 2; ++i) {
            const double r = b(i) - Ax[i];
            r2 += r * r;
            b2 += b(i) * b(i);
        }
        const double independent_relative = std::sqrt(r2 / b2);
        EXPECT_TRUE(result.status != SolverStatus::CONVERGED);
        EXPECT_TRUE(independent_relative > 1e-8);
        EXPECT_NEAR(result.residual_relative, independent_relative, 1e-12);
    });

    run_case("fgmres_nonzero_rhs_uses_rhs_relative_tolerance", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.finalize();
        Vector b(2, 1.0);
        Vector x(2, 1e6);
        const auto result = solve_fgmres(A, b, x, 1, 1, 1e-8);
        const auto Ax = A.matvec(x);
        double r2 = 0.0, b2 = 0.0;
        for (std::size_t i = 0; i < 2; ++i) {
            const double r = b(i) - Ax[i];
            r2 += r * r;
            b2 += b(i) * b(i);
        }
        const double independent_relative = std::sqrt(r2 / b2);
        EXPECT_TRUE(result.status != SolverStatus::CONVERGED);
        EXPECT_TRUE(independent_relative > 1e-8);
        EXPECT_NEAR(result.residual_relative, independent_relative, 1e-12);
    });

    run_case("gmres_rejects_nonfinite_cycle_residual", []() {
        int applications = 0;
        LinearOperator op;
        op.size = 1;
        op.apply = [&applications](const Vector& in, Vector& out) {
            if (out.size() != 1) out.resize(1);
            ++applications;
            out(0) = applications == 2
                ? std::numeric_limits<double>::quiet_NaN() : in(0);
        };
        Vector b(1, 1.0), x(1, 0.0);
        const auto result = solve_gmres(op, b, x, 1, 4, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_NEAR(result.residual, 1.0, 1e-14);
        EXPECT_NEAR(x(0), 0.0, 1e-14);
    });

    run_case("gmres_rejects_nonfinite_arnoldi_action", []() {
        int applications = 0;
        LinearOperator op;
        op.size = 1;
        op.apply = [&applications](const Vector& in, Vector& out) {
            if (out.size() != 1) out.resize(1);
            ++applications;
            out(0) = applications == 3
                ? std::numeric_limits<double>::quiet_NaN() : in(0);
        };
        Vector b(1, 1.0), x(1, 0.0);
        const auto result = solve_gmres(op, b, x, 1, 4, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_NEAR(result.residual, 1.0, 1e-14);
        EXPECT_NEAR(x(0), 0.0, 1e-14);
    });

    run_case("gmres_adaptive_restart_never_exceeds_operator_dimension", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.finalize();
        LinearOperator op;
        op.size = 2;
        op.apply = [&A](const Vector& in, Vector& out) {
            const auto values = A.matvec(in);
            if (out.size() != values.size()) out.resize(values.size());
            for (std::size_t i = 0; i < values.size(); ++i) out(i) = values[i];
        };
        Vector b(2, 1.0), x(2, 0.0);
        KrylovControls controls;
        controls.restart_min = 1;
        controls.restart_max = 512;
        controls.adaptive_restart = true;
        GmresWorkspace workspace;
        // max_iter counts total Arnoldi iterations, not restart cycles:
        // one step in the first cycle at restart=1, then two steps in the
        // second cycle at the dimension-clamped restart=2.
        const auto result = solve_gmres(op, b, x, 1, 8, 1e-14, nullptr, controls, &workspace);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(workspace.n == 2);
        EXPECT_TRUE(workspace.restart <= 2);
    });

    run_case("gmres_rejects_invalid_owned_entry_range", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();
        Vector b(2, 1.0), x(2, 0.0);
        KrylovControls controls;
        controls.reduction.local_begin = 2;
        controls.reduction.local_end = 1;
        const auto result = solve_gmres(A, b, x, 1, 4, 1e-12, nullptr, controls);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_rejects_owned_range_past_operator_dimension", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();
        Vector b(2, 1.0), x(2, 0.0);
        KrylovControls controls;
        controls.reduction.local_begin = 3;
        controls.reduction.local_end = 4;
        const auto result = solve_gmres(A, b, x, 1, 4, 1e-12, nullptr, controls);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_rejects_inverted_restart_bounds_before_clamp", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();
        Vector b(2, 1.0);
        Vector x(2, 0.0);
        KrylovControls controls;
        controls.restart_min = 8;
        controls.restart_max = 4;
        const auto result = solve_gmres(A, b, x, 4, 10, 1e-12, nullptr, controls);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_rejects_nonfinite_restart_policy", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0);
        A.finalize();
        Vector b(1, 1.0), x(1, 0.0);
        KrylovControls controls;
        controls.target_reduction_per_cycle = std::numeric_limits<double>::quiet_NaN();
        const auto result = solve_gmres(A, b, x, 1, 4, 1e-12, nullptr, controls);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_rejects_nonpositive_tolerance", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0);
        A.finalize();
        Vector b(1, 1.0), x(1, 0.0);
        const auto zero_tol = solve_gmres(A, b, x, 1, 4, 0.0);
        EXPECT_TRUE(zero_tol.status == SolverStatus::NOT_APPLICABLE);
        const auto negative_tol = solve_gmres(A, b, x, 1, 4, -1e-12);
        EXPECT_TRUE(negative_tol.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("gmres_rejects_zero_iteration_budget", []() {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 1.0);
        A.finalize();
        Vector b(1, 1.0), x(1, 0.0);
        const auto result = solve_gmres(A, b, x, 1, 0, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    return run_all();
}
