#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_tridiagonal_nonsymmetric(std::size_t n) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) A.push_back(i, i - 1, -1.2);
        A.push_back(i, i, 4.0);
        if (i + 1 < n) A.push_back(i, i + 1, -0.8);
    }
    A.finalize();
    return A;
}

static SparseMatrix make_fvm_diffusion_2d(std::size_t nx, std::size_t ny,
                                          double kx, double ky, double reaction = 1.0) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diag = reaction;
            if (x > 0) { A.push_back(i, i - 1, -kx); diag += kx; }
            if (x + 1 < nx) { A.push_back(i, i + 1, -kx); diag += kx; }
            if (y > 0) { A.push_back(i, i - nx, -ky); diag += ky; }
            if (y + 1 < ny) { A.push_back(i, i + nx, -ky); diag += ky; }
            A.push_back(i, i, diag);
        }
    }
    A.finalize();
    return A;
}

// Cell-centred convection-diffusion stencil. Positive east/west bias makes
// the operator nonsymmetric while retaining diagonal dominance typical of
// a stable upwind finite-volume discretisation.
static SparseMatrix make_convection_diffusion_2d(std::size_t nx, std::size_t ny,
                                                 double diffusion, double convection) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diag = 1.0 + 2.0 * diffusion + convection;
            if (x > 0) {
                A.push_back(i, i - 1, -diffusion - convection);
            }
            if (x + 1 < nx) {
                A.push_back(i, i + 1, -diffusion);
            }
            if (y > 0) {
                A.push_back(i, i - nx, -diffusion);
            }
            if (y + 1 < ny) {
                A.push_back(i, i + nx, -diffusion);
            }
            A.push_back(i, i, diag);
        }
    }
    A.finalize();
    return A;
}

static Vector make_rhs_for_ones(const SparseMatrix& A) {
    Vector x(A.n_rows(), 1.0);
    const auto values = A.matvec(x);
    Vector b(A.n_rows(), 0.0);
    for (std::size_t i = 0; i < b.size(); ++i) b(i) = values[i];
    return b;
}

static double true_residual(const SparseMatrix& A, const Vector& x, const Vector& b) {
    const auto Ax = A.matvec(x);
    double sum = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        const double r = Ax[i] - b(i);
        sum += r * r;
    }
    return std::sqrt(sum);
}

static double relative_true_residual(const SparseMatrix& A, const Vector& x,
                                     const Vector& b) {
    const double denom = std::max(std::sqrt(std::inner_product(
        b.data(), b.data() + b.size(), b.data(), 0.0)), 1e-30);
    return true_residual(A, x, b) / denom;
}

int main() {
    run_case("gmres_preconditioner_e2e_iteration_and_true_residual", [] {
        const auto A = make_tridiagonal_nonsymmetric(20);
        const auto b = make_rhs_for_ones(A);

        Vector x_plain(20, 0.0);
        const auto plain = solve_gmres(A, b, x_plain, 5, 200, 1e-12);
        EXPECT_TRUE(plain.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(A, x_plain, b) < 1e-10);

        JacobiPreconditioner jacobi;
        Vector x_jacobi(20, 0.0);
        const auto jacobi_result = solve_gmres(A, b, x_jacobi, 5, 200, 1e-12, &jacobi);
        EXPECT_TRUE(jacobi_result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(A, x_jacobi, b) < 1e-10);

        GaussSeidelPreconditioner gs;
        Vector x_gs(20, 0.0);
        const auto gs_result = solve_gmres(A, b, x_gs, 5, 200, 1e-12, &gs);
        EXPECT_TRUE(gs_result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(A, x_gs, b) < 1e-10);

        ILU0Preconditioner ilu;
        Vector x_ilu(20, 0.0);
        const auto ilu_result = solve_gmres(A, b, x_ilu, 5, 200, 1e-12, &ilu);
        EXPECT_TRUE(ilu_result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(true_residual(A, x_ilu, b) < 1e-10);
        EXPECT_TRUE(ilu_result.iterations < plain.iterations);

        std::cout << "GMRES E2E iterations: plain=" << plain.iterations
                  << " Jacobi=" << jacobi_result.iterations
                  << " GS=" << gs_result.iterations
                  << " ILU0=" << ilu_result.iterations << "\n";
    });

    run_case("bicgstab_ilu0_e2e_true_residual", [] {
        const auto A = make_tridiagonal_nonsymmetric(20);
        const auto b = make_rhs_for_ones(A);

        Vector x_plain(20, 0.0);
        const auto plain = solve_bicgstab(A, b, x_plain, 200, 1e-12);
        EXPECT_TRUE(plain.status == SolverStatus::CONVERGED);
        const double plain_residual = true_residual(A, x_plain, b);
        EXPECT_TRUE(plain_residual < 1e-10);
        EXPECT_TRUE(std::abs(plain.residual - plain_residual) < 1e-10);

        ILU0Preconditioner ilu;
        Vector x_ilu(20, 0.0);
        const auto ilu_result = solve_bicgstab(A, b, x_ilu, 200, 1e-12, &ilu);
        EXPECT_TRUE(ilu_result.status == SolverStatus::CONVERGED);
        const double ilu_residual = true_residual(A, x_ilu, b);
        EXPECT_TRUE(ilu_residual < 1e-10);
        EXPECT_TRUE(std::abs(ilu_result.residual - ilu_residual) < 1e-10);
        EXPECT_TRUE(ilu_result.iterations <= plain.iterations);

        std::cout << "BiCGStab E2E iterations: plain=" << plain.iterations
                  << " ILU0=" << ilu_result.iterations << "\n";
    });

    run_case("gmres_aborts_before_iteration_on_preconditioner_setup_failure", [] {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 0.0);
        A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(2, 2, 1.0);
        A.finalize();

        Vector b(3, 1.0), x(3, 0.0);
        ILU0Preconditioner ilu;
        EXPECT_TRUE(!ilu.setup(A));

        const auto result = solve_gmres(A, b, x, 5, 100, 1e-12, &ilu);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(result.iterations == 0);
        EXPECT_TRUE(true_residual(A, x, b) > 0.0);
    });

    run_case("bicgstab_aborts_before_iteration_on_preconditioner_setup_failure", [] {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 0.0);
        A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 1.0);
        A.push_back(1, 1, 2.0);
        A.push_back(2, 2, 1.0);
        A.finalize();

        Vector b(3, 1.0), x(3, 0.0);
        ILU0Preconditioner ilu;
        const auto result = solve_bicgstab(A, b, x, 100, 1e-12, &ilu);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(result.iterations == 0);
        EXPECT_TRUE(true_residual(A, x, b) > 0.0);
    });

    run_case("cg_fvm_diffusion_matrix_robustness", [] {
        const auto A = make_fvm_diffusion_2d(12, 10, 1.0, 1.0);
        const auto b = make_rhs_for_ones(A);

        Vector x(120, 0.0);
        const auto result = solve_cg(A, b, x, 500, 1e-10);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_TRUE(std::isfinite(result.residual_relative));
        EXPECT_TRUE(relative_true_residual(A, x, b) < 1e-8);
        for (std::size_t i = 0; i < x.size(); ++i) EXPECT_TRUE(std::isfinite(x(i)));
    });

    run_case("cg_strong_anisotropy_fvm_matrix_robustness", [] {
        const auto A = make_fvm_diffusion_2d(12, 12, 1.0, 1e3);
        const auto b = make_rhs_for_ones(A);

        Vector x(144, 0.0);
        JacobiPreconditioner jacobi;
        const auto result = solve_cg(A, b, x, 1000, 1e-10, PrecisionPolicy{});
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(A, x, b) < 1e-8);
        EXPECT_TRUE(jacobi.setup(A));
        Vector z(144, 0.0);
        EXPECT_TRUE(jacobi.apply(b, z));
        for (std::size_t i = 0; i < z.size(); ++i) EXPECT_TRUE(std::isfinite(z(i)));
    });

    run_case("gmres_convection_diffusion_fvm_matrix_robustness", [] {
        const auto A = make_convection_diffusion_2d(20, 10, 0.2, 1.5);
        const auto b = make_rhs_for_ones(A);

        Vector x(200, 0.0);
        ILU0Preconditioner ilu;
        const auto result = solve_gmres(A, b, x, 10, 500, 1e-10, &ilu);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        EXPECT_TRUE(relative_true_residual(A, x, b) < 1e-8);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_TRUE(std::isfinite(result.residual_relative));
    });

    run_case("bicgstab_convection_diffusion_fvm_matrix_robustness", [] {
        const auto A = make_convection_diffusion_2d(20, 10, 0.2, 1.5);
        const auto b = make_rhs_for_ones(A);

        Vector x(200, 0.0);
        ILU0Preconditioner ilu;
        const auto result = solve_bicgstab(A, b, x, 500, 1e-10, &ilu);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        EXPECT_TRUE(relative_true_residual(A, x, b) < 1e-8);
        EXPECT_TRUE(std::isfinite(result.residual));
        EXPECT_TRUE(std::isfinite(result.residual_relative));
    });

    run_case("gmres_bad_scaling_remains_finite_and_converges", [] {
        SparseMatrix A(8, 8);
        Vector x_exact(8);
        for (std::size_t i = 0; i < 8; ++i) {
            const double scale = std::pow(10.0, -4.0 + static_cast<double>(i));
            A.push_back(i, i, scale);
            if (i > 0) A.push_back(i, i - 1, -0.1 * scale);
            if (i + 1 < 8) A.push_back(i, i + 1, -0.05 * scale);
            x_exact(i) = (i % 2 == 0) ? 1.0 : -0.5;
        }
        A.finalize();

        const auto b_values = A.matvec(x_exact);
        Vector b(8);
        for (std::size_t i = 0; i < 8; ++i) b(i) = b_values[i];

        Vector x(8, 0.0);
        ILU0Preconditioner ilu;
        const auto result = solve_gmres(A, b, x, 8, 500, 1e-10, &ilu);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(relative_true_residual(A, x, b) < 1e-8);
        for (std::size_t i = 0; i < x.size(); ++i) {
            EXPECT_TRUE(std::isfinite(x(i)));
            EXPECT_NEAR(x(i), x_exact(i), 1e-6);
        }
    });

    return run_all();
}
