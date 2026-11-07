// Qualification of CFDX-native algorithms inspired by HYPRE/PETSc concepts.
// This test does not link to either upstream library. It checks the physical
// CSR true residual, manufactured solution, and line/AMG setup contracts.
#include "cfdx/core/linalg/algebraic_line_preconditioner.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix laplacian(std::size_t nx, std::size_t ny,
                       double ratio, bool strong_x) {
    SparseMatrix A(nx * ny, nx * ny);
    const double kx = strong_x ? ratio : 1.0;
    const double ky = strong_x ? 1.0 : ratio;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t k = j * nx + i;
            double d = 1.0;
            if (i) { A.push_back(k, k - 1, -kx); d += kx; }
            if (i + 1 < nx) { A.push_back(k, k + 1, -kx); d += kx; }
            if (j) { A.push_back(k, k - nx, -ky); d += ky; }
            if (j + 1 < ny) { A.push_back(k, k + nx, -ky); d += ky; }
            A.push_back(k, k, d);
        }
    }
    A.finalize();
    return A;
}

// Symmetric graph diffusion predominantly along northeast diagonal edges.
// Each edge contributes k (u_i-u_j)^2 and the +I reaction anchors the
// operator. This is a rotated, non-axis-aligned anisotropy stress case.
SparseMatrix diagonal_laplacian(std::size_t nx, std::size_t ny, double ratio) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    std::vector<double> diagonal(n, 1.0);
    const auto edge = [&](std::size_t i, std::size_t j, double weight) {
        A.push_back(i, j, -weight);
        A.push_back(j, i, -weight);
        diagonal[i] += weight;
        diagonal[j] += weight;
    };
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t k = j * nx + i;
            if (i + 1 < nx) edge(k, k + 1, 1.0);
            if (j + 1 < ny) edge(k, k + nx, 1.0);
            if (i + 1 < nx && j + 1 < ny) edge(k, k + nx + 1, ratio);
        }
    for (std::size_t i = 0; i < n; ++i) A.push_back(i, i, diagonal[i]);
    A.finalize();
    return A;
}

Vector exact_solution(std::size_t nx, std::size_t ny) {
    Vector x(nx * ny, 0.0);
    constexpr double pi = 3.14159265358979323846264;
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            x(j * nx + i) = 0.25 +
                std::sin(pi * (static_cast<double>(i) + 0.5) / static_cast<double>(nx)) *
                std::sin(2.0 * pi * (static_cast<double>(j) + 0.5) / static_cast<double>(ny));
    return x;
}

Vector rhs_of(const SparseMatrix& A, const Vector& x) {
    const auto ax = A.matvec(x);
    Vector rhs(ax.size(), 0.0);
    for (std::size_t i = 0; i < ax.size(); ++i) rhs(i) = ax[i];
    return rhs;
}

double relative_error(const Vector& x, const Vector& reference) {
    double d2 = 0.0, r2 = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double e = x(i) - reference(i);
        d2 += e * e;
        r2 += reference(i) * reference(i);
    }
    return std::sqrt(d2 / std::max(r2, 1e-300));
}

void verify(std::size_t nx, std::size_t ny, double ratio, bool x_direction,
            KrylovModel krylov, PreconditionerModel pc) {
    const auto A = laplacian(nx, ny, ratio, x_direction);
    const auto exact = exact_solution(nx, ny);
    const auto rhs = rhs_of(A, exact);
    Vector out(A.n_rows(), 0.0);
    LinearSolverRequest request;
    request.krylov = krylov;
    request.preconditioner = pc;
    request.gmres_restart = 35;
    const auto before = std::chrono::steady_clock::now();
    const auto report = solve_linear_system(A, rhs, out, LinearProblemKind::PressurePoisson,
                                            request, 250, 1e-9);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - before).count();
    const double err = relative_error(out, exact);
    if (report.result.status != SolverStatus::CONVERGED) {
        auto diagnostic_pc = make_scalar_preconditioner(pc);
        if (auto* amg = dynamic_cast<NativeAMGPreconditioner*>(diagnostic_pc.get())) {
            const bool ok = amg->setup(A);
            std::cerr << "AMG setup retry=" << ok << " reason=" << amg->last_error() << '\n';
        }
    }
    std::cout << std::setprecision(12) << "native_port_csv,"
              << nx << ',' << ny << ',' << ratio << ',' << (x_direction ? 'x' : 'y')
              << ',' << to_string(krylov) << ',' << to_string(pc)
              << ',' << static_cast<int>(report.result.status)
              << ',' << report.result.iterations
              << ',' << report.physical_residual_relative
              << ',' << err << ',' << elapsed << '\n';
    EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
    EXPECT_TRUE(std::isfinite(report.physical_residual_relative) &&
                report.physical_residual_relative <= 1e-9);
    EXPECT_TRUE(std::isfinite(err) && err <= 1e-3);
}
} // namespace

int main(int argc, char** argv) {
    const bool benchmark = argc > 1 &&
        (std::string(argv[1]) == "--benchmark" ||
         std::string(argv[1]) == "--ablation");
    const std::size_t nx = benchmark && argc > 2 ? std::stoul(argv[2]) : 32;
    const std::size_t ny = benchmark && argc > 3 ? std::stoul(argv[3]) : 24;
    if (nx < 4 || ny < 4 || nx * ny > 65536) return 2;

    if (argc > 1 && std::string(argv[1]) == "--ablation") {
        // Evaluate the independently selectable PETSc-inspired filtering and
        // line-block relaxation before claiming which component saves time.
        for (const double ratio : {1e3, 1e6}) {
            const auto A = laplacian(nx, ny, ratio, true);
            const auto exact = exact_solution(nx, ny);
            const auto rhs = rhs_of(A, exact);
            for (const bool lines : {false, true}) {
                for (const double filter : {0.0, 0.05}) {
                    NativeAMGPreconditioner pc(
                        NativeAMGMethod::SmoothedAggregation, false, false,
                        lines, filter);
                    Vector x(A.n_rows(), 0.0);
                    const auto start = std::chrono::steady_clock::now();
                    const bool prepared = pc.setup(A);
                    const auto setup_end = std::chrono::steady_clock::now();
                    SolverResult result;
                    result.status = SolverStatus::NOT_APPLICABLE;
                    if (prepared)
                        result = solve_fgmres(A, rhs, x, 35, 250, 1e-9, &pc);
                    const auto solve_end = std::chrono::steady_clock::now();
                    const double elapsed = std::chrono::duration<double>(
                        solve_end - start).count();
                    const double setup_seconds = std::chrono::duration<double>(
                        setup_end - start).count();
                    const double solve_seconds = std::chrono::duration<double>(
                        solve_end - setup_end).count();
                    const double residual = physical_relative_residual(A, x, rhs);
                    const double error = relative_error(x, exact);
                    std::cout << std::setprecision(12) << "ablation_csv,"
                              << nx * ny << ',' << ratio << ',' << lines << ','
                              << filter << ',' << prepared << ','
                              << static_cast<int>(result.status) << ','
                              << result.iterations << ',' << residual << ','
                              << error << ',' << elapsed << ','
                              << setup_seconds << ',' << solve_seconds << ','
                              << pc.first_prolongation_nnz() << ','
                              << pc.hierarchy_level_sizes().size() << ','
                              << pc.finest_line_coverage() << ','
                              << pc.line_factor_fallback_levels() << '\n';
                    if (!prepared || result.status != SolverStatus::CONVERGED ||
                        residual > 1e-9 || error > 2e-3) return 1;
                }
            }
        }
        return 0;
    }

    run_case("csr_line_graph_is_anisotropic_not_geometric", [&] {
        for (const bool strong_x : {false, true}) {
            const auto A = laplacian(nx, ny, 1e6, strong_x);
            AlgebraicLinePreconditioner lines;
            EXPECT_TRUE(lines.setup(A));
            EXPECT_TRUE(lines.covered_rows() >= A.n_rows() / 2);
            EXPECT_TRUE(lines.longest_line() == (strong_x ? nx : ny));
            Vector rhs(A.n_rows(), 1.0), z(A.n_rows(), 0.0);
            EXPECT_TRUE(lines.apply(rhs, z));
            const auto saved = z;
            EXPECT_TRUE(lines.apply(rhs, z));
            EXPECT_TRUE(relative_error(z, saved) < 1e-14);
            // Thomas forward/backward substitutions must work without a
            // temporary heap vector and when RHS/output alias.
            Vector in_place = rhs;
            EXPECT_TRUE(lines.apply(in_place, in_place));
            EXPECT_TRUE(relative_error(in_place, saved) < 1e-14);
            Vector varying(rhs.size(), 0.0), varied(rhs.size(), 0.0);
            for (std::size_t i = 0; i < varying.size(); ++i)
                varying(i) = 0.125 + static_cast<double>(i % 23);
            EXPECT_TRUE(lines.apply(varying, varied));
            Vector varied_alias = varying;
            EXPECT_TRUE(lines.apply(varied_alias, varied_alias));
            EXPECT_TRUE(relative_error(varied_alias, varied) < 1e-14);
            EXPECT_TRUE(lines.update_values(A));
            EXPECT_TRUE(lines.apply(rhs, z));
        }
    });

    run_case("anisotropy_portable_krylov_and_amg", [&] {
        for (const bool strong_x : {false, true}) {
            for (const double ratio : {1e3, 1e6}) {
                for (const auto pc : {
                       PreconditionerModel::AnisotropicLine,
                       PreconditionerModel::AnisotropicBoomerAMG,
                       PreconditionerModel::AnisotropicSmoothedAggregationAMG}) {
                    verify(nx, ny, ratio, strong_x, KrylovModel::FGMRES, pc);
                }
            }
        }
    });
    run_case("petsc_inspired_interpolation_preserves_constant_mode", [&] {
        const auto A = laplacian(nx, ny, 1e6, false);
        NativeAnisotropicSmoothedAggregationAMGPreconditioner amg;
        EXPECT_TRUE(amg.setup(A));
        EXPECT_TRUE(amg.finest_line_coverage() > A.n_rows() / 2);
        EXPECT_TRUE(amg.finest_longest_line() == ny);
        std::cout << "line_factor_fallback_levels="
                  << amg.line_factor_fallback_levels() << '\n';
        EXPECT_TRUE(amg.prolongation_row_sum_min() > 1.0 - 1e-9);
        EXPECT_TRUE(amg.prolongation_row_sum_max() < 1.0 + 1e-9);
    });
    run_case("rotated_direction_strong_graph_preserves_physical_solution", [&] {
        const auto A = diagonal_laplacian(nx, ny, 1e4);
        const auto exact = exact_solution(nx, ny);
        const auto rhs = rhs_of(A, exact);
        AlgebraicLinePreconditioner line_pc;
        EXPECT_TRUE(line_pc.setup(A));
        EXPECT_TRUE(line_pc.covered_rows() > A.n_rows() / 2);
        EXPECT_TRUE(line_pc.longest_line() >= std::min(nx, ny) / 2);
        for (auto kind : {PreconditionerModel::AnisotropicBoomerAMG,
                          PreconditionerModel::AnisotropicSmoothedAggregationAMG,
                          PreconditionerModel::AnisotropicLine}) {
            Vector x(A.n_rows(), 0.0);
            LinearSolverRequest request;
            request.krylov = KrylovModel::FGMRES;
            request.preconditioner = kind;
            request.gmres_restart = 40;
            const auto report = solve_linear_system(A, rhs, x,
                LinearProblemKind::PressurePoisson, request, 300, 1e-9);
            const double error = relative_error(x, exact);
            std::cout << "rotated_csv," << nx * ny << ',' << to_string(kind)
                      << ',' << report.result.iterations << ','
                      << report.physical_residual_relative << ',' << error << '\n';
            EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
            EXPECT_TRUE(std::isfinite(report.physical_residual_relative) &&
                        report.physical_residual_relative <= 1e-9);
            EXPECT_TRUE(std::isfinite(error) && error <= 1e-4);
        }
    });
    run_case("non_anisotropic_line_fallback_and_numeric_reuse", [&] {
        const auto A = laplacian(nx, ny, 1.0, true);
        AlgebraicLinePreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.covered_rows() == 0);
        EXPECT_TRUE(pc.longest_line() == 1);
        const auto B = laplacian(nx, ny, 100.0, true);
        // Coefficients changed the strong-direction graph while keeping the
        // same CSR. This update is rejected without poisoning the old PC.
        EXPECT_TRUE(!pc.update_values(B));
        Vector rhs(B.n_rows(), 1.0), z(B.n_rows(), 0.0);
        EXPECT_TRUE(pc.apply(rhs, z));
        EXPECT_TRUE(z.is_valid());
        EXPECT_TRUE(pc.covered_rows() == 0);
        EXPECT_TRUE(pc.setup(B));
        EXPECT_TRUE(pc.covered_rows() >= B.n_rows() / 2);
    });
    run_case("automatic_krylov_pairing_and_spd_safety", [&] {
        const auto A = laplacian(nx, ny, 1.0, true);
        const auto features = measure_matrix_characteristics(A);
        LinearSolverRequest request;
        request.preconditioner = PreconditionerModel::AnisotropicBoomerAMG;
        const auto plan = select_linear_solver(
            LinearProblemKind::PressurePoisson, features, request);
        EXPECT_TRUE(plan.krylov == KrylovModel::FGMRES);
        request.krylov = KrylovModel::CG;
        EXPECT_THROW(select_linear_solver(
            LinearProblemKind::PressurePoisson, features, request),
            std::invalid_argument);
        request.krylov = KrylovModel::Auto;
        request.preconditioner = PreconditionerModel::AnisotropicSmoothedAggregationAMG;
        EXPECT_TRUE(select_linear_solver(LinearProblemKind::PressurePoisson,
                       features, request).krylov == KrylovModel::FGMRES);
    });
    run_case("numeric_refresh_failure_does_not_mutate_line_factorization", [&] {
        const auto A = laplacian(nx, ny, 1e3, true);
        AlgebraicLinePreconditioner lines;
        EXPECT_TRUE(lines.setup(A));
        Vector rhs(A.n_rows(), 1.0), original(A.n_rows(), 0.0), after(A.n_rows(), 0.0);
        EXPECT_TRUE(lines.apply(rhs, original));
        SparseMatrix unusable(A.n_rows(), A.n_cols());
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i)
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                unusable.push_back(i, col[k], col[k] == i ? 0.0 : val[k]);
        unusable.finalize();
        EXPECT_TRUE(!lines.update_values(unusable));
        EXPECT_TRUE(lines.apply(rhs, after));
        EXPECT_TRUE(relative_error(after, original) < 1e-14);
    });
    run_case("native_amg_numeric_update_is_still_available", [&] {
        const auto A = laplacian(nx, ny, 1e3, true);
        const auto B = laplacian(nx, ny, 1e4, true);
        NativeAnisotropicBoomerAMGPreconditioner amg;
        EXPECT_TRUE(amg.setup(A));
        EXPECT_TRUE(amg.update_values(B));
        EXPECT_TRUE(amg.hierarchy_builds() == 1);
        EXPECT_TRUE(amg.numeric_updates() == 1);
        EXPECT_TRUE(amg.finest_line_coverage() > A.n_rows() / 2);
        const auto exact = exact_solution(nx, ny);
        const auto rhs = rhs_of(B, exact);
        Vector solved(B.n_rows(), 0.0);
        const auto result = solve_fgmres(B, rhs, solved, 35, 250, 1e-9, &amg);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(physical_relative_residual(B, solved, rhs) <= 1e-9);
        EXPECT_TRUE(relative_error(solved, exact) <= 1e-3);
        std::cout << "amg_numeric_refresh_fallbacks="
                  << amg.numeric_line_fallbacks() << '\n';
    });
    return run_all();
}