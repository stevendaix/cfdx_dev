// Anisotropic pressure/diffusion manufactured-solution qualification.
// Standalone: --benchmark [nx ny] runs a larger timing matrix; no timing
// threshold in CTest (wall-clock timing is hardware-dependent).
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/linear_solver_context.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix make_stretched_diffusion(std::size_t nx, std::size_t ny,
                                      double contrast, bool strong_x) {
    // FV-style axis-aligned diffusion on a stretched Cartesian grid. A > 0
    // (reaction) anchors pressure and prevents an unqualified null space.
    const double ax = strong_x ? contrast : 1.0;
    const double ay = strong_x ? 1.0 : contrast;
    const double reaction = 1.0;
    SparseMatrix A(nx * ny, nx * ny);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t p = j * nx + i;
            double diagonal = reaction;
            auto add_face = [&](std::size_t q, double conductivity) {
                A.push_back(p, q, -conductivity);
                diagonal += conductivity;
            };
            if (i > 0) add_face(p - 1, ax);
            if (i + 1 < nx) add_face(p + 1, ax);
            if (j > 0) add_face(p - nx, ay);
            if (j + 1 < ny) add_face(p + nx, ay);
            A.push_back(p, p, diagonal);
        }
    }
    A.finalize();
    return A;
}

Vector manufactured_solution(std::size_t nx, std::size_t ny) {
    Vector x(nx * ny, 0.0);
    constexpr double pi = 3.141592653589793;
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const double xx = (static_cast<double>(i) + 0.5) / static_cast<double>(nx);
            const double yy = (static_cast<double>(j) + 0.5) / static_cast<double>(ny);
            x(j * nx + i) = 0.25 + std::sin(pi * xx) * std::sin(2 * pi * yy);
        }
    return x;
}

Vector rhs_from_exact(const SparseMatrix& A, const Vector& exact) {
    const auto y = A.matvec(exact);
    Vector rhs(y.size(), 0.0);
    for (std::size_t i = 0; i < y.size(); ++i) rhs(i) = y[i];
    return rhs;
}

double error_norm(const Vector& x, const Vector& exact) {
    double err2 = 0.0, ref2 = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        err2 += (x(i) - exact(i)) * (x(i) - exact(i));
        ref2 += exact(i) * exact(i);
    }
    return std::sqrt(err2 / std::max(ref2, 1e-300));
}

struct Outcome {
    SolverStatus status;
    double true_residual;
    double error;
    std::size_t iterations;
    double elapsed;
};

Outcome measure(std::size_t nx, std::size_t ny, double ratio, bool strong_x,
                KrylovModel krylov, PreconditionerModel preconditioner,
                LinearScalingModel scaling = LinearScalingModel::None) {
    const auto A = make_stretched_diffusion(nx, ny, ratio, strong_x);
    const auto exact = manufactured_solution(nx, ny);
    const auto rhs = rhs_from_exact(A, exact);
    Vector solution(A.n_rows(), 0.0);
    LinearSolverRequest request;
    request.krylov = krylov;
    request.preconditioner = preconditioner;
    request.scaling = scaling;
    request.gmres_restart = 35;
    const auto begin = std::chrono::steady_clock::now();
    const auto report = solve_linear_system(
        A, rhs, solution, LinearProblemKind::PressurePoisson,
        request, 350, 1e-9);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    const auto out = Outcome{report.result.status, report.physical_residual_relative,
                             error_norm(solution, exact),
                             report.result.iterations, elapsed};
    std::cout << std::setprecision(12)
              << "aniso_csv," << (strong_x ? "x" : "y")
              << ',' << ratio << ',' << A.n_rows() << ',' << to_string(krylov)
              << ',' << to_string(preconditioner) << ',' << to_string(scaling)
              << ',' << static_cast<int>(out.status) << ',' << out.iterations
              << ',' << out.true_residual << ',' << out.error
              << ',' << out.elapsed << '\n';
    return out;
}

void assert_qualified(const Outcome& o) {
    EXPECT_TRUE(o.status == SolverStatus::CONVERGED);
    EXPECT_TRUE(std::isfinite(o.true_residual) && o.true_residual <= 1e-9);
    EXPECT_TRUE(std::isfinite(o.error) && o.error <= 2e-3);
}
} // namespace

int main(int argc, char** argv) {
    // A shifted isotropic operator can make strength coarsening stop at
    // thousands of unknowns. This regression must never use a dense solve on
    // that terminal level. A separate CTest timeout catches cubic fallbacks.
    if (argc > 1 && std::string(argv[1]) == "--coarse-stall") {
        run_case("sparse_terminal_level_8192_unknowns", [&] {
            constexpr std::size_t sx = 128, sy = 64;
            const auto A = make_stretched_diffusion(sx, sy, 1.0, true);
            NativeSmoothedAggregationAMGPreconditioner amg;
            EXPECT_TRUE(amg.setup(A));
            const auto levels = amg.hierarchy_level_sizes();
            EXPECT_TRUE(levels.size() >= 2);
            EXPECT_TRUE(levels.back() > 16);
            assert_qualified(measure(sx, sy, 1.0, true,
                                     KrylovModel::FGMRES,
                                     PreconditionerModel::SmoothedAggregationAMG));
            assert_qualified(measure(sx, sy, 1.0, true,
                                     KrylovModel::FGMRES,
                                     PreconditionerModel::NativeAMG));
        });
        return run_all();
    }
    const bool large = argc > 1 && std::string(argv[1]) == "--benchmark";
    const std::size_t nx = large && argc > 2 ? std::stoul(argv[2]) : 32;
    const std::size_t ny = large && argc > 3 ? std::stoul(argv[3]) : 24;
    if (nx < 4 || ny < 4 || nx * ny > 262144)
        throw std::invalid_argument("benchmark dimensions outside qualification bounds");

    run_case("strongly_directional_auto_dispatch_uses_flexible_amg", [&] {
        const auto A = make_stretched_diffusion(nx, ny, 1e6, true);
        const auto c = measure_matrix_characteristics(A);
        EXPECT_TRUE(c.numerically_symmetric);
        EXPECT_TRUE(c.diagonally_dominant);
        EXPECT_TRUE(c.positive_diagonal);
        EXPECT_TRUE(c.strongly_directional_coupling);
        const auto plan = select_linear_solver(LinearProblemKind::PressurePoisson, c);
        EXPECT_TRUE(plan.krylov == KrylovModel::FGMRES);
        EXPECT_TRUE(plan.preconditioner == PreconditionerModel::NativeAMG);
        assert_qualified(measure(nx, ny, 1e6, true,
                                 KrylovModel::Auto, PreconditionerModel::Auto));
    });
    run_case("anisotropic_native_amg_and_smoothed_aggregation", [&] {
        for (const bool strong_x : {true, false})
            for (const double ratio : {1.0, 1e3, 1e6})
                for (const auto pc : {PreconditionerModel::NativeAMG,
                                      PreconditionerModel::SmoothedAggregationAMG})
                    assert_qualified(measure(nx, ny, ratio, strong_x,
                                             KrylovModel::FGMRES, pc));
    });
    run_case("cg_amg_spd_and_ilu0_gmres", [&] {
        // CG+AMG is verified on the certified SPD hierarchy, not assumed valid
        // for every possible anisotropy or nonsymmetric factorization.
        assert_qualified(measure(nx, ny, 1e3, true,
                                 KrylovModel::CG, PreconditionerModel::NativeAMG));
        assert_qualified(measure(nx, ny, 1e3, false,
                                 KrylovModel::CG, PreconditionerModel::SmoothedAggregationAMG));
        assert_qualified(measure(nx, ny, 1e3, true,
                                 KrylovModel::GMRES, PreconditionerModel::ILU0));
        assert_qualified(measure(nx, ny, 1e3, false,
                                 KrylovModel::FGMRES, PreconditionerModel::ILU0));
        assert_qualified(measure(nx, ny, 1e3, true,
                                 KrylovModel::GMRES, PreconditionerModel::NativeAMG));
        assert_qualified(measure(nx, ny, 1e3, false,
                                 KrylovModel::GMRES, PreconditionerModel::SmoothedAggregationAMG));
    });
    run_case("scaling_and_invalid_krylov_pairs_are_fail_closed", [&] {
        const auto A = make_stretched_diffusion(nx, ny, 1e4, true);
        const auto c = measure_matrix_characteristics(A);
        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::ILU0;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson, c, request),
                     std::invalid_argument);
        request.preconditioner = PreconditionerModel::NativeAMG;
        request.scaling = LinearScalingModel::Row;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson, c, request),
                     std::invalid_argument);
        assert_qualified(measure(nx, ny, 1e3, false,
                                 KrylovModel::FGMRES, PreconditionerModel::NativeAMG,
                                 LinearScalingModel::SymmetricDiagonal));
    });
    return run_all();
}