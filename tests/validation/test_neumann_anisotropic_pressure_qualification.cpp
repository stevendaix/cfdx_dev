// Qualification of a pressure Poisson operator with a true constant null mode.
// Nonuniform Cartesian finite volumes: the 2-point interior transmissibilities
// are assembled independently of the face-flux balance used to verify them.
// No reaction term, no pinned pressure row and no relaxation of true residuals.
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
constexpr double pi = 3.14159265358979323846264338327950288;
struct Face { std::size_t a, b; double transmissibility; };
struct VolumeMesh {
    std::size_t nx, ny;
    std::vector<double> xf, yf, xc, yc, volume;
    std::vector<Face> faces;
    SparseMatrix matrix;
    explicit VolumeMesh(std::size_t nx_, std::size_t ny_, double ratio,
                        bool strong_x, double stretch = 1.5)
        : nx(nx_), ny(ny_), xf(nx_ + 1), yf(ny_ + 1),
          xc(nx_), yc(ny_), volume(nx_ * ny_),
          matrix(nx_ * ny_, nx_ * ny_) {
        if (nx < 3 || ny < 3 || !(ratio > 0.0))
            throw std::invalid_argument("invalid pressure qualification grid");
        for (std::size_t i = 0; i <= nx; ++i)
            xf[i] = std::pow(static_cast<double>(i) / static_cast<double>(nx), stretch);
        for (std::size_t j = 0; j <= ny; ++j)
            yf[j] = std::pow(static_cast<double>(j) / static_cast<double>(ny), stretch);
        for (std::size_t i = 0; i < nx; ++i) xc[i] = (xf[i] + xf[i + 1]) * 0.5;
        for (std::size_t j = 0; j < ny; ++j) yc[j] = (yf[j] + yf[j + 1]) * 0.5;
        const double kx = strong_x ? ratio : 1.0;
        const double ky = strong_x ? 1.0 : ratio;
        std::vector<double> diag(nx * ny, 0.0);
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t k = j * nx + i;
                const double dx = xf[i + 1] - xf[i];
                const double dy = yf[j + 1] - yf[j];
                volume[k] = dx * dy;
                if (i + 1 < nx) {
                    const double t = kx * dy / (xc[i + 1] - xc[i]);
                    faces.push_back({k, k + 1, t});
                }
                if (j + 1 < ny) {
                    const double t = ky * dx / (yc[j + 1] - yc[j]);
                    faces.push_back({k, k + nx, t});
                }
            }
        for (const auto& f : faces) {
            diag[f.a] += f.transmissibility;
            diag[f.b] += f.transmissibility;
            matrix.push_back(f.a, f.b, -f.transmissibility);
            matrix.push_back(f.b, f.a, -f.transmissibility);
        }
        for (std::size_t k = 0; k < diag.size(); ++k)
            matrix.push_back(k, k, diag[k]);
        matrix.finalize();
    }

    Vector analytic_pressure() const {
        Vector pressure(nx * ny);
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i)
                pressure(j * nx + i) = std::cos(pi * xc[i]) * std::cos(2.0 * pi * yc[j]);
        NullSpaceProjector::constant(nx * ny).remove(pressure);
        return pressure;
    }

    Vector source_quadrature(double ratio, bool strong_x) const {
        Vector source(nx * ny);
        const double kx = strong_x ? ratio : 1.0;
        const double ky = strong_x ? 1.0 : ratio;
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t k = j * nx + i;
                source(k) = pi * pi * (kx + 4.0 * ky) *
                    std::cos(pi * xc[i]) * std::cos(2.0 * pi * yc[j]) * volume[k];
            }
        return source;
    }

    Vector face_balance(const Vector& pressure) const {
        Vector balance(nx * ny, 0.0);
        for (const auto& face : faces) {
            const double flux = face.transmissibility *
                (pressure(face.a) - pressure(face.b));
            balance(face.a) += flux;
            balance(face.b) -= flux;
        }
        return balance;
    }
};

Vector apply(const SparseMatrix& A, const Vector& x) {
    const auto result = A.matvec(x);
    Vector b(result.size());
    for (std::size_t i = 0; i < result.size(); ++i) b(i) = result[i];
    return b;
}

double gauge_free_error(const Vector& x, const Vector& exact) {
    Vector error = x - exact;
    NullSpaceProjector::constant(x.size()).remove(error);
    return error.norm2() / std::max(exact.norm2(), 1e-300);
}

double balance_error(const Vector& a, const Vector& b) {
    return (a - b).norm2() / std::max(b.norm2(), 1e-300);
}

struct Observation {
    std::size_t iterations{};
    double residual{}, solution_error{}, flux_error{}, seconds{};
    SolverStatus status = SolverStatus::NOT_APPLICABLE;
};

Observation run_solve(const VolumeMesh& grid, const Vector& rhs,
                      const Vector& exact, KrylovModel krylov,
                      PreconditionerModel preconditioner,
                      double ratio, bool strong_x, double tolerance = 1e-9,
                      bool qualification_required = true) {
    Vector solution(grid.nx * grid.ny, 15.0); // arbitrary initial pressure gauge
    LinearSolverRequest request;
    request.krylov = krylov;
    request.preconditioner = preconditioner;
    request.null_space = NullSpaceModel::Constant;
    request.gmres_restart = 40;
    const auto start = std::chrono::steady_clock::now();
    const auto report = solve_linear_system(grid.matrix, rhs, solution,
        LinearProblemKind::PressurePoisson, request, 350, tolerance);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    const double flux_err = balance_error(grid.face_balance(solution), rhs);
    const double err = gauge_free_error(solution, exact);
    const double gauge = NullSpaceProjector::constant(solution.size()).component_norm(solution);
    std::cout << std::setprecision(12) << "neumann_qualification_csv,"
              << grid.nx << ',' << grid.ny << ',' << ratio << ','
              << (strong_x ? "x" : "y") << ','
              << to_string(krylov) << ',' << to_string(preconditioner) << ','
              << static_cast<int>(report.result.status) << ','
              << report.result.iterations << ','
              << report.physical_residual_relative << ',' << err << ','
              << flux_err << ',' << gauge << ',' << seconds << '\n';
    // A line-only preconditioner is a baseline, not a multigrid guarantee;
    // at large N it may hit the strict iteration cap. Record its failure
    // rather than suppressing data or relaxing the qualification threshold.
    if (qualification_required) {
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(report.physical_residual_relative) &&
                    report.physical_residual_relative <= tolerance);
        EXPECT_TRUE(gauge < 1e-7);
        EXPECT_TRUE(err < 2e-3);
        EXPECT_TRUE(flux_err < 2e-8);
    } else {
        EXPECT_TRUE(std::isfinite(report.physical_residual_relative));
        EXPECT_TRUE(report.result.status == SolverStatus::MAX_ITER_REACHED ||
                    report.result.status == SolverStatus::CONVERGED);
    }
    return {report.result.iterations, report.physical_residual_relative,
            err, flux_err, seconds, report.result.status};
}

} // namespace

int main(int argc, char** argv) {
    const bool benchmark = argc > 1 && std::string(argv[1]) == "--benchmark";
    const std::size_t nx = benchmark && argc > 2 ? std::stoul(argv[2]) : 24;
    const std::size_t ny = benchmark && argc > 3 ? std::stoul(argv[3]) : 18;
    if (nx < 3 || ny < 3 || nx * ny > 65536) return 2;

    run_case("singular_strongly_stretched_pressure_both_directions", [&] {
        for (const bool strong_x : {true, false}) {
            for (const double ratio : {1e3, 1e6}) {
                const VolumeMesh grid(nx, ny, ratio, strong_x);
                const auto exact = grid.analytic_pressure();
                const auto rhs = apply(grid.matrix, exact);
                const auto projector = NullSpaceProjector::constant(nx * ny);
                EXPECT_TRUE(projector.is_null_space(grid.matrix));
                EXPECT_TRUE(projector.is_compatible(rhs));
                EXPECT_TRUE(balance_error(grid.face_balance(exact), rhs) < 1e-12);
                double net_source = 0.0;
                for (std::size_t k = 0; k < rhs.size(); ++k)
                    net_source += rhs(k);
                EXPECT_TRUE(std::abs(net_source) <
                    1e-10 * std::max(1.0, rhs.norm2()));
                for (const auto pc : {
                    PreconditionerModel::AnisotropicLine,
                    PreconditionerModel::AnisotropicBoomerAMG,
                    PreconditionerModel::AnisotropicSmoothedAggregationAMG}) {
                    run_solve(grid, rhs, exact, KrylovModel::FGMRES, pc,
                              ratio, strong_x, 1e-9,
                              pc != PreconditionerModel::AnisotropicLine ||
                              nx * ny <= 600);
                }
                run_solve(grid, rhs, exact, KrylovModel::CG,
                          PreconditionerModel::NativeAMG, ratio, strong_x);
            }
        }
    });

    run_case("rotated_strong_coupling_neumann_preserves_gauge", [&] {
        // Graph edges along northeast lines emulate a rotated anisotropic
        // diffusion tensor. No artificial reaction or gauge row is added.
        constexpr std::size_t n = 30;
        SparseMatrix A(n * n, n * n);
        std::vector<double> diag(n * n, 0.0);
        const auto edge = [&](std::size_t i, std::size_t j, double w) {
            A.push_back(i, j, -w);
            A.push_back(j, i, -w);
            diag[i] += w;
            diag[j] += w;
        };
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t cell = j * n + i;
                if (i + 1 < n) edge(cell, cell + 1, 1.0);
                if (j + 1 < n) edge(cell, cell + n, 1.0);
                if (i + 1 < n && j + 1 < n)
                    edge(cell, cell + n + 1, 1e4);
            }
        for (std::size_t i = 0; i < diag.size(); ++i)
            A.push_back(i, i, diag[i]);
        A.finalize();
        const auto space = NullSpaceProjector::constant(n * n);
        EXPECT_TRUE(is_bilateral_null_space(A, space));
        Vector exact(n * n);
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                exact(j * n + i) = std::sin(pi * (static_cast<double>(i) + 0.5) /
                    static_cast<double>(n)) *
                    std::cos(2.0 * pi * (static_cast<double>(j) + 0.5) /
                    static_cast<double>(n));
        space.remove(exact);
        const Vector rhs = apply(A, exact);
        for (const auto pc : {
            PreconditionerModel::AnisotropicBoomerAMG,
            PreconditionerModel::AnisotropicSmoothedAggregationAMG}) {
            LinearSolverRequest req;
            req.null_space = NullSpaceModel::Constant;
            req.krylov = KrylovModel::FGMRES;
            req.preconditioner = pc;
            Vector x(n * n, 3.0);
            const auto report = solve_linear_system(A, rhs, x,
                LinearProblemKind::PressurePoisson, req, 350, 1e-9);
            const double error = gauge_free_error(x, exact);
            std::cout << std::setprecision(12) << "rotated_neumann_csv,"
                      << to_string(pc) << ',' << report.result.iterations << ','
                      << report.physical_residual_relative << ',' << error << '\n';
            EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
            EXPECT_TRUE(report.physical_residual_relative <= 1e-9);
            EXPECT_TRUE(error <= 1e-3);
            EXPECT_TRUE(space.component_norm(x) < 1e-9);
        }
    });

    run_case("projected_gmres_uses_singular_amg", [&] {
        const VolumeMesh grid(16, 12, 1e3, true);
        const auto exact = grid.analytic_pressure();
        const auto rhs = apply(grid.matrix, exact);
        run_solve(grid, rhs, exact, KrylovModel::GMRES,
                  PreconditionerModel::SmoothedAggregationAMG, 1e3, true);
    });

    run_case("nonzero_source_mean_is_rejected_without_mutation", [&] {
        const VolumeMesh grid(10, 9, 1e3, true);
        auto rhs = apply(grid.matrix, grid.analytic_pressure());
        rhs(0) += 1.0;
        Vector solution(rhs.size(), 9.5);
        LinearSolverRequest request;
        request.null_space = NullSpaceModel::Constant;
        request.krylov = KrylovModel::FGMRES;
        request.preconditioner = PreconditionerModel::AnisotropicSmoothedAggregationAMG;
        const auto report = solve_linear_system(
            grid.matrix, rhs, solution, LinearProblemKind::PressurePoisson,
            request, 100, 1e-9);
        EXPECT_TRUE(report.result.status == SolverStatus::NOT_APPLICABLE);
        for (std::size_t i = 0; i < solution.size(); ++i)
            EXPECT_NEAR(solution(i), 9.5, 0.0);
    });

    run_case("homogeneous_pressure_residual_is_truly_certified", [&] {
        const VolumeMesh grid(12, 10, 1e4, true);
        Vector rhs(grid.nx * grid.ny, 0.0);
        Vector solution = grid.analytic_pressure();
        LinearSolverRequest request;
        request.null_space = NullSpaceModel::Constant;
        request.krylov = KrylovModel::FGMRES;
        request.preconditioner = PreconditionerModel::AnisotropicSmoothedAggregationAMG;
        const double initial_norm = physical_residual_norm(grid.matrix, solution, rhs);
        const auto report = solve_linear_system(grid.matrix, rhs, solution,
            LinearProblemKind::PressurePoisson, request, 350, 1e-9);
        const double final_norm = physical_residual_norm(grid.matrix, solution, rhs);
        std::cout << "homogeneous_neumann_csv," << initial_norm << ','
                  << final_norm << ',' << report.physical_residual_relative << '\n';
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(final_norm <= 1e-9 * initial_norm);
        EXPECT_TRUE(report.physical_residual_relative <= 1e-9);
        EXPECT_TRUE(NullSpaceProjector::constant(solution.size()).component_norm(solution) < 1e-10);
        solution.fill(0.0);
        const auto already_solved = solve_linear_system(grid.matrix, rhs, solution,
            LinearProblemKind::PressurePoisson, request, 200, 1e-9);
        EXPECT_TRUE(already_solved.result.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(already_solved.physical_residual_relative, 0.0, 0.0);
    });

    run_case("incorrect_right_and_left_nullspaces_fail_closed", [&] {
        const VolumeMesh grid(10, 9, 1e3, true);
        const auto exact = grid.analytic_pressure();
        const auto rhs = apply(grid.matrix, exact);
        LinearSolverRequest request;
        request.null_space = NullSpaceModel::Constant;
        request.krylov = KrylovModel::FGMRES;
        request.preconditioner = PreconditionerModel::Jacobi;

        SparseMatrix anchored(10 * 9, 10 * 9);
        const auto* row = grid.matrix.row_offsets_data();
        const auto* col = grid.matrix.columns_data();
        const auto* val = grid.matrix.values_data();
        for (std::size_t i = 0; i < anchored.n_rows(); ++i)
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                anchored.push_back(i, col[k], val[k] + (i == 0 && col[k] == i ? 1.0 : 0.0));
        anchored.finalize();
        Vector x(rhs.size(), 2.0);
        auto result = solve_linear_system(anchored, rhs, x,
            LinearProblemKind::PressurePoisson, request, 100, 1e-9);
        EXPECT_TRUE(result.result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_NEAR(x(0), 2.0, 0.0);

        // Multiplying each matrix row by a positive nonuniform factor leaves
        // the right constant mode intact but destroys the left mode.
        SparseMatrix wrong_left(10 * 9, 10 * 9);
        for (std::size_t i = 0; i < wrong_left.n_rows(); ++i)
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                wrong_left.push_back(i, col[k], val[k] * (i % 2 ? 3.0 : 1.0));
        wrong_left.finalize();
        const auto nullspace = NullSpaceProjector::constant(x.size());
        EXPECT_TRUE(nullspace.is_null_space(wrong_left));
        EXPECT_FALSE(is_bilateral_null_space(wrong_left, nullspace));
        x.fill(2.0);
        result = solve_linear_system(wrong_left, rhs, x,
            LinearProblemKind::PressurePoisson, request, 100, 1e-9);
        EXPECT_TRUE(result.result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_NEAR(x(0), 2.0, 0.0);
    });

    run_case("automatic_selection_with_line_amg_nullspace", [&] {
        const VolumeMesh grid(16, 12, 1e4, false);
        auto exact = grid.analytic_pressure();
        const auto rhs = apply(grid.matrix, exact);
        LinearSolverRequest request;
        request.null_space = NullSpaceModel::Constant;
        request.preconditioner = PreconditionerModel::AnisotropicSmoothedAggregationAMG;
        Vector x(rhs.size(), 0.0);
        const auto result = solve_linear_system(grid.matrix, rhs, x,
            LinearProblemKind::PressurePoisson, request, 350, 1e-9);
        EXPECT_TRUE(result.plan.krylov == KrylovModel::FGMRES);
        EXPECT_TRUE(result.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(gauge_free_error(x, exact) < 2e-3);
        request.krylov = KrylovModel::CG;
        EXPECT_THROW(select_linear_solver(LinearProblemKind::PressurePoisson,
            measure_matrix_characteristics(grid.matrix), request), std::invalid_argument);
    });

    run_case("manufactured_continuum_operator_refinement", [&] {
        // Independent continuum MMS: -div(K grad p) = pi^2*(Kx+4*Ky)*p
        // with exact homogeneous Neumann conditions. The midpoint source is
        // not globally compatible to roundoff on a nonuniform grid, so this
        // test compares the discretization defect without altering the RHS.
        double coarse_error = 0.0;
        double coarse_solution_error = 0.0;
        for (const std::size_t n : {12U, 24U, 48U}) {
            const VolumeMesh grid(n, n, 100.0, true, 1.3);
            const Vector exact = grid.analytic_pressure();
            const Vector discrete = apply(grid.matrix, exact);
            const Vector continuum = grid.source_quadrature(100.0, true);
            const double defect = balance_error(discrete, continuum);
            std::cout << std::setprecision(12) << "continuum_mms_csv,"
                      << n << ',' << defect << '\n';
            EXPECT_TRUE(std::isfinite(defect));
            if (coarse_error > 0.0)
                EXPECT_TRUE(defect < 0.65 * coarse_error);
            coarse_error = defect;

            // The continuum quadrature has a nonzero net discrete source on
            // graded cells. Compute that compatibility defect explicitly,
            // then use the physically necessary zero-integral forcing, not
            // a silently modified RHS inside the linear solver.
            Vector compatible = continuum;
            const double source_mean = [&]() {
                double sum = 0.0;
                for (std::size_t k = 0; k < compatible.size(); ++k)
                    sum += compatible(k);
                return sum / static_cast<double>(compatible.size());
            }();
            const double source_balance_adjustment =
                std::abs(source_mean) * std::sqrt(static_cast<double>(compatible.size())) /
                std::max(continuum.norm2(), 1e-300);
            for (std::size_t k = 0; k < compatible.size(); ++k)
                compatible(k) -= source_mean;
            LinearSolverRequest request;
            request.krylov = KrylovModel::FGMRES;
            request.preconditioner = PreconditionerModel::AnisotropicSmoothedAggregationAMG;
            request.null_space = NullSpaceModel::Constant;
            Vector result(n * n, 0.0);
            const auto report = solve_linear_system(grid.matrix, compatible, result,
                LinearProblemKind::PressurePoisson, request, 400, 1e-10);
            const double sol_error = gauge_free_error(result, exact);
            std::cout << std::setprecision(12) << "continuum_solution_csv,"
                      << n << ',' << source_balance_adjustment << ','
                      << report.physical_residual_relative << ','
                      << sol_error << ',' << report.result.iterations << '\n';
            EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
            EXPECT_TRUE(report.physical_residual_relative <= 1e-10);
            EXPECT_TRUE(std::isfinite(sol_error));
            if (coarse_solution_error > 0.0)
                EXPECT_TRUE(sol_error < 0.75 * coarse_solution_error);
            coarse_solution_error = sol_error;
        }
    });

    return run_all();
}