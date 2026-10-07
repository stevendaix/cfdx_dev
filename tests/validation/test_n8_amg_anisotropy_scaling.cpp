#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_anisotropic_diffusion(std::size_t nx,
                                        std::size_t ny,
                                        double ax,
                                        double ay,
                                        double reaction) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diagonal = reaction;
            if (x > 0) {
                A.push_back(i, i - 1, -ax);
                diagonal += ax;
            }
            if (x + 1 < nx) {
                A.push_back(i, i + 1, -ax);
                diagonal += ax;
            }
            if (y > 0) {
                A.push_back(i, i - nx, -ay);
                diagonal += ay;
            }
            if (y + 1 < ny) {
                A.push_back(i, i + nx, -ay);
                diagonal += ay;
            }
            A.push_back(i, i, diagonal);
        }
    }
    A.finalize();
    return A;
}

Vector make_rhs(std::size_t n) {
    Vector b(n);
    for (std::size_t i = 0; i < n; ++i)
        b(i) = 1.0 + 0.25 * std::sin(0.13 * static_cast<double>(i));
    return b;
}

double relative_true_residual(const SparseMatrix& A,
                              const Vector& x,
                              const Vector& b) {
    const auto ax = A.matvec(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        const double r = b(i) - ax[i];
        r2 += r * r;
        b2 += b(i) * b(i);
    }
    return std::sqrt(r2 / std::max(b2, 1e-300));
}

struct QualificationRow {
    double anisotropy;
    double scale;
    std::size_t iterations;
    double true_residual;
};

template <typename AMG>
QualificationRow run_amg_case(const SparseMatrix& A,
                              const Vector& b,
                              const char* method,
                              double anisotropy,
                              double scale,
                              double tolerance,
                              bool& qualification_ok) {
    AMG amg;
    const bool setup_ok = amg.setup(A);
    const auto levels = amg.hierarchy_level_sizes();

    std::cout << "n8_amg_sweep"
              << " method=" << method
              << " anisotropy=" << anisotropy
              << " scale=" << scale
              << " setup=" << (setup_ok ? "PASS" : "FAIL")
              << " levels=" << levels.size()
              << " fine=" << A.n_rows();

    if (!setup_ok) {
        std::cout << " qualification=FAIL reason=setup
";
        qualification_ok = false;
        return {anisotropy, scale, 0, std::numeric_limits<double>::infinity()};
    }

    Vector z(A.n_rows(), 0.0);
    const bool apply_ok = amg.apply(b, z);
    const double apply_residual = relative_true_residual(A, z, b);

    Vector x(A.n_rows(), 0.0);
    const auto result = solve_cg(A, b, x, amg, 5000, tolerance);
    const double true_residual = relative_true_residual(A, x, b);

    const bool finite_evidence = std::isfinite(apply_residual) &&
                                 std::isfinite(true_residual);
    const bool hierarchy_ok = levels.size() >= 2 &&
                              amg.first_prolongation_nnz() > 0;
    const bool converged = result.status == SolverStatus::CONVERGED;
    // The solve tolerance is the acceptance contract; no independent
    // iteration-count threshold or arbitrary residual margin is introduced.
    const bool residual_ok = finite_evidence && true_residual <= tolerance;
    const bool case_ok = apply_ok && hierarchy_ok && finite_evidence &&
                         converged && residual_ok;

    std::cout << " apply=" << (apply_ok ? "PASS" : "FAIL")
              << " hierarchy=" << (hierarchy_ok ? "PASS" : "FAIL")
              << " P_nnz=" << amg.first_prolongation_nnz()
              << " apply_true_residual=" << apply_residual
              << " cg_status=" << static_cast<int>(result.status)
              << " cg_iterations=" << result.iterations
              << " cg_true_residual=" << true_residual
              << " tolerance=" << tolerance
              << " qualification=" << (case_ok ? "PASS" : "FAIL")
              << '
';

    qualification_ok = qualification_ok && case_ok;
    return {anisotropy, scale, result.iterations, true_residual};
}

} // namespace

int main() {
    std::cout << std::setprecision(17);

    constexpr std::size_t nx = 32;
    constexpr std::size_t ny = 32;
    constexpr double reaction = 1.0;
    constexpr double anisotropies[] = {1.0, 10.0, 100.0, 1000.0};
    constexpr double scales[] = {1e-3, 1.0, 1e3};

    constexpr double cg_tolerance = 1e-9;
    bool qualification_ok = true;
    std::vector<QualificationRow> rows;

    for (const double anisotropy : anisotropies) {
        for (const double scale : scales) {
            const double ax = scale * anisotropy;
            const double ay = scale;
            const auto A = make_anisotropic_diffusion(
                nx, ny, ax, ay, reaction * scale);
            const auto b = make_rhs(A.n_rows());

            std::cout << "n8_amg_operator"
                      << " anisotropy=" << anisotropy
                      << " scale=" << scale
                      << " ax=" << ax
                      << " ay=" << ay
                      << " rows=" << A.n_rows()
                      << " cols=" << A.n_cols() << '
';

            rows.push_back(run_amg_case<NativeBoomerAMGPreconditioner>(
                A, b, "NativeBoomerAMG", anisotropy, scale, cg_tolerance,
                qualification_ok));
        }
    }

    EXPECT_TRUE(structural_ok);
    std::cout << "N8 AMG anisotropy/scaling qualification: "
              << (structural_ok ? "PASS" : "FAIL") << '
';
    return 0;
}
