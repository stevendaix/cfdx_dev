#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/chebyshev_smoother.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

using cfdx::core::FunctionalLinearOperator;
using cfdx::core::MatrixFreeVcyclePreconditioner;
using cfdx::core::SparseMatrix;
using cfdx::core::Vector;



static SparseMatrix make_poisson_1d(std::size_t n) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

static SparseMatrix make_poisson_2d(std::size_t nx, std::size_t ny) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            A.push_back(i, i, 4.0);
            if (x > 0) A.push_back(i, i - 1, -1.0);
            if (x + 1 < nx) A.push_back(i, i + 1, -1.0);
            if (y > 0) A.push_back(i, i - nx, -1.0);
            if (y + 1 < ny) A.push_back(i, i + nx, -1.0);
        }
    }
    A.finalize();
    return A;
}

static double true_residual_ratio(const SparseMatrix& A, const Vector& rhs,
                                  const Vector& correction) {
    const auto Az = A.matvec(correction);
    double before2 = 0.0;
    double after2 = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        before2 += rhs(i) * rhs(i);
        const double ri = rhs(i) - Az[i];
        after2 += ri * ri;
    }
    return std::sqrt(after2 / before2);
}

static bool check_amg(
    const SparseMatrix& A,
    Vector rhs,
    double max_ratio,
    cfdx::core::AMGInterpolationPolicy interpolation =
        cfdx::core::AMGInterpolationPolicy::DirectCF) {
    FunctionalLinearOperator op(
        A.n_rows(),
        [&A](const Vector& x, Vector& y) {
            const auto result = A.matvec(x);
            for (std::size_t i = 0; i < result.size(); ++i) y(i) = result[i];
        });

    MatrixFreeVcyclePreconditioner amg(
        op, 0.7, 4, 4, 0.25, 25, interpolation);
    if (!amg.setup(A)) return false;

    Vector correction;
    if (!amg.apply(rhs, correction) || !correction.is_valid()) return false;
    const double ratio = true_residual_ratio(A, rhs, correction);
    std::cerr << "AMG residual ratio=" << ratio << " threshold=" << max_ratio << "\\n";
    return std::isfinite(ratio) && ratio < max_ratio;
}

static SparseMatrix make_anisotropic_diffusion_2d(std::size_t nx, std::size_t ny,
                                                   double ax, double ay) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diag = 0.0;
            if (x > 0) { A.push_back(i, i - 1, -ax); diag += ax; }
            if (x + 1 < nx) { A.push_back(i, i + 1, -ax); diag += ax; }
            if (y > 0) { A.push_back(i, i - nx, -ay); diag += ay; }
            if (y + 1 < ny) { A.push_back(i, i + nx, -ay); diag += ay; }
            A.push_back(i, i, diag + 1.0);
        }
    }
    // Small positive reaction term removes the pure-Neumann nullspace while
    // retaining strong directional anisotropy.
    A.finalize();
    return A;
}

// Finite-volume-style cell-centred diffusion assembly: each internal face
// contributes equal/opposite off-diagonal flux coefficients to its two cells.
static SparseMatrix make_fvm_diffusion_2d(std::size_t nx, std::size_t ny,
                                           double kx, double ky) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y) {
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diag = 0.0;
            if (x > 0) { A.push_back(i, i - 1, -kx); diag += kx; }
            if (x + 1 < nx) { A.push_back(i, i + 1, -kx); diag += kx; }
            if (y > 0) { A.push_back(i, i - nx, -ky); diag += ky; }
            if (y + 1 < ny) { A.push_back(i, i + nx, -ky); diag += ky; }
            A.push_back(i, i, diag + 1.0);
        }
    }
    A.finalize();
    return A;
}

int main() {
    using namespace cfdx::core;

    SparseMatrix A(4, 4);
    A.push_back(0, 0, 2.0);
    A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -1.0);
    A.push_back(1, 1, 2.0);
    A.push_back(1, 2, -1.0);
    A.push_back(2, 1, -1.0);
    A.push_back(2, 2, 2.0);
    A.push_back(2, 3, -1.0);
    A.push_back(3, 2, -1.0);
    A.push_back(3, 3, 2.0);
    A.finalize();

    FunctionalLinearOperator op(
        4,
        [&A](const Vector& x, Vector& y) {
            const auto result = A.matvec(x);
            for (std::size_t i = 0; i < result.size(); ++i) {
                y(i) = result[i];
            }
        });

    MatrixFreeVcyclePreconditioner amg(op);
    if (!amg.setup(A)) {
        return 1;
    }

    if (amg.coarse_size() != 2 ||
        amg.aggregate_of(0) != amg.aggregate_of(1) ||
        amg.aggregate_of(2) != amg.aggregate_of(3) ||
        amg.aggregate_of(0) == amg.aggregate_of(2)) {
        return 2;
    }

    Vector rhs(4, 1.0);
    Vector correction;
    if (!amg.apply(rhs, correction) || !correction.is_valid()) {
        return 3;
    }
    const auto Az = A.matvec(correction);
    double before2 = 0.0;
    double after2 = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double before = rhs(i);
        const double after = rhs(i) - Az[i];
        before2 += before * before;
        after2 += after * after;
    }
    const double before_norm = std::sqrt(before2);
    const double after_norm = std::sqrt(after2);
    if (!(after_norm < before_norm) || !std::isfinite(after_norm)) {
        std::cerr << "AMG did not reduce the true residual: "
                  << before_norm << " -> " << after_norm << "\n";
        return 6;
    }

    MatrixFreeVcyclePreconditioner smoothed_aggregation(
        op, 0.7, 4, 4, 0.25, 25,
        AMGInterpolationPolicy::SmoothedAggregation);
    if (!smoothed_aggregation.setup(A) ||
        smoothed_aggregation.coarse_size() == 0 ||
        smoothed_aggregation.coarse_size() >= A.n_rows()) {
        return 30;
    }
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        if (std::abs(smoothed_aggregation.prolongation_row_sum(i) - 1.0) >
            1e-12) {
            std::cerr << "Smoothed aggregation did not preserve constants\n";
            return 31;
        }
    }
    Vector sa_correction;
    if (!smoothed_aggregation.apply(rhs, sa_correction) ||
        true_residual_ratio(A, rhs, sa_correction) >= 1.0) {
        std::cerr << "Smoothed aggregation did not reduce the true residual\n";
        return 32;
    }

    // Quantitative multilevel checks requested by Phase 4.7.
    const SparseMatrix poisson1d = make_poisson_1d(4);
    Vector rhs1d(4);
    for (std::size_t i = 0; i < rhs1d.size(); ++i) rhs1d(i) = (i % 2 == 0) ? 1.0 : -1.0;
    if (!check_amg(poisson1d, rhs1d, 1.0)) {
        std::cerr << "1D Poisson AMG residual reduction failed\\n";
        return 9;
    }

    const SparseMatrix poisson2d = make_poisson_2d(4, 4);
    Vector rhs2d(16);
    for (std::size_t i = 0; i < rhs2d.size(); ++i) rhs2d(i) = ((i % 4 + i / 4) % 2 == 0) ? 1.0 : -1.0;
    if (!check_amg(poisson2d, rhs2d, 0.999)) {
        std::cerr << "2D Poisson AMG residual reduction failed\\n";
        return 10;
    }
    if (!check_amg(
            poisson2d, rhs2d, 0.999,
            AMGInterpolationPolicy::SmoothedAggregation)) {
        std::cerr << "2D Poisson smoothed aggregation residual reduction failed\n";
        return 33;
    }

    const SparseMatrix anisotropic = make_anisotropic_diffusion_2d(16, 16, 1.0, 1000.0);
    Vector rhs_aniso(256, 1.0);
    if (!check_amg(anisotropic, rhs_aniso, 0.99)) {
        std::cerr << "Strongly anisotropic AMG regression failed\\n";
        return 12;
    }
    if (!check_amg(
            anisotropic, rhs_aniso, 0.99,
            AMGInterpolationPolicy::SmoothedAggregation)) {
        std::cerr << "Anisotropic smoothed aggregation regression failed\n";
        return 34;
    }

    // Exercise repeated V-cycles rather than only a single residual reduction.
    FunctionalLinearOperator aniso_op(
        anisotropic.n_rows(),
        [&anisotropic](const Vector& x, Vector& y) {
            const auto result = anisotropic.matvec(x);
            for (std::size_t i = 0; i < result.size(); ++i) y(i) = result[i];
        });
    MatrixFreeVcyclePreconditioner repeated(aniso_op);
    if (!repeated.setup(anisotropic)) return 13;
    Vector repeated_rhs(256, 1.0);
    Vector repeated_x(256, 0.0);
    const double initial = std::sqrt(256.0);
    double previous = initial;
    for (std::size_t cycle = 0; cycle < 4; ++cycle) {
        const auto Ax = anisotropic.matvec(repeated_x);
        Vector residual(repeated_rhs.size());
        double rhs_norm2 = 0.0;
        double residual_norm2 = 0.0;
        for (std::size_t i = 0; i < repeated_rhs.size(); ++i) {
            residual(i) = repeated_rhs(i) - Ax[i];
            rhs_norm2 += repeated_rhs(i) * repeated_rhs(i);
            residual_norm2 += residual(i) * residual(i);
        }

        Vector correction;
        if (!repeated.apply(residual, correction)) return 14;
        for (std::size_t i = 0; i < repeated_x.size(); ++i) {
            repeated_x(i) += correction(i);
        }

        const auto updated = anisotropic.matvec(repeated_x);
        double current2 = 0.0;
        for (std::size_t i = 0; i < repeated_rhs.size(); ++i) {
            const double ri = repeated_rhs(i) - updated[i];
            current2 += ri * ri;
        }
        const double current = std::sqrt(current2);
        if (!std::isfinite(current) || current >= previous) return 15;
        previous = current;

        // The first cycle must actually reduce the original residual, while
        // subsequent cycles must continue solving the current residual.
        if (cycle == 0 && !(current < std::sqrt(rhs_norm2))) return 15;
        (void)residual_norm2;
    }

    const SparseMatrix fvm_diffusion = make_fvm_diffusion_2d(16, 16, 1.0, 20.0);
    Vector rhs_fvm(256, 1.0);
    if (!check_amg(fvm_diffusion, rhs_fvm, 0.95)) {
        std::cerr << "FVM diffusion AMG regression failed\\n";
        return 16;
    }

    // Duplicate diagonal entries are a normal FVM assembly pattern and must
    // be accumulated by the Galerkin construction.
    SparseMatrix duplicate_diag(4, 4);
    duplicate_diag.push_back(0, 0, 1.0);
    duplicate_diag.push_back(0, 0, 1.0);
    duplicate_diag.push_back(0, 1, -1.0);
    duplicate_diag.push_back(1, 0, -1.0);
    duplicate_diag.push_back(1, 1, 2.0);
    duplicate_diag.push_back(1, 2, -1.0);
    duplicate_diag.push_back(2, 1, -1.0);
    duplicate_diag.push_back(2, 2, 2.0);
    duplicate_diag.push_back(2, 3, -1.0);
    duplicate_diag.push_back(3, 2, -1.0);
    duplicate_diag.push_back(3, 3, 2.0);
    duplicate_diag.finalize();
    Vector rhs_dup(4, 1.0);
    if (!check_amg(duplicate_diag, rhs_dup, 0.95)) {
        std::cerr << "Duplicate-diagonal AMG regression failed\\n";
        return 11;
    }

    ChebyshevSmoother::Controls controls;
    controls.lambda_max = 0.0;
    controls.spectral_iterations = 4;

    ChebyshevSmoother smoother(controls);
    const double lambda_max = smoother.lambda_max(op);

    if (!(lambda_max > 0.0) || !std::isfinite(lambda_max)) {
        return 4;
    }

    Vector x(4, 0.0);
    smoother.apply(op, rhs, x);

    if (!x.is_valid()) {
        return 5;
    }

    // Malformed/non-finite matrix rejection is covered by the core SparseMatrix
    // and preconditioner regression suites; keep this AMG test focused on hierarchy
    // construction and quantitative convergence.

    SparseMatrix missing_diag(4, 4);
    missing_diag.push_back(0, 1, -1.0);
    missing_diag.push_back(1, 0, -1.0);
    missing_diag.push_back(1, 1, 2.0);
    missing_diag.push_back(1, 2, -1.0);
    missing_diag.push_back(2, 1, -1.0);
    missing_diag.push_back(2, 2, 2.0);
    missing_diag.push_back(2, 3, -1.0);
    missing_diag.push_back(3, 2, -1.0);
    missing_diag.push_back(3, 3, 2.0);
    missing_diag.finalize();
    if (amg.setup(missing_diag)) {
        return 8;
    }

    std::cout << "AMG/Chebyshev runtime: PASS\n";
    return 0;
}
