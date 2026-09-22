#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/chebyshev_smoother.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/sparse_matrix.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>



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

    // Setup must reject malformed matrix data rather than silently falling
    // back to an identity diagonal.
    SparseMatrix bad(4, 4);
    bad.push_back(0, 0, std::numeric_limits<double>::quiet_NaN());
    bad.push_back(0, 1, -1.0);
    bad.push_back(1, 0, -1.0);
    bad.push_back(1, 1, 2.0);
    bad.push_back(1, 2, -1.0);
    bad.push_back(2, 1, -1.0);
    bad.push_back(2, 2, 2.0);
    bad.push_back(2, 3, -1.0);
    bad.push_back(3, 2, -1.0);
    bad.push_back(3, 3, 2.0);
    bad.finalize();
    if (amg.setup(bad)) {
        return 7;
    }

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
