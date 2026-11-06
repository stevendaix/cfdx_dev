// Regression gates for anisotropic line preconditioners. A solver may reuse a
// factorization only while its directional graph remains representative.
#include "cfdx/core/linalg/algebraic_line_preconditioner.h"
#include "cfdx/core/linalg/anisotropic_line_amg.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix shifted_diffusion(std::size_t nx, std::size_t ny,
                               double x_strength, double y_strength,
                               double oblique_strength) {
    const std::size_t n = nx * ny;
    SparseMatrix A(n, n);
    for (std::size_t y = 0; y < ny; ++y)
        for (std::size_t x = 0; x < nx; ++x) {
            const std::size_t i = y * nx + x;
            double diag = 1.0;
            const auto add = [&](std::size_t j, double coefficient) {
                A.push_back(i, j, -coefficient);
                diag += coefficient;
            };
            if (x > 0) add(i - 1, x_strength);
            if (x + 1 < nx) add(i + 1, x_strength);
            if (y > 0) add(i - nx, y_strength);
            if (y + 1 < ny) add(i + nx, y_strength);
            if (x > 0 && y > 0) add(i - nx - 1, oblique_strength);
            if (x + 1 < nx && y + 1 < ny)
                add(i + nx + 1, oblique_strength);
            if (x + 1 < nx && y > 0)
                add(i - nx + 1, oblique_strength);
            if (x > 0 && y + 1 < ny)
                add(i + nx - 1, oblique_strength);
            A.push_back(i, i, diag);
        }
    A.finalize();
    return A;
}

double relative_inf(const Vector& x, const Vector& y) {
    if (x.size() != y.size()) return 1e100;
    double worst = 0.0;
    double magnitude = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        worst = std::max(worst, std::abs(x(i) - y(i)));
        magnitude = std::max(magnitude, std::abs(y(i)));
    }
    return worst / std::max(1.0, magnitude);
}
} // namespace

int main() {
    run_case("line_reuse_rejects_rotated_dominant_strength_even_with_weak_edge", [] {
        // The tiny oblique edge keeps the min/max row ratio above 32, while
        // x/y strengths swap and the original x-line is no longer dominant.
        const auto A = shifted_diffusion(12, 9, 100.0, 1.0, 0.01);
        const auto B = shifted_diffusion(12, 9, 10.0, 100.0, 0.01);
        AlgebraicLinePreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        EXPECT_TRUE(pc.covered_rows() >= A.n_rows() / 2);
        Vector rhs(A.n_rows(), 1.0), before(A.n_rows(), 0.0), after(A.n_rows(), 0.0);
        EXPECT_TRUE(pc.apply(rhs, before));
        EXPECT_TRUE(!pc.update_values(B));
        EXPECT_TRUE(pc.apply(rhs, after));
        EXPECT_TRUE(relative_inf(after, before) <= 1e-14);
    });
    run_case("failed_full_setup_keeps_previously_certified_line_pc", [] {
        const auto A = shifted_diffusion(12, 9, 1000.0, 1.0, 0.01);
        AlgebraicLinePreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        Vector rhs(A.n_rows(), 1.0), before(A.n_rows(), 0.0), after(A.n_rows(), 0.0);
        EXPECT_TRUE(pc.apply(rhs, before));
        const SparseMatrix invalid;
        EXPECT_TRUE(!pc.setup(invalid));
        EXPECT_TRUE(pc.apply(rhs, after));
        EXPECT_TRUE(relative_inf(after, before) <= 1e-14);
    });
    run_case("line_amg_failed_setup_preserves_previous_hierarchy", [] {
        const auto A = shifted_diffusion(16, 12, 1000.0, 1.0, 0.01);
        AnisotropicLineAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        Vector rhs(A.n_rows(), 1.0), before(A.n_rows(), 0.0), after(A.n_rows(), 0.0);
        EXPECT_TRUE(pc.apply(rhs, before));
        const SparseMatrix invalid;
        EXPECT_TRUE(!pc.setup(invalid));
        EXPECT_TRUE(pc.apply(rhs, after));
        EXPECT_TRUE(relative_inf(after, before) <= 1e-14);
    });
    run_case("line_amg_apply_supports_in_place_rhs", [] {
        const auto A = shifted_diffusion(16, 12, 1000.0, 1.0, 0.01);
        AnisotropicLineAMGPreconditioner pc;
        EXPECT_TRUE(pc.setup(A));
        Vector rhs(A.n_rows(), 0.0), out(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i)
            rhs(i) = 1.0 + std::sin(static_cast<double>(i) * 0.13);
        EXPECT_TRUE(pc.apply(rhs, out));
        Vector alias = rhs;
        EXPECT_TRUE(pc.apply(alias, alias));
        EXPECT_TRUE(relative_inf(alias, out) <= 3e-12);
    });
    return run_all();
}