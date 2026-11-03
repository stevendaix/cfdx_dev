#pragma once

#include "cfdx/physics/wall_face_reconstruction.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// A sample is a cell-centred value or an exact boundary value. Include the
// wall-face sample with prescribed no-slip velocity in the stencil.
struct WallQuadraticSample {
    cfdx::core::Vec3 position;
    double pressure = 0.0;
    cfdx::core::Vec3 velocity;
};

namespace detail {
inline std::array<double, 10> wall_quadratic_basis(
    const cfdx::core::Vec3& x, const cfdx::core::Vec3& origin, double scale)
{
    const double a = (x.x - origin.x) / scale;
    const double b = (x.y - origin.y) / scale;
    const double c = (x.z - origin.z) / scale;
    return {1.0, a, b, c, a*a, a*b, a*c, b*b, b*c, c*c};
}

inline std::array<double, 10> solve_wall_quadratic_fit(
    const std::vector<WallQuadraticSample>& samples,
    const cfdx::core::Vec3& origin, double scale, int component)
{
    constexpr std::size_t n = 10;
    // Solve the weighted least-squares problem with column-pivoted
    // Householder QR. Forming A^T A would square the stencil condition
    // number and make rank decisions unnecessarily fragile.
    std::vector<std::array<double, n>> a(samples.size());
    std::vector<double> b(samples.size(), 0.0);
    std::array<std::size_t, n> permutation{};
    for (std::size_t j = 0; j < n; ++j) permutation[j] = j;

    for (std::size_t row = 0; row < samples.size(); ++row) {
        const auto& sample = samples[row];
        const auto basis = wall_quadratic_basis(sample.position, origin, scale);
        const double value = component == 0 ? sample.pressure :
            (component == 1 ? sample.velocity.x :
             (component == 2 ? sample.velocity.y : sample.velocity.z));
        const double r2 = (sample.position - origin).mag2() / (scale*scale);
        const double sqrt_weight = 1.0 / std::sqrt(1.0 + r2);
        for (std::size_t j = 0; j < n; ++j)
            a[row][j] = sqrt_weight * basis[j];
        b[row] = sqrt_weight * value;
    }

    std::array<double, n> diagonal{};
    for (std::size_t k = 0; k < n; ++k) {
        // Column pivoting selects the largest remaining column norm.
        std::size_t pivot = k;
        double pivot_norm2 = -1.0;
        for (std::size_t j = k; j < n; ++j) {
            double norm2 = 0.0;
            for (std::size_t i = k; i < samples.size(); ++i)
                norm2 += a[i][j] * a[i][j];
            if (norm2 > pivot_norm2) {
                pivot_norm2 = norm2;
                pivot = j;
            }
        }
        if (!(pivot_norm2 > 0.0) || !std::isfinite(pivot_norm2))
            throw std::invalid_argument(
                "wall quadratic reconstruction: rank-deficient stencil");
        if (pivot != k) {
            for (std::size_t i = 0; i < samples.size(); ++i)
                std::swap(a[i][k], a[i][pivot]);
            std::swap(permutation[k], permutation[pivot]);
        }

        double norm = std::sqrt(pivot_norm2);
        const double alpha = a[k][k] >= 0.0 ? -norm : norm;
        std::vector<double> reflector(samples.size() - k);
        for (std::size_t i = k; i < samples.size(); ++i)
            reflector[i-k] = a[i][k];
        reflector[0] -= alpha;
        double reflector_norm2 = 0.0;
        for (const double value : reflector)
            reflector_norm2 += value * value;
        if (!(reflector_norm2 > 0.0) || !std::isfinite(reflector_norm2))
            throw std::invalid_argument(
                "wall quadratic reconstruction: unstable QR reflector");
        const double beta = 2.0 / reflector_norm2;

        for (std::size_t j = k; j < n; ++j) {
            double projection = 0.0;
            for (std::size_t i = k; i < samples.size(); ++i)
                projection += reflector[i-k] * a[i][j];
            projection *= beta;
            for (std::size_t i = k; i < samples.size(); ++i)
                a[i][j] -= projection * reflector[i-k];
        }
        double rhs_projection = 0.0;
        for (std::size_t i = k; i < samples.size(); ++i)
            rhs_projection += reflector[i-k] * b[i];
        rhs_projection *= beta;
        for (std::size_t i = k; i < samples.size(); ++i)
            b[i] -= rhs_projection * reflector[i-k];

        a[k][k] = alpha;
        for (std::size_t i = k + 1; i < samples.size(); ++i)
            a[i][k] = 0.0;
        diagonal[k] = std::abs(alpha);
        norm = diagonal[k];
        if (!(norm > 0.0) || !std::isfinite(norm))
            throw std::invalid_argument(
                "wall quadratic reconstruction: rank-deficient stencil");
    }

    const auto [min_diagonal, max_diagonal] = std::minmax_element(
        diagonal.begin(), diagonal.end());
    // This is a conservative diagonal-ratio condition proxy for the
    // pivoted triangular factor, not a claim of a full SVD condition number.
    if (!(*max_diagonal > 0.0) ||
        *min_diagonal <= 1e-10 * *max_diagonal)
        throw std::invalid_argument(
            "wall quadratic reconstruction: ill-conditioned stencil");

    std::array<double, n> permuted_solution{};
    for (std::size_t reverse = n; reverse-- > 0;) {
        double value = b[reverse];
        for (std::size_t j = reverse + 1; j < n; ++j)
            value -= a[reverse][j] * permuted_solution[j];
        permuted_solution[reverse] = value / a[reverse][reverse];
        if (!std::isfinite(permuted_solution[reverse]))
            throw std::runtime_error(
                "wall quadratic reconstruction: non-finite fit coefficient");
    }
    std::array<double, n> result{};
    for (std::size_t j = 0; j < n; ++j)
        result[permutation[j]] = permuted_solution[j];
    return result;
}
} // namespace detail

inline WallFaceReconstruction reconstruct_wall_face_quadratic(
    const std::vector<WallQuadraticSample>& samples,
    const cfdx::core::Vec3& face_center,
    const cfdx::core::Vec3& area_vector_fluid_outward,
    double dynamic_viscosity)
{
    using cfdx::core::Vec3;
    const auto finite = [](const Vec3& v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    if (samples.size() < 10 || !finite(face_center) ||
        !finite(area_vector_fluid_outward) ||
        !(dynamic_viscosity >= 0.0) || !std::isfinite(dynamic_viscosity))
        throw std::invalid_argument("wall quadratic reconstruction: invalid input");
    double scale = 0.0;
    for (const auto& s : samples) {
        if (!finite(s.position) || !finite(s.velocity) || !std::isfinite(s.pressure))
            throw std::invalid_argument("wall quadratic reconstruction: non-finite sample");
        scale = std::max(scale, (s.position - face_center).mag());
    }
    const double area = area_vector_fluid_outward.mag();
    if (!(scale > 0.0) || !(area > 0.0) || !std::isfinite(scale) || !std::isfinite(area))
        throw std::invalid_argument("wall quadratic reconstruction: degenerate geometry");

    const auto p = detail::solve_wall_quadratic_fit(samples, face_center, scale, 0);
    std::array<std::array<double, 10>, 3> u{};
    for (int component = 0; component < 3; ++component)
        u[static_cast<std::size_t>(component)] =
            detail::solve_wall_quadratic_fit(samples, face_center, scale, component + 1);

    // Coefficients are in coordinates normalized by scale. At the expansion
    // point (the face centre), the linear coefficients are the exact gradient.
    WallFaceReconstruction result;
    result.pressure = p[0];
    const Vec3 grad_p{p[1]/scale, p[2]/scale, p[3]/scale};
    (void)grad_p; // The fitted polynomial is centred at the face itself.
    for (std::size_t i = 0; i < 3; ++i)
        result.velocity_gradient[i] = Vec3{
            u[i][1]/scale, u[i][2]/scale, u[i][3]/scale};
    const Vec3 n = area_vector_fluid_outward * (1.0/area);
    const auto& g = result.velocity_gradient;
    const Vec3 tau_n{
        dynamic_viscosity*(2.0*g[0].x*n.x + (g[0].y+g[1].x)*n.y + (g[0].z+g[2].x)*n.z),
        dynamic_viscosity*((g[1].x+g[0].y)*n.x + 2.0*g[1].y*n.y + (g[1].z+g[2].y)*n.z),
        dynamic_viscosity*((g[2].x+g[0].z)*n.x + (g[2].y+g[1].z)*n.y + 2.0*g[2].z*n.z)};
    result.force_per_area = n*result.pressure - tau_n;
    if (!std::isfinite(result.pressure) || !finite(result.force_per_area))
        throw std::runtime_error("wall quadratic reconstruction: non-finite result");
    return result;
}

} // namespace cfdx::physics
