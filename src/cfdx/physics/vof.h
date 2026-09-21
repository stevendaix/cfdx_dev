#pragma once

#include "cfdx/core/field/field.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct VofControls {
    double surface_tension = 0.0;
    double contact_angle = 90.0;
    double compression = 0.0;
};

inline void validate_volume_fraction(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& alpha,
    double tolerance = 1e-12)
{
    if (alpha.dimension() != 1)
        throw std::invalid_argument("VOF volume fraction must be scalar");
    for (std::size_t c = 0; c < alpha.size(); ++c) {
        if (!std::isfinite(alpha(c)) ||
            alpha(c) < -tolerance || alpha(c) > 1.0 + tolerance)
            throw std::runtime_error("VOF volume fraction outside [0,1]");
    }
}

// Conservative first-order geometric transport in 1-D.
// The returned field preserves the integral of alpha for periodic boundaries
// and is bounded by construction. It is a reference kernel for later
// polyhedral face-flux reconstruction/PLIC implementations.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
advect_volume_fraction_1d(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& alpha,
    double velocity,
    double dt,
    double dx)
{
    using namespace cfdx::core;
    if (!(dt >= 0.0) || !(dx > 0.0) || !std::isfinite(velocity))
        throw std::invalid_argument("invalid VOF advection controls");
    validate_volume_fraction(alpha);
    const std::size_t n = alpha.size();
    Field<double, Location::CELL> next(n, alpha.name(), alpha.metadata().unit, 1);
    if (n == 0) return next;

    const double cfl = velocity * dt / dx;
    if (std::abs(cfl) > 1.0 + 1e-14)
        throw std::invalid_argument("VOF reference advection requires CFL <= 1");

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t left = (i + n - 1) % n;
        const std::size_t right = (i + 1) % n;
        double value;
        if (cfl >= 0.0)
            value = alpha(i) - cfl * (alpha(i) - alpha(left));
        else
            value = alpha(i) - cfl * (alpha(right) - alpha(i));
        next(i) = std::clamp(value, 0.0, 1.0);
    }
    return next;
}

inline double interface_curvature_from_normals(
    const std::vector<double>& normal_x,
    const std::vector<double>& normal_y,
    std::size_t i,
    double dx)
{
    if (dx <= 0.0 || i == 0 || i + 1 >= normal_x.size() ||
        normal_x.size() != normal_y.size())
        throw std::invalid_argument("invalid normal stencil");
    return -((normal_x[i + 1] - normal_x[i - 1]) +
             (normal_y[i + 1] - normal_y[i - 1])) / (2.0 * dx);
}

inline double continuum_surface_force(double sigma, double curvature, double grad_alpha)
{
    if (!(sigma >= 0.0) || !std::isfinite(sigma) ||
        !std::isfinite(curvature) || !std::isfinite(grad_alpha))
        throw std::invalid_argument("invalid VOF surface-force parameters");
    return sigma * curvature * grad_alpha;
}

inline double contact_angle_wall_normal_component(double angle_degrees)
{
    if (!(angle_degrees >= 0.0 && angle_degrees <= 180.0))
        throw std::invalid_argument("contact angle must be in [0,180] degrees");
    const double radians = angle_degrees * 3.14159265358979323846 / 180.0;
    return std::cos(radians);
}

} // namespace cfdx::physics
