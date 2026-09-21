#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct AdaptiveCflControls {
    double target_cfl = 1.0;
    double min_dt = 1e-12;
    double max_dt = 1e12;
    double growth_limit = 1.25;
    double shrink_limit = 0.5;
};

inline void validate_adaptive_cfl_controls(const AdaptiveCflControls& c)
{
    if (!std::isfinite(c.target_cfl) || !(c.target_cfl > 0.0) ||
        !std::isfinite(c.min_dt) || !(c.min_dt > 0.0) ||
        !std::isfinite(c.max_dt) || c.max_dt < c.min_dt ||
        !std::isfinite(c.growth_limit) || !(c.growth_limit >= 1.0) ||
        !std::isfinite(c.shrink_limit) || !(c.shrink_limit > 0.0) ||
        c.shrink_limit > 1.0 ||
        c.shrink_limit > c.growth_limit) {
        throw std::invalid_argument("invalid adaptive CFL controls");
    }
}

inline double adaptive_time_step(
    double dt, double measured_cfl, const AdaptiveCflControls& c = {})
{
    validate_adaptive_cfl_controls(c);
    if (!std::isfinite(dt) || !(dt > 0.0) ||
        !std::isfinite(measured_cfl) || !(measured_cfl > 0.0)) {
        throw std::invalid_argument("invalid CFL state");
    }

    double factor = std::sqrt(c.target_cfl / measured_cfl);
    factor = std::clamp(factor, c.shrink_limit, c.growth_limit);
    return std::clamp(dt * factor, c.min_dt, c.max_dt);
}

inline double pseudo_transient_cfl(
    std::size_t iteration,
    double cfl0 = 0.5,
    double growth = 1.2,
    double cfl_max = 100.0)
{
    if (!std::isfinite(cfl0) || !(cfl0 > 0.0) ||
        !std::isfinite(growth) || !(growth > 1.0) ||
        !std::isfinite(cfl_max) || cfl_max < cfl0) {
        throw std::invalid_argument("invalid pseudo-CFL controls");
    }
    return std::min(cfl_max,
                    cfl0 * std::pow(growth, static_cast<double>(iteration)));
}

}  // namespace cfdx::physics
