#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct InexactNewtonControls {
    double eta_min = 1e-4;
    double eta_max = 0.9;
    double gamma = 0.9;
    double alpha = 1.5;
    double residual_floor = 1e-14;
};

inline double eisenstat_walker_tolerance(
    double nonlinear_residual,
    double previous_residual,
    const InexactNewtonControls& c = {}) {
    if (!(nonlinear_residual >= 0.0) || !(previous_residual > 0.0) ||
        c.eta_min <= 0.0 || c.eta_max >= 1.0 || c.eta_min >= c.eta_max ||
        c.gamma <= 0.0 || c.gamma > 1.0 || c.alpha <= 0.0) {
        throw std::invalid_argument("invalid inexact Newton controls");
    }
    const double ratio = std::max(nonlinear_residual, c.residual_floor) /
                         std::max(previous_residual, c.residual_floor);
    const double eta = c.gamma * std::pow(ratio, c.alpha);
    return std::clamp(eta, c.eta_min, c.eta_max);
}

struct AdaptivePressureCorrectionControls {
    std::size_t min_correctors = 1;
    std::size_t max_correctors = 3;
    double high_continuity_residual = 1e-3;
    double low_continuity_residual = 1e-6;
};

inline std::size_t choose_pressure_correctors(
    double continuity_residual,
    const AdaptivePressureCorrectionControls& c = {}) {
    if (c.min_correctors == 0 || c.min_correctors > c.max_correctors ||
        !(c.high_continuity_residual > c.low_continuity_residual) ||
        c.low_continuity_residual < 0.0) {
        throw std::invalid_argument("invalid pressure-correction controls");
    }
    if (continuity_residual >= c.high_continuity_residual) return c.max_correctors;
    if (continuity_residual <= c.low_continuity_residual) return c.min_correctors;
    const double t = (continuity_residual - c.low_continuity_residual) /
                     (c.high_continuity_residual - c.low_continuity_residual);
    const auto span = c.max_correctors - c.min_correctors;
    return c.min_correctors + static_cast<std::size_t>(std::ceil(t * static_cast<double>(span)));
}

struct PhysicsStoppingControls {
    double nonlinear_tolerance = 1e-7;
    double linear_tolerance_min = 1e-10;
    double linear_tolerance_max = 1e-3;
};

inline double physics_based_linear_tolerance(
    double nonlinear_residual,
    const PhysicsStoppingControls& c = {}) {
    if (!(nonlinear_residual >= 0.0) || c.nonlinear_tolerance <= 0.0 ||
        c.linear_tolerance_min <= 0.0 ||
        c.linear_tolerance_min > c.linear_tolerance_max) {
        throw std::invalid_argument("invalid physics stopping controls");
    }
    const double target = std::sqrt(std::max(nonlinear_residual, 1e-30));
    return std::clamp(target, c.linear_tolerance_min, c.linear_tolerance_max);
}

} // namespace cfdx::physics
