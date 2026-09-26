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

inline void validate_inexact_newton_controls(const InexactNewtonControls& c) {
    if (!std::isfinite(c.eta_min) || !std::isfinite(c.eta_max) ||
        !std::isfinite(c.gamma) || !std::isfinite(c.alpha) ||
        !std::isfinite(c.residual_floor) || c.eta_min <= 0.0 ||
        c.eta_max >= 1.0 || c.eta_min >= c.eta_max ||
        c.gamma <= 0.0 || c.gamma > 1.0 || c.alpha <= 0.0 ||
        c.residual_floor <= 0.0) {
        throw std::invalid_argument("invalid inexact Newton controls");
    }
}

inline double eisenstat_walker_tolerance(
    double nonlinear_residual,
    double previous_residual,
    const InexactNewtonControls& c = {}) {
    validate_inexact_newton_controls(c);
    if (!std::isfinite(nonlinear_residual) || nonlinear_residual < 0.0 ||
        !std::isfinite(previous_residual) || previous_residual <= 0.0)
        throw std::invalid_argument("invalid nonlinear residual history");
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

inline void validate_adaptive_pressure_correction_controls(
    const AdaptivePressureCorrectionControls& c) {
    if (c.min_correctors == 0 || c.min_correctors > c.max_correctors ||
        !std::isfinite(c.high_continuity_residual) ||
        !std::isfinite(c.low_continuity_residual) ||
        !(c.high_continuity_residual > c.low_continuity_residual) ||
        c.low_continuity_residual < 0.0) {
        throw std::invalid_argument("invalid pressure-correction controls");
    }
}

inline std::size_t choose_pressure_correctors(
    double continuity_residual,
    const AdaptivePressureCorrectionControls& c = {}) {
    validate_adaptive_pressure_correction_controls(c);
    if (!std::isfinite(continuity_residual) || continuity_residual < 0.0)
        throw std::invalid_argument("invalid continuity residual");
    if (continuity_residual >= c.high_continuity_residual) return c.max_correctors;
    if (continuity_residual <= c.low_continuity_residual) return c.min_correctors;
    const double t = (continuity_residual - c.low_continuity_residual) /
                     (c.high_continuity_residual - c.low_continuity_residual);
    const auto span = c.max_correctors - c.min_correctors;
    return c.min_correctors + static_cast<std::size_t>(std::ceil(t * static_cast<double>(span)));
}

struct ConvergenceAccelerationControls {
    bool adaptive_linear_tolerance = false;
    bool adaptive_pressure_correctors = false;
    // CFD outer iterations rarely benefit from solving every inner system to
    // machine accuracy. These tighter bounds than generic inexact Newton are
    // deliberately conservative for segregated pressure-velocity coupling.
    InexactNewtonControls linear_forcing{1e-6, 5e-2, 0.9, 1.5, 1e-14};
    AdaptivePressureCorrectionControls pressure_correctors;
};

inline void validate_convergence_acceleration_controls(
    const ConvergenceAccelerationControls& c) {
    validate_inexact_newton_controls(c.linear_forcing);
    validate_adaptive_pressure_correction_controls(c.pressure_correctors);
}

inline double nonlinear_forcing_tolerance(
    double nonlinear_residual,
    double previous_residual,
    const InexactNewtonControls& c) {
    const double eisenstat = eisenstat_walker_tolerance(
        nonlinear_residual, previous_residual, c);
    const double residual_target = std::sqrt(
        std::max(nonlinear_residual, c.residual_floor));
    return std::clamp(
        std::min(eisenstat, residual_target), c.eta_min, c.eta_max);
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
