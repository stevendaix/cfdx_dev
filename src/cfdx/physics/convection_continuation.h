#pragma once
// Optional conservative deferred high-order convection homotopy for steady
// SIMPLE. It is NOT a convergence shortcut: only the full-order operator
// after a prescribed hold period may satisfy the physical convergence gate.
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct ConvectionContinuationControls {
    bool enabled = false;
    std::size_t upwind_warmup_iterations = 40;
    std::size_t high_order_ramp_iterations = 80;
    std::size_t full_order_hold_iterations = 4;
};

inline std::size_t convection_full_order_iteration(
    const ConvectionContinuationControls& controls) {
    if (controls.high_order_ramp_iterations == 0 ||
        controls.full_order_hold_iterations < 2 ||
        controls.upwind_warmup_iterations >
            std::numeric_limits<std::size_t>::max() - controls.high_order_ramp_iterations)
        throw std::invalid_argument("convection continuation: invalid ramp or overflow");
    return controls.upwind_warmup_iterations + controls.high_order_ramp_iterations;
}

inline void validate_convection_continuation(
    const ConvectionContinuationControls& controls,
    std::size_t max_iterations) {
    if (!controls.enabled) return;
    const std::size_t full_at = convection_full_order_iteration(controls);
    if (full_at > max_iterations ||
        controls.full_order_hold_iterations > max_iterations - full_at + 1)
        throw std::invalid_argument(
            "convection continuation: insufficient budget for full-order certification");
}

inline double convection_high_order_fraction(
    const ConvectionContinuationControls& controls, std::size_t iteration) {
    if (!controls.enabled) return 1.0;
    const std::size_t full_at = convection_full_order_iteration(controls);
    if (iteration <= controls.upwind_warmup_iterations) return 0.0;
    if (iteration >= full_at) return 1.0;
    return static_cast<double>(iteration - controls.upwind_warmup_iterations) /
        static_cast<double>(controls.high_order_ramp_iterations);
}

// A growth or adaptive-degradation check may compare nonlinear residuals
// only when the *same discrete convection operator* is active. At a ramp
// boundary the newly evaluated residual includes a changed deferred-order
// correction; comparing it with a residual of the previous operator can
// incorrectly reject all trials before the ramp begins. This DOES NOT relax
// the final full-order physical convergence gates or permit nonfinite states.
inline bool convection_residuals_share_operator(
    const ConvectionContinuationControls& controls, std::size_t iteration) {
    if (iteration == 0)
        throw std::invalid_argument("convection operator comparison: invalid iteration");
    return iteration == 1 ||
        convection_high_order_fraction(controls, iteration) ==
        convection_high_order_fraction(controls, iteration - 1);
}

// The nonlinear monitor receives a dimensionless, RHS-normalized synthetic
// residual. A reference from an intermediate discretization cannot be used as
// the denominator for a full-order convergence certificate. Preserve the
// ordinary absolute criterion when continuation is disabled.
inline double convection_monitor_residual_threshold(
    const ConvectionContinuationControls& controls,
    double absolute_tolerance, double relative_tolerance) {
    if (!(absolute_tolerance > 0.0) || !(relative_tolerance > 0.0) ||
        !std::isfinite(absolute_tolerance) || !std::isfinite(relative_tolerance))
        throw std::invalid_argument("convection monitor: invalid tolerance");
    return controls.enabled ? relative_tolerance : absolute_tolerance;
}

inline bool convection_full_order_eligible(
    const ConvectionContinuationControls& controls, std::size_t iteration) {
    if (!controls.enabled) return true;
    const std::size_t full_at = convection_full_order_iteration(controls);
    return iteration >= full_at &&
        iteration - full_at + 1 >= controls.full_order_hold_iterations;
}

} // namespace cfdx::physics
