#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/flux.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

enum class PressureVelocityAlgorithm {
    SIMPLE,
    SIMPLEC,
    PISO,
    PIMPLE,
    // Fractional-step projection: momentum predictor followed by one or
    // more pressure projections. This is intentionally distinct from PISO:
    // each projection operates on the current conservative face flux.
    FRACTIONAL_STEP,
    // Fully coupled pressure-based solve. Momentum and continuity are
    // assembled in one block system with an implicit Rhie-Chow pressure block.
    COUPLED
};

struct CouplingControls {
    double alpha_u = 0.7;
    double alpha_p = 0.3;
    int n_outer_correctors = 1;
    int n_pressure_correctors = 2;
    int n_fractional_steps = 2;
    // Coupled linear solve controls. The coupled matrix is indefinite but the
    // Rhie-Chow pressure block makes the collocated formulation nonsingular
    // after pressure gauge fixing.
    std::size_t coupled_max_iterations = 2000;
    double coupled_linear_tolerance = 1e-10;
};

inline void validate_coupling_controls(const CouplingControls& c)
{
    if (!(c.alpha_u > 0.0 && c.alpha_u <= 1.0) ||
        !(c.alpha_p > 0.0 && c.alpha_p <= 1.0) ||
        c.n_outer_correctors < 1 || c.n_pressure_correctors < 1 ||
        c.n_fractional_steps < 1 || c.coupled_max_iterations == 0 ||
        !(c.coupled_linear_tolerance > 0.0))
        throw std::invalid_argument("pressure-velocity controls: invalid relaxation/corrector count");
}

struct AdaptiveRelaxationControls {
    bool enabled = false;
    double min_alpha_u = 0.2;
    double max_alpha_u = 0.9;
    double min_alpha_p = 0.1;
    double max_alpha_p = 0.5;
    double improvement_threshold = 0.05;
    double degradation_threshold = 0.10;
    // Additive updates avoid the large multiplicative jumps that can collapse
    // the controller to its floor and create a deterministic limit cycle.
    double increase_step = 0.01;
    double decrease_step = 0.02;
    // Adapt against a short history window rather than consecutive noisy
    // iterations. This prevents the controller from reacting to the
    // pressure/velocity coupling oscillation of SIMPLE itself.
    std::size_t adaptation_window = 4;
};

inline void validate_adaptive_relaxation_controls(const AdaptiveRelaxationControls& c)
{
    if (c.min_alpha_u <= 0.0 || c.max_alpha_u < c.min_alpha_u || c.max_alpha_u > 1.0 ||
        c.min_alpha_p <= 0.0 || c.max_alpha_p < c.min_alpha_p || c.max_alpha_p > 1.0 ||
        c.improvement_threshold < 0.0 || c.degradation_threshold < c.improvement_threshold ||
        !(c.increase_step > 0.0) || !(c.decrease_step > 0.0) ||
        c.increase_step > (c.max_alpha_u - c.min_alpha_u) ||
        c.decrease_step > (c.max_alpha_u - c.min_alpha_u) ||
        c.increase_step > (c.max_alpha_p - c.min_alpha_p) ||
        c.decrease_step > (c.max_alpha_p - c.min_alpha_p) ||
        c.adaptation_window == 0)
        throw std::invalid_argument("invalid adaptive relaxation controls");
}

inline double adapt_relaxation_factor(
    double alpha,
    double previous_residual,
    double residual,
    double min_alpha,
    double max_alpha,
    const AdaptiveRelaxationControls& controls)
{
    validate_adaptive_relaxation_controls(controls);
    if (!std::isfinite(alpha) || !std::isfinite(previous_residual) ||
        !std::isfinite(residual) || previous_residual < 0.0 || residual < 0.0)
        throw std::invalid_argument("invalid adaptive relaxation state");
    const double safe_alpha = std::clamp(alpha, min_alpha, max_alpha);
    if (!controls.enabled || previous_residual == 0.0)
        return safe_alpha;
    const double ratio = residual / previous_residual;
    if (ratio <= 1.0 - controls.improvement_threshold)
        return std::min(max_alpha, safe_alpha + controls.increase_step);
    if (ratio >= 1.0 + controls.degradation_threshold)
        return std::max(min_alpha, safe_alpha - controls.decrease_step);
    // Hysteresis band: retain the current relaxation and filter small residual noise.
    return safe_alpha;
}

inline double relaxed_value(double old_value, double computed_value, double alpha)
{
    if (!std::isfinite(alpha) || !std::isfinite(old_value) || !std::isfinite(computed_value) ||
        alpha <= 0.0 || alpha > 1.0)
        throw std::invalid_argument("relaxed_value: alpha must be in (0,1]");
    return old_value + alpha * (computed_value - old_value);
}

// Discrete coefficient convention used throughout the pressure-velocity path:
//
//   A_P   = integrated FV momentum diagonal
//   rAU   = 1/A_P
//   dAU   = V*rAU = V/A_P
//
// The volume is intentionally NOT folded into rAU. rAU belongs to the
// momentum algebra; dAU is the physical velocity response to a pressure
// gradient. Pressure correction, Rhie-Chow, and coupled Schur coefficients
// must use dAU (or an equivalent V*rAU expression) so they remain consistent
// with the cell velocity reconstruction and continuity flux. Keep this
// distinction explicit because 1/A_P and V/A_P are both valid quantities,
// but they belong to different discrete operators.

inline double piso_correction_gain(double diagonal, double neighbor_sum)
{
    if (!std::isfinite(diagonal) || !std::isfinite(neighbor_sum) || diagonal <= 0.0)
        throw std::invalid_argument("piso_correction_gain: diagonal must be positive");
    return 1.0 / std::max(diagonal - neighbor_sum, diagonal * 1e-12);
}


// The fractional-step pressure projection is deliberately explicit about its
// contract. It is a projection of the conservative face flux, not an alias for
// PISO. Keeping this helper here also makes the algorithm choice visible to
// callers without embedding policy in the finite-volume transport layer.
inline bool is_fractional_step_algorithm(PressureVelocityAlgorithm a)
{
    return a == PressureVelocityAlgorithm::FRACTIONAL_STEP;
}

inline bool is_coupled_algorithm(PressureVelocityAlgorithm a)
{
    return a == PressureVelocityAlgorithm::COUPLED;
}

}  // namespace cfdx::physics
