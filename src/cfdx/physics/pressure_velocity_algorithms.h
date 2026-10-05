#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/flux.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

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

enum class CoupledSchurModel { BlockLocal, PCD, LSC, BFBT, SIMPLE, SIMPLEC };

inline const char* to_string(CoupledSchurModel model) noexcept {
    switch (model) {
        case CoupledSchurModel::BlockLocal: return "block_local";
        case CoupledSchurModel::PCD: return "pcd";
        case CoupledSchurModel::LSC: return "lsc";
        case CoupledSchurModel::BFBT: return "bfbt";
        case CoupledSchurModel::SIMPLE: return "simple";
        case CoupledSchurModel::SIMPLEC: return "simplec";
    }
    return "unknown";
}

inline bool parse_coupled_schur_model(const std::string& name, CoupledSchurModel& out) {
    if (name == "block_local") { out = CoupledSchurModel::BlockLocal; return true; }
    if (name == "pcd") { out = CoupledSchurModel::PCD; return true; }
    if (name == "lsc") { out = CoupledSchurModel::LSC; return true; }
    if (name == "bfbt") { out = CoupledSchurModel::BFBT; return true; }
    if (name == "simple") { out = CoupledSchurModel::SIMPLE; return true; }
    if (name == "simplec") { out = CoupledSchurModel::SIMPLEC; return true; }
    return false;
}

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
    // Explicit pressure-Schur approximation for the coupled production path.
    // Unsupported models must be rejected; there is no implicit fallback.
    CoupledSchurModel schur_model = CoupledSchurModel::PCD;
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
    // Small additive increases make the controller conservative near a stable
    // fixed point; degradation remains twice as fast so overshoot is recovered
    // more quickly than it is introduced.
    double increase_step = 0.005;
    double decrease_step = 0.02;
    // A severe residual jump is treated as an overshoot event rather than
    // waiting for the normal multi-window hysteresis. This guard is checked
    // on every completed nonlinear iteration.
    double severe_degradation_ratio = 2.0;
    // Number of adaptation windows during which increases are frozen after a
    // severe-degradation rollback.
    std::size_t recovery_cooldown_windows = 2;
    // Adapt against a short history window rather than consecutive noisy
    // iterations. This prevents the controller from reacting to the
    // pressure/velocity coupling oscillation of SIMPLE itself.
    std::size_t adaptation_window = 4;
    // Number of consecutive window decisions required before changing a
    // relaxation factor. This is deliberately explicit: the controller is
    // hysteretic and must not react to a single SIMPLE oscillation.
    std::size_t required_trend_windows = 2;
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
        !(c.severe_degradation_ratio > 1.0 + c.degradation_threshold) ||
        c.recovery_cooldown_windows == 0 ||
        c.adaptation_window == 0 || c.required_trend_windows == 0)
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

inline double adapt_relaxation_factor_windowed(
    double alpha,
    double reference_metric,
    double current_metric,
    double min_alpha,
    double max_alpha,
    const AdaptiveRelaxationControls& controls,
    std::size_t& improvement_streak,
    std::size_t& degradation_streak)
{
    validate_adaptive_relaxation_controls(controls);
    if (!std::isfinite(alpha) || !std::isfinite(reference_metric) ||
        !std::isfinite(current_metric) || reference_metric <= 0.0 ||
        current_metric < 0.0)
        throw std::invalid_argument("invalid windowed adaptive relaxation state");

    const double safe_alpha = std::clamp(alpha, min_alpha, max_alpha);
    if (!controls.enabled)
        return safe_alpha;

    const double ratio = current_metric / reference_metric;
    if (ratio <= 1.0 - controls.improvement_threshold) {
        ++improvement_streak;
        degradation_streak = 0;
    } else if (ratio >= 1.0 + controls.degradation_threshold) {
        ++degradation_streak;
        improvement_streak = 0;
    } else {
        improvement_streak = 0;
        degradation_streak = 0;
    }

    if (improvement_streak >= controls.required_trend_windows) {
        improvement_streak = 0;
        return std::min(max_alpha, safe_alpha + controls.increase_step);
    }
    if (degradation_streak >= controls.required_trend_windows) {
        degradation_streak = 0;
        return std::max(min_alpha, safe_alpha - controls.decrease_step);
    }
    return safe_alpha;
}

inline bool adaptive_relaxation_severe_degradation(
    double previous_metric,
    double current_metric,
    const AdaptiveRelaxationControls& controls)
{
    validate_adaptive_relaxation_controls(controls);
    if (!std::isfinite(previous_metric) || !std::isfinite(current_metric) ||
        previous_metric <= 0.0 || current_metric < 0.0)
        throw std::invalid_argument("invalid adaptive relaxation guard state");
    return current_metric / previous_metric >= controls.severe_degradation_ratio;
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
