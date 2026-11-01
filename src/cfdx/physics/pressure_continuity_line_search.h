#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

// Experimental, opt-in globalisation of a segregated pressure update.
// Evaluations MUST use the production conservative face flux, including all
// boundary conditions and pressure-response terms, not a matrix surrogate.
struct PressureContinuityLineSearchControls {
    bool enabled = false;
    std::size_t max_trials = 5;
    double reduction = 0.5;
    double minimum_relative_improvement = 1e-8;
    double absolute_floor = 1e-13;
};

struct PressureContinuityLineSearchResult {
    double alpha = 0.0;
    double baseline_norm = 0.0;
    double accepted_norm = 0.0;
    std::size_t trials = 0;
};

inline void validate_pressure_continuity_line_search(
    const PressureContinuityLineSearchControls& c)
{
    if (c.max_trials == 0 || c.max_trials > 32 ||
        !(c.reduction > 0.0 && c.reduction < 1.0) ||
        !(c.minimum_relative_improvement >= 0.0 &&
          c.minimum_relative_improvement < 1.0) ||
        !(c.absolute_floor > 0.0) ||
        !std::isfinite(c.reduction) ||
        !std::isfinite(c.minimum_relative_improvement) ||
        !std::isfinite(c.absolute_floor))
        throw std::invalid_argument("pressure continuity line search: invalid controls");
}

// score(alpha) is a pure evaluation: this routine does not commit the state.
// A failed selection is explicit, never silently accepted as converged.
template <class Score>
PressureContinuityLineSearchResult select_pressure_continuity_step(
    double nominal_alpha, const PressureContinuityLineSearchControls& c,
    Score&& score)
{
    validate_pressure_continuity_line_search(c);
    if (!(nominal_alpha > 0.0 && nominal_alpha <= 1.0) ||
        !std::isfinite(nominal_alpha))
        throw std::invalid_argument("pressure continuity line search: invalid alpha");
    const double baseline = score(0.0);
    if (!(baseline >= 0.0) || !std::isfinite(baseline))
        throw std::runtime_error("pressure continuity line search: invalid baseline");
    const double required = baseline * (1.0 - c.minimum_relative_improvement);
    double alpha = nominal_alpha;
    for (std::size_t i = 0; i < c.max_trials; ++i) {
        const double value = score(alpha);
        if (std::isfinite(value) && value >= 0.0 &&
            (baseline <= c.absolute_floor || value <= c.absolute_floor ||
             value < required))
            return {alpha, baseline, value, i + 1};
        alpha *= c.reduction;
    }
    throw std::runtime_error("pressure continuity line search: no decreasing step");
}

}  // namespace cfdx::physics
