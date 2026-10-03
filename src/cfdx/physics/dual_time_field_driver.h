#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/dual_time_adaptive.h"
#include <cmath>
#include <functional>
#include <stdexcept>

namespace cfdx::physics {

using ScalarField = cfdx::core::Field<double, cfdx::core::Location::CELL>;

struct DualTimeFieldHistory {
    ScalarField current;
    ScalarField previous;
    bool has_previous = false;
    double dt_previous = 0.0;

    DualTimeFieldHistory() = default;

    explicit DualTimeFieldHistory(const ScalarField& initial)
        : current(initial),
          previous(initial),
          has_previous(false),
          dt_previous(0.0) {}

    void initialize(const ScalarField& initial) {
        current = initial;
        previous = initial;
        has_previous = false;
        dt_previous = 0.0;
    }

    void accept(const ScalarField& next, double dt) {
        if (!(dt > 0.0) || !std::isfinite(dt))
            throw std::invalid_argument("dual-time field history requires positive finite dt");
        if (next.size() != current.size() || next.dimension() != current.dimension())
            throw std::invalid_argument("dual-time field history shape mismatch");
        previous = current;
        current = next;
        has_previous = true;
        dt_previous = dt;
    }
};

struct DualTimeFieldStepResult {
    ScalarField state;
    DualTimeStepAcceptance acceptance;
    DualTimeTemporalErrorEstimate temporal_error;
    double dt_used = 0.0;
};

using DualTimeFieldSolve = std::function<ScalarField(
    const ScalarField& initial,
    const ScalarField& physical_previous,
    double dt)>;

inline DualTimeFieldStepResult accept_dual_time_field_step(
    const ScalarField& accepted_state,
    const ScalarField& lower_order_state,
    const ScalarField& dual_time_state,
    double dt,
    const DualTimeStepAcceptance& gates,
    double absolute_tolerance,
    double relative_tolerance)
{
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("dual-time field step requires positive finite dt");
    if (accepted_state.size() != lower_order_state.size() ||
        accepted_state.size() != dual_time_state.size() ||
        accepted_state.dimension() != lower_order_state.dimension() ||
        accepted_state.dimension() != dual_time_state.dimension())
        throw std::invalid_argument("dual-time field state shape mismatch");

    std::vector<double> high;
    std::vector<double> low;
    high.reserve(accepted_state.size() * accepted_state.dimension());
    low.reserve(accepted_state.size() * accepted_state.dimension());
    for (std::size_t d = 0; d < accepted_state.dimension(); ++d) {
        const double* hp = accepted_state.component_data(d);
        const double* lp = lower_order_state.component_data(d);
        for (std::size_t i = 0; i < accepted_state.size(); ++i) {
            high.push_back(hp[i]);
            low.push_back(lp[i]);
        }
    }

    const auto error = estimate_dual_time_temporal_error(
        high, low, absolute_tolerance, relative_tolerance);
    DualTimeFieldStepResult result{dual_time_state, gates, error, dt};
    if (!gates.accepted())
        throw DualTimeConvergenceFailure(
            "dual-time physical step rejected by temporal/nonlinear/admissibility gate");
    result.state = accepted_state;
    return result;
}

} // namespace cfdx::physics
