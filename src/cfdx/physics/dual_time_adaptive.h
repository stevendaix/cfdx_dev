#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct DualTimeTemporalErrorEstimate {
    double absolute_norm = 0.0;
    double relative_error = 0.0;
    double normalized_error = 0.0;
};

inline DualTimeTemporalErrorEstimate estimate_dual_time_temporal_error(
    const std::vector<double>& high_order_state,
    const std::vector<double>& low_order_state,
    double absolute_tolerance,
    double relative_tolerance)
{
    if (high_order_state.size() != low_order_state.size())
        throw std::invalid_argument("dual-time temporal-error state size mismatch");
    if (absolute_tolerance < 0.0 || !std::isfinite(absolute_tolerance) ||
        relative_tolerance < 0.0 || !std::isfinite(relative_tolerance) ||
        (absolute_tolerance == 0.0 && relative_tolerance == 0.0))
        throw std::invalid_argument("invalid dual-time temporal-error tolerances");

    long double diff2 = 0.0L;
    long double state2 = 0.0L;
    for (std::size_t i = 0; i < high_order_state.size(); ++i) {
        const double d = high_order_state[i] - low_order_state[i];
        if (!std::isfinite(high_order_state[i]) || !std::isfinite(low_order_state[i]))
            throw std::invalid_argument("non-finite dual-time temporal-error state");
        diff2 += static_cast<long double>(d) * static_cast<long double>(d);
        state2 += static_cast<long double>(high_order_state[i]) *
                  static_cast<long double>(high_order_state[i]);
    }

    const double absolute_norm = std::sqrt(static_cast<double>(diff2));
    const double state_norm = std::sqrt(static_cast<double>(state2));
    const double scale = absolute_tolerance + relative_tolerance * state_norm;
    return {
        absolute_norm,
        state_norm > 0.0 ? absolute_norm / state_norm : absolute_norm,
        scale > 0.0 ? absolute_norm / scale : std::numeric_limits<double>::infinity()
    };
}

struct DualTimePhysicalStepProposal {
    double dt = 0.0;
    bool accept = false;
    double error_ratio = 0.0;
};

inline DualTimePhysicalStepProposal propose_dual_time_physical_step(
    double current_dt,
    double normalized_error,
    int method_order,
    double dt_min,
    double dt_max,
    double safety = 0.9)
{
    if (!(current_dt > 0.0) || !std::isfinite(current_dt) ||
        !(dt_min > 0.0) || !std::isfinite(dt_min) ||
        !(dt_max >= dt_min) || !std::isfinite(dt_max) ||
        method_order <= 0 ||
        !(safety > 0.0 && safety < 1.0) || !std::isfinite(safety) ||
        normalized_error < 0.0 || !std::isfinite(normalized_error))
        throw std::invalid_argument("invalid dual-time physical-step controller input");

    const bool accept = normalized_error <= 1.0;
    const double bounded_error = std::max(normalized_error, 1.0e-14);
    const double exponent = 1.0 / static_cast<double>(method_order + 1);
    double factor = safety * std::pow(1.0 / bounded_error, exponent);
    factor = std::clamp(factor, accept ? 1.0 : 0.1, 5.0);

    return {
        std::clamp(current_dt * factor, dt_min, dt_max),
        accept,
        normalized_error
    };
}

struct DualTimePseudoTimeDecision {
    double pseudo_dt = 0.0;
    double contraction = 0.0;
    bool contracted = false;
};

inline DualTimePseudoTimeDecision adapt_dual_time_pseudo_dt(
    double pseudo_dt,
    double previous_residual,
    double current_residual,
    double pseudo_dt_min,
    double pseudo_dt_max,
    double target_contraction = 0.3,
    double safety = 0.9)
{
    if (!(pseudo_dt > 0.0) || !std::isfinite(pseudo_dt) ||
        !(pseudo_dt_min > 0.0) || !std::isfinite(pseudo_dt_min) ||
        !(pseudo_dt_max >= pseudo_dt_min) || !std::isfinite(pseudo_dt_max) ||
        !(previous_residual > 0.0) || !std::isfinite(previous_residual) ||
        current_residual < 0.0 || !std::isfinite(current_residual) ||
        !(target_contraction > 0.0 && target_contraction < 1.0) ||
        !(safety > 0.0 && safety < 1.0))
        throw std::invalid_argument("invalid dual-time pseudo-time controller input");

    const double contraction = current_residual / previous_residual;
    const bool contracted = contraction < 1.0;
    const double response = target_contraction / std::max(contraction, 1.0e-14);
    const double factor = std::clamp(safety * std::sqrt(response), 0.5, 2.0);

    return {
        std::clamp(pseudo_dt * factor, pseudo_dt_min, pseudo_dt_max),
        contraction,
        contracted
    };
}

struct DualTimeStepAcceptance {
    bool temporal_error_ok = false;
    bool nonlinear_converged = false;
    bool physically_admissible = false;

    bool accepted() const {
        return temporal_error_ok && nonlinear_converged && physically_admissible;
    }
};

} // namespace cfdx::physics
