#pragma once

// Experimental safeguard for a COMPLETE returned U,p,k,omega,nu_t state.
// This is an ACCEPTANCE policy for a nonlinear *trial*, NOT a physical
// convergence tolerance. It compares reassembled residuals from the SAME
// final convection operator and only after the SST state has been updated.
// The solver's nonlinear rollback must restore U,p,phi,p_low,SST together.

#include "cfdx/physics/global_nonlinear_residual_v8.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

struct GlobalNonlinearAcceptanceV13 {
    bool allowed = false;
    double largest_growth = 0.0;
    std::string failed_equation;
};

inline GlobalNonlinearAcceptanceV13 evaluate_global_nonlinear_trial_v13(
    const GlobalNonlinearResidualV8& previous,
    const GlobalNonlinearResidualV8& trial,
    const ConvergenceCriteria& physical_gates,
    double allowed_growth = 1.20,
    double noise_guard_multiplier = 10.0)
{
    if (!(allowed_growth >= 1.0 && allowed_growth <= 2.0) ||
        !std::isfinite(allowed_growth) ||
        !(noise_guard_multiplier >= 1.0 && noise_guard_multiplier <= 100.0) ||
        !std::isfinite(noise_guard_multiplier) ||
        !(physical_gates.relative_tolerance > 0.) ||
        !(physical_gates.continuity_tolerance > 0.) ||
        !std::isfinite(physical_gates.relative_tolerance) ||
        !std::isfinite(physical_gates.continuity_tolerance))
        throw std::invalid_argument("V13 joint trial: invalid acceptance controls");
    if (!previous.all_finite || !trial.all_finite ||
        !previous.sst || !trial.sst ||
        !previous.bounded || !trial.bounded ||
        !previous.sst_transport_committed || !trial.sst_transport_committed ||
        previous.algorithm != trial.algorithm ||
        previous.schur != trial.schur)
        throw std::invalid_argument("V13 joint trial: missing comparable accepted U,p,SST state");

    const std::array<double, 6> old_values = {{
        previous.momentum.max_active_component_relative, previous.continuity,
        previous.k, previous.omega, previous.eddy_viscosity_closure,
        previous.momentum.global_relative}};
    const std::array<double, 6> new_values = {{
        trial.momentum.max_active_component_relative, trial.continuity,
        trial.k, trial.omega, trial.eddy_viscosity_closure,
        trial.momentum.global_relative}};
    constexpr std::array<const char*, 6> names = {{
        "U_component", "p_continuity", "k", "omega", "nu_t_closure", "U_global"}};
    GlobalNonlinearAcceptanceV13 result;
    result.allowed = true;
    for (std::size_t i = 0; i < old_values.size(); ++i) {
        const double old_r = old_values[i], new_r = new_values[i];
        if (!std::isfinite(old_r) || !std::isfinite(new_r) || old_r < 0 || new_r < 0)
            throw std::invalid_argument("V13 joint trial: nonfinite/negative physical residual");
        const double floor = noise_guard_multiplier * (i == 1
            ? physical_gates.continuity_tolerance
            : physical_gates.relative_tolerance);
        const double base = std::max(old_r, floor);
        const double ratio = new_r / base;
        if (ratio > result.largest_growth) result.largest_growth = ratio;
        // Growth may be tolerated only within a RETRY envelope. The ORIGINAL
        // physical convergence checks are still evaluated separately, unchanged.
        if (new_r > allowed_growth * base && result.allowed) {
            result.allowed = false;
            result.failed_equation = names[i];
        }
    }
    return result;
}

} // namespace cfdx::physics
