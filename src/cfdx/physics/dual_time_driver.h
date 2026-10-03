#pragma once

#include "cfdx/physics/dual_time_stepping.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <utility>

namespace cfdx::physics {

using DualTimeUpdate = std::function<std::vector<double>(
    const std::vector<double>&,
    const std::vector<double>&,
    double)>;

inline std::pair<std::vector<double>, DualTimeStepReport> run_dual_time_step(
    const std::vector<double>& initial_state,
    const std::vector<double>& state_n,
    const std::vector<double>& state_nm1,
    double dt,
    const DualTimeStepControls& controls,
    const std::function<std::vector<double>(const std::vector<double>&)>& rhs,
    const DualTimeUpdate& update)
{
    validate_dual_time_controls(controls);
    if (initial_state.size() != state_n.size())
        throw std::invalid_argument("dual-time initial/history size mismatch");
    if (controls.physical_scheme == DualTimePhysicalScheme::BDF2 &&
        initial_state.size() != state_nm1.size())
        throw std::invalid_argument("dual-time BDF2 initial/history size mismatch");

    std::vector<double> state = initial_state;
    double pseudo_dt = controls.pseudo_dt_initial;
    DualTimeStepReport report;
    double previous_residual = std::numeric_limits<double>::infinity();

    for (std::size_t iteration = 1; iteration <= controls.max_pseudo_iterations; ++iteration) {
        const auto residual = dual_time_physical_residual(
            state, state_n, state_nm1, dt, controls.physical_scheme, rhs);
        const double norm = dual_time_norm(residual);
        if (iteration == 1) report.initial_residual = norm;
        const double relative = report.initial_residual > 0.0
            ? norm / report.initial_residual : 0.0;
        report.history.push_back({iteration, norm, relative, pseudo_dt});
        if (norm <= controls.absolute_tolerance ||
            relative <= controls.relative_tolerance) {
            report.converged = true;
            report.pseudo_iterations = iteration;
            report.final_residual = norm;
            return {state, report};
        }

        const auto candidate = update(state, residual, pseudo_dt);
        if (candidate.size() != state.size())
            throw std::runtime_error("dual-time update changed state size");
        for (double value : candidate)
            if (!std::isfinite(value))
                throw std::runtime_error("dual-time update produced a non-finite state");

        const double next_residual = dual_time_norm(
            dual_time_physical_residual(candidate, state_n, state_nm1, dt,
                                         controls.physical_scheme, rhs));
        state = candidate;
        if (std::isfinite(previous_residual)) {
            if (next_residual < previous_residual)
                pseudo_dt = std::min(controls.pseudo_dt_max,
                                     pseudo_dt * controls.pseudo_dt_growth);
            else
                pseudo_dt = std::max(controls.pseudo_dt_min,
                                     pseudo_dt * controls.pseudo_dt_shrink);
        }
        previous_residual = next_residual;
    }

    report.pseudo_iterations = controls.max_pseudo_iterations;
    report.final_residual = report.history.back().residual_norm;
    throw DualTimeConvergenceFailure(
        "dual-time physical step did not converge within the pseudo-time iteration budget");
}

} // namespace cfdx::physics
