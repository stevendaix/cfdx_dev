#pragma once

#include "cfdx/core/numerics/temporal.h"
#include "cfdx/physics/dual_time_adaptive.h"
#include "cfdx/physics/dual_time_driver.h"
#include "cfdx/physics/dual_time_field_driver.h"
#include "cfdx/physics/dual_time_stepping.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

using DualTimeField = cfdx::core::Field<double, cfdx::core::Location::CELL>;
using DualTimeFieldRhs = std::function<void(
    const DualTimeField&,
    DualTimeField&)>;

struct DualTimeFieldSolveControls {
    DualTimeStepControls pseudo_time;
    double physical_dt_min = 1.0e-8;
    double physical_dt_max = 1.0;
    double physical_dt_growth = 1.5;
    double physical_dt_shrink = 0.5;
    double temporal_absolute_tolerance = 1.0e-6;
    double temporal_relative_tolerance = 1.0e-3;
    bool use_embedded_be_estimator = true;
    std::function<bool(const DualTimeField&)> physical_admissibility;
};

struct DualTimeFieldStepReport {
    DualTimeField state;
    DualTimeField lower_order_state;
    DualTimeStepReport nonlinear_report;
    DualTimeStepAcceptance acceptance;
    DualTimeTemporalErrorEstimate temporal_error;
    double dt_requested = 0.0;
    double dt_used = 0.0;
    double dt_proposed = 0.0;
};

struct DualTimeFieldIntegrator {
    DualTimeFieldHistory history;

    DualTimeFieldIntegrator() = default;
    explicit DualTimeFieldIntegrator(const DualTimeField& initial) : history(initial) {}

    void initialize(const DualTimeField& initial) {
        history.initialize(initial);
    }

    bool initialized() const {
        return history.current.size() != 0;
    }

private:
    static std::vector<double> flatten(const DualTimeField& field) {
        std::vector<double> values;
        values.reserve(field.size() * field.dimension());
        for (std::size_t d = 0; d < field.dimension(); ++d) {
            const double* data = field.component_data(d);
            for (std::size_t i = 0; i < field.size(); ++i) {
                if (!std::isfinite(data[i]))
                    throw std::invalid_argument("dual-time field contains a non-finite value");
                values.push_back(data[i]);
            }
        }
        return values;
    }

    static DualTimeField unflatten(
        const DualTimeField& prototype,
        const std::vector<double>& values,
        const std::string& name)
    {
        const std::size_t expected = prototype.size() * prototype.dimension();
        if (values.size() != expected)
            throw std::invalid_argument("dual-time field vector shape mismatch");

        DualTimeField result(
            prototype.size(), name, prototype.metadata().unit, prototype.dimension());
        std::size_t k = 0;
        for (std::size_t d = 0; d < prototype.dimension(); ++d) {
            double* data = result.component_data(d);
            for (std::size_t i = 0; i < prototype.size(); ++i)
                data[i] = values[k++];
        }
        return result;
    }

    static void validate_controls(const DualTimeFieldSolveControls& c) {
        validate_dual_time_controls(c.pseudo_time);
        if (!(c.physical_dt_min > 0.0) ||
            !(c.physical_dt_max >= c.physical_dt_min) ||
            !std::isfinite(c.physical_dt_min) ||
            !std::isfinite(c.physical_dt_max))
            throw std::invalid_argument("invalid physical time-step bounds");
        if (!(c.physical_dt_growth > 1.0) ||
            !std::isfinite(c.physical_dt_growth))
            throw std::invalid_argument("invalid physical time-step growth factor");
        if (!(c.physical_dt_shrink > 0.0 && c.physical_dt_shrink < 1.0) ||
            !std::isfinite(c.physical_dt_shrink))
            throw std::invalid_argument("invalid physical time-step shrink factor");
        if (!(c.temporal_absolute_tolerance >= 0.0) ||
            !(c.temporal_relative_tolerance >= 0.0) ||
            !std::isfinite(c.temporal_absolute_tolerance) ||
            !std::isfinite(c.temporal_relative_tolerance) ||
            (c.temporal_absolute_tolerance == 0.0 &&
             c.temporal_relative_tolerance == 0.0))
            throw std::invalid_argument("invalid temporal error tolerances");
    }

    std::pair<DualTimeField, DualTimeStepReport> solve_one(
        const DualTimeField& initial,
        const DualTimeField& state_n,
        const DualTimeField& state_nm1,
        double dt,
        const DualTimeFieldSolveControls& controls,
        const DualTimeFieldRhs& rhs) const
    {
        const auto initial_values = flatten(initial);
        const auto n_values = flatten(state_n);
        const auto nm1_values = flatten(state_nm1);

        auto rhs_vector = [&](const std::vector<double>& state) {
            const auto field = unflatten(initial, state, "dual_time_state");
            DualTimeField rhs_field(
                initial.size(), "dual_time_rhs",
                initial.metadata().unit + "/s", initial.dimension());
            rhs(field, rhs_field);
            return flatten(rhs_field);
        };

        auto update = [](
            const std::vector<double>& state,
            const std::vector<double>& residual,
            double pseudo_dt) {
            std::vector<double> next(state.size());
            for (std::size_t i = 0; i < state.size(); ++i) {
                next[i] = state[i] - pseudo_dt * residual[i];
                if (!std::isfinite(next[i]))
                    throw std::runtime_error("dual-time field update produced a non-finite value");
            }
            return next;
        };

        auto solved = run_dual_time_step(
            initial_values, n_values, nm1_values, dt,
            controls.pseudo_time, rhs_vector, update);
        return {
            unflatten(initial, solved.first, "dual_time_state"),
            std::move(solved.second)
        };
    }

public:
    DualTimeFieldStepReport advance(
        double dt_requested,
        const DualTimeFieldSolveControls& controls,
        const DualTimeFieldRhs& rhs)
    {
        validate_controls(controls);
        if (!initialized())
            throw std::logic_error("dual-time field integrator is not initialized");
        if (!rhs)
            throw std::invalid_argument("dual-time field RHS callback must be valid");
        if (!(dt_requested > 0.0) || !std::isfinite(dt_requested))
            throw std::invalid_argument("requested physical dt must be finite and positive");

        const double dt = std::clamp(
            dt_requested, controls.physical_dt_min, controls.physical_dt_max);

        DualTimeField state_nm1 = history.has_previous
            ? history.previous : history.current;

        DualTimeFieldSolveControls high_controls = controls;
        high_controls.pseudo_time.physical_scheme = history.has_previous
            ? controls.pseudo_time.physical_scheme
            : DualTimePhysicalScheme::BACKWARD_EULER;
        const auto high = solve_one(
            history.current, history.current, state_nm1, dt, high_controls, rhs);

        DualTimeStepControls low_controls = controls.pseudo_time;
        low_controls.physical_scheme = DualTimePhysicalScheme::BACKWARD_EULER;
        const auto low = solve_one(
            history.current, history.current, state_nm1, dt, 
            DualTimeFieldSolveControls{
                low_controls,
                controls.physical_dt_min,
                controls.physical_dt_max,
                controls.physical_dt_growth,
                controls.physical_dt_shrink,
                controls.temporal_absolute_tolerance,
                controls.temporal_relative_tolerance,
                controls.use_embedded_be_estimator,
                controls.physical_admissibility
            }, rhs);

        const auto high_values = flatten(high.first);
        const auto low_values = flatten(low.first);
        const auto error = estimate_dual_time_temporal_error(
            high_values, low_values,
            controls.temporal_absolute_tolerance,
            controls.temporal_relative_tolerance);

        const bool temporal_ok = !controls.use_embedded_be_estimator || error.accepted();
        const bool nonlinear_ok = high.second.converged;
        const bool physical_ok = controls.physical_admissibility
            ? controls.physical_admissibility(high.first)
            : true;
        if (!physical_ok)
            throw DualTimeConvergenceFailure("dual-time physical step rejected by admissibility callback");

        DualTimeStepAcceptance gates{temporal_ok, nonlinear_ok, physical_ok};
        const auto accepted = accept_dual_time_field_step(
            high.first, low.first, high.first, dt, gates,
            controls.temporal_absolute_tolerance,
            controls.temporal_relative_tolerance);

        const int method_order = high_controls.pseudo_time.physical_scheme ==
            DualTimePhysicalScheme::BDF2 ? 2 : 1;
        const auto proposal = propose_dual_time_physical_step(
            dt, error.normalized_error, method_order,
            controls.physical_dt_min, controls.physical_dt_max);
        const double dt_proposed = proposal.dt;

        history.accept(accepted.state, dt);

        return {
            accepted.state,
            low.first,
            high.second,
            accepted.acceptance,
            accepted.temporal_error,
            dt_requested,
            dt,
            dt_proposed
        };
    }
};

} // namespace cfdx::physics
