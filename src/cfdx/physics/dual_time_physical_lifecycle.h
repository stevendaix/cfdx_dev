#pragma once

#include "cfdx/physics/dual_time_field_integrator.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

/**
 * Transactional physical-time lifecycle around DualTimeFieldIntegrator.
 *
 * The wrapped field integrator owns the pseudo-time nonlinear solve. This
 * layer owns physical-time acceptance/rejection, deterministic dt retry,
 * transient callback orchestration, and restartable accepted history.
 *
 * A rejected physical step never modifies accepted history.
 */
struct DualTimePhysicalRetryControls {
    int max_retries = 8;
    double retry_shrink = 0.5;
};

struct DualTimePhysicalCheckpoint {
    DualTimeField current;
    DualTimeField previous;
    bool has_previous = false;
    double dt_previous = 0.0;
    double physical_time = 0.0;
    std::size_t step = 0;
};

struct DualTimePhysicalStepReport {
    DualTimeFieldStepReport dual_time;
    double physical_time_before = 0.0;
    double physical_time_after = 0.0;
    double dt_requested = 0.0;
    double dt_accepted = 0.0;
    int retries = 0;
};

struct DualTimePhysicalLifecycle {
    DualTimeFieldIntegrator integrator;
    double physical_time = 0.0;
    std::size_t step = 0;

    using StepCallback = std::function<void(
        double physical_time,
        double dt,
        const DualTimeField& state)>;

    DualTimePhysicalLifecycle() = default;

    explicit DualTimePhysicalLifecycle(const DualTimeField& initial)
        : integrator(initial) {}

    void initialize(const DualTimeField& initial, double time = 0.0) {
        if (!std::isfinite(time))
            throw std::invalid_argument("dual-time physical time must be finite");
        integrator.initialize(initial);
        physical_time = time;
        step = 0;
    }

    bool initialized() const { return integrator.initialized(); }

    DualTimePhysicalCheckpoint checkpoint() const {
        if (!initialized())
            throw std::logic_error("cannot checkpoint an uninitialized dual-time lifecycle");
        return {
            integrator.history.current,
            integrator.history.previous,
            integrator.history.has_previous,
            integrator.history.dt_previous,
            physical_time,
            step
        };
    }

    void restore(const DualTimePhysicalCheckpoint& checkpoint) {
        if (checkpoint.current.size() == 0 ||
            checkpoint.current.dimension() == 0)
            throw std::invalid_argument("dual-time checkpoint has an empty current field");
        if (checkpoint.previous.size() != checkpoint.current.size() ||
            checkpoint.previous.dimension() != checkpoint.current.dimension())
            throw std::invalid_argument("dual-time checkpoint history shape mismatch");
        if (!std::isfinite(checkpoint.physical_time))
            throw std::invalid_argument("dual-time checkpoint has invalid physical time");
        if (checkpoint.has_previous &&
            (!(checkpoint.dt_previous > 0.0) ||
             !std::isfinite(checkpoint.dt_previous)))
            throw std::invalid_argument("dual-time checkpoint has invalid previous dt");

        integrator.history.current = checkpoint.current;
        integrator.history.previous = checkpoint.previous;
        integrator.history.has_previous = checkpoint.has_previous;
        integrator.history.dt_previous = checkpoint.dt_previous;
        physical_time = checkpoint.physical_time;
        step = checkpoint.step;
    }

    DualTimePhysicalStepReport advance(
        double dt_requested,
        const DualTimeFieldSolveControls& controls,
        const DualTimePhysicalRetryControls& retry_controls,
        const DualTimeFieldRhs& rhs,
        StepCallback before_attempt = {},
        StepCallback after_accept = {})
    {
        if (!(dt_requested > 0.0) || !std::isfinite(dt_requested))
            throw std::invalid_argument("requested physical dt must be finite and positive");
        if (retry_controls.max_retries < 0 ||
            !(retry_controls.retry_shrink > 0.0 &&
              retry_controls.retry_shrink < 1.0) ||
            !std::isfinite(retry_controls.retry_shrink))
            throw std::invalid_argument("invalid physical retry controls");

        const double time_before = physical_time;
        double dt = std::clamp(
            dt_requested, controls.physical_dt_min, controls.physical_dt_max);

        for (int retry = 0; retry <= retry_controls.max_retries; ++retry) {
            const DualTimeField accepted_before = integrator.history.current;
            if (before_attempt)
                before_attempt(physical_time, dt, accepted_before);

            try {
                const auto report = integrator.advance(dt, controls, rhs);

                const double time_after = physical_time + dt;
                if (!std::isfinite(time_after))
                    throw std::runtime_error(
                        "dual-time physical step produced a non-finite time");

                ++step;
                physical_time = time_after;

                if (after_accept)
                    after_accept(physical_time, dt, report.state);

                return {
                    report,
                    time_before,
                    physical_time,
                    dt_requested,
                    dt,
                    retry
                };
            } catch (const DualTimeConvergenceFailure&) {
                // DualTimeFieldIntegrator commits history only after all
                // acceptance gates pass. Keep this invariant explicit here:
                // a rejected attempt is transactional at physical time.
                integrator.history.current = accepted_before;
                // The embedded 2nd-order temporal estimator assumes the stored
                // history is consistent with the step that rejected. After a
                // shrink, the old history (spacing dt) is inconsistent with the
                // retried dt (spacing dt*shrink), so BDF2 would never clear
                // the temporal gate. Retrying therefore RESTARTS the temporal
                // history from the accepted state (self-start with the
                // first-order backward-Euler scheme): the embedded BE-vs-BE
                // estimate is exact and the first retry accepts deterministically.
                integrator.history.has_previous = false;
                integrator.history.previous = accepted_before;
                if (retry == retry_controls.max_retries)
                    throw;
                dt = std::max(
                    controls.physical_dt_min,
                    dt * retry_controls.retry_shrink);
            }
        }

        throw std::logic_error("unreachable dual-time physical retry state");
    }
};

} // namespace cfdx::physics
