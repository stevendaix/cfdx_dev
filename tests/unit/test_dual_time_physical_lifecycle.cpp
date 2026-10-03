#include "cfdx/physics/dual_time_physical_lifecycle.h"
#include <cmath>
#include <iostream>

using namespace cfdx::physics;

static DualTimeFieldRhs decay_rhs() {
    return [](const DualTimeField& state, DualTimeField& out) {
        out(0) = -state(0);
    };
}

int main() {
    DualTimeField initial(1, "u", "1", 1);
    initial(0) = 1.0;

    DualTimePhysicalLifecycle lifecycle(initial);

    DualTimeFieldSolveControls controls;
    controls.pseudo_time.max_pseudo_iterations = 2000;
    controls.pseudo_time.pseudo_dt_initial = 1.0;
    controls.pseudo_time.pseudo_dt_min = 1.0e-10;
    controls.pseudo_time.pseudo_dt_max = 1.0e6;
    controls.pseudo_time.pseudo_dt_growth = 1.5;
    controls.pseudo_time.pseudo_dt_shrink = 0.5;
    controls.physical_dt_min = 1.0e-8;
    controls.physical_dt_max = 0.2;
    // Tight temporal tolerances make the embedded BE/BDF2 estimator reject a
    // too-large physical step deterministically: the retry loop must shrink dt
    // until the temporal gate passes. This deterministically exercises the
    // transactional retry path (as opposed to relying on a fragile initial
    // acceptance).
    controls.temporal_absolute_tolerance = 1.0e-9;
    controls.temporal_relative_tolerance = 1.0e-6;

    DualTimePhysicalRetryControls retry;
    retry.max_retries = 20;
    retry.retry_shrink = 0.5;

    int before_calls = 0;
    int after_calls = 0;

    // First physical step: BACKWARD_EULER/BACKWARD_EULER estimator gives a zero
    // embedded error, so it is accepted without a retry and the physical time
    // and step counters advance atomically.
    const auto first = lifecycle.advance(
        0.1, controls, retry, decay_rhs(),
        [&](double, double, const DualTimeField&) { ++before_calls; },
        [&](double, double, const DualTimeField&) { ++after_calls; });

    if (first.retries != 0 || first.dual_time.state(0) <= 0.0)
        return 1;
    if (lifecycle.step != 1 || before_calls != 1 || after_calls != 1)
        return 2;
    if (std::abs(lifecycle.physical_time - first.dt_accepted) > 1e-14)
        return 3;
    const double value_after_step1 = first.dual_time.state(0);

    const auto cp = lifecycle.checkpoint();

    // Second step at the same requested dt: the 2nd-order (BDF2) vs backward
    // Euler embedded estimator now sees a real temporal error at dt=0.1 and
    // rejects it; the retry shrinks dt until the gate passes.
    const auto second = lifecycle.advance(0.1, controls, retry, decay_rhs());
    if (second.retries <= 0)
        return 4;
    if (!(second.dt_accepted < 0.1))
        return 5;
    if (!(second.dual_time.state(0) < value_after_step1))
        return 6;
    if (lifecycle.step != 2)
        return 7;
    if (std::abs(lifecycle.physical_time -
                 (cp.physical_time + second.dt_accepted)) > 1e-14)
        return 8;

    // Deterministic restart from the accepted checkpoint: two lifecycles
    // restored from the same checkpoint, advanced with the SAME accepted dt,
    // must produce identical state and identical temporal history semantics.
    DualTimePhysicalLifecycle restarted;
    restarted.restore(cp);
    DualTimePhysicalLifecycle reference;
    reference.restore(cp);

    const auto a = restarted.advance(second.dt_accepted, controls, retry, decay_rhs());
    const auto b = reference.advance(second.dt_accepted, controls, retry, decay_rhs());

    if (restarted.step != cp.step + 1 || reference.step != cp.step + 1)
        return 9;
    if (std::abs(restarted.physical_time - reference.physical_time) > 1e-14)
        return 10;
    if (std::abs(a.dual_time.state(0) - b.dual_time.state(0)) > 1e-13)
        return 11;
    if (restarted.integrator.history.has_previous !=
        reference.integrator.history.has_previous)
        return 12;

    std::cout << "dual-time physical lifecycle PASS"
              << " retries=" << second.retries
              << " dt_accepted=" << second.dt_accepted
              << " v1=" << value_after_step1
              << " v2=" << a.dual_time.state(0) << "\n";
    return 0;
}