#include "cfdx/cfdx/core/utils/convergence_monitor.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include "common/live_assert.hpp"

int main() {
    using namespace cfdx::core;

    ConvergenceMonitorControls controls;
    controls.residual.absolute_tolerance = 1e-8;
    controls.residual.relative_tolerance = 1e-4;
    controls.conservation_tolerance = 1e-8;
    controls.minimum_iterations = 2;
    controls.max_iterations = 20;
    controls.stagnation_window = 3;
    controls.stagnation_relative_improvement = 1e-3;
    controls.divergence_window = 2;
    controls.divergence_growth_factor = 5.0;
    controls.qoi_gates.push_back(QoIGate{"drag", 1e-8, 1e-4});

    ConvergenceMonitor monitor(controls);

    auto r1 = monitor.update({1, 1.0, 1e-2, {{"drag", 1.0}}});
    assert(r1.status == ConvergenceStatus::CONTINUE);
    auto r2 = monitor.update({2, 1e-2, 1e-3, {{"drag", 1.0001}}});
    assert(r2.status == ConvergenceStatus::CONTINUE);
    auto r3 = monitor.update({3, 1e-7, 1e-9, {{"drag", 1.000100001}}});
    assert(r3.status == ConvergenceStatus::CONVERGED);
    assert(r3.residual_gate);
    assert(r3.conservation_gate);
    assert(r3.qoi_gate);
    assert(monitor.deterministic_report().find("status=CONVERGED") == 0);

    monitor.reset();
    ConvergenceMonitorControls stagnation = controls;
    stagnation.qoi_gates.clear();
    stagnation.stagnation_window = 2;
    stagnation.stagnation_relative_improvement = 0.1;
    ConvergenceMonitor stagnant(stagnation);
    stagnant.update({1, 1.0, 1.0, {}});
    stagnant.update({2, 0.99, 1.0, {}});
    auto sr = stagnant.update({3, 0.989, 1.0, {}});
    assert(sr.status == ConvergenceStatus::STAGNATED);

    // Regression (#906): a materially rising residual is the divergence
    // detector's concern, not stagnation. The stagnation window must only
    // count flat-envelope samples (relative change below the threshold in
    // magnitude); counting rising samples made a healthy transient overshoot
    // of a synthetic max-residual abort the solve as STAGNATED while every
    // underlying metric was still decaying geometrically (channel n64,
    // diagnostic run 38690061509).
    ConvergenceMonitorControls transient_rise = controls;
    transient_rise.qoi_gates.clear();
    transient_rise.stagnation_window = 2;
    transient_rise.stagnation_relative_improvement = 1e-3;
    transient_rise.divergence_window = 3;
    transient_rise.divergence_growth_factor = 10.0;
    ConvergenceMonitor rising(transient_rise);
    rising.update({1, 1.0, 1.0, {}});
    rising.update({2, 1.005, 1.0, {}});
    auto rr = rising.update({3, 1.01, 1.0, {}});
    assert(rr.status == ConvergenceStatus::CONTINUE);

    // A flat envelope that includes micro-rises (|change| below the
    // threshold) still stagnates: 0.05% per iteration is not progress.
    ConvergenceMonitorControls micro_rise_controls = controls;
    micro_rise_controls.qoi_gates.clear();
    micro_rise_controls.stagnation_window = 2;
    micro_rise_controls.stagnation_relative_improvement = 1e-3;
    ConvergenceMonitor micro_rise(micro_rise_controls);
    micro_rise.update({1, 1.0, 1.0, {}});
    micro_rise.update({2, 1.0005, 1.0, {}});
    auto mr = micro_rise.update({3, 1.001, 1.0, {}});
    assert(mr.status == ConvergenceStatus::STAGNATED);

    ConvergenceMonitorControls divergence = controls;
    divergence.qoi_gates.clear();
    divergence.divergence_window = 2;
    divergence.divergence_growth_factor = 2.0;
    ConvergenceMonitor diverging(divergence);
    diverging.update({1, 1.0, 1.0, {}});
    diverging.update({2, 3.0, 1.0, {}});
    auto dr = diverging.update({3, 10.0, 1.0, {}});
    assert(dr.status == ConvergenceStatus::DIVERGED);

    bool threw = false;
    try {
        monitor.update({2, 1.0, 1.0, {}});
        monitor.update({1, 1.0, 1.0, {}});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    return 0;
}
