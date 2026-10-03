#include "cfdx/cfdx/core/utils/convergence_monitor.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

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
