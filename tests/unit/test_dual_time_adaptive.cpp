#include "cfdx/physics/dual_time_adaptive.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    using namespace cfdx::physics;

    const std::vector<double> high{1.0, 2.0};
    const std::vector<double> low{0.99, 2.0};
    const auto error = estimate_dual_time_temporal_error(high, low, 1.0e-3, 1.0e-2);
    require(error.absolute_norm > 0.0, "temporal error estimate must be positive");
    require(error.normalized_error > 0.0, "normalized temporal error must be positive");

    const auto accepted = propose_dual_time_physical_step(
        0.1, 0.1, 2, 1.0e-3, 1.0);
    require(accepted.accept, "small temporal error should accept the physical step");
    require(accepted.dt > 0.1, "accepted accurate step should grow");

    const auto rejected = propose_dual_time_physical_step(
        0.1, 10.0, 2, 1.0e-3, 1.0);
    require(!rejected.accept, "large temporal error must reject the physical step");
    require(rejected.dt < 0.1, "rejected step should shrink");

    const auto pseudo = adapt_dual_time_pseudo_dt(
        1.0, 1.0, 0.1, 1.0e-3, 10.0);
    require(pseudo.contracted, "residual contraction must be detected");
    require(pseudo.contraction < 1.0, "contraction ratio must be below one");
    require(pseudo.pseudo_dt > 1.0, "strong contraction should increase pseudo-time step");

    const DualTimeStepAcceptance gate{true, true, false};
    require(!gate.accepted(), "inadmissible state must reject the physical step");

    std::cout << "dual-time adaptive controller tests: PASS\n";
    return 0;
}
