#include "cfdx/physics/momentum_residual_certification.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
int main() {
    using cfdx::physics::certify_momentum_residual_components;
    // Old single maximum RHS would accept y even when its own equation has
    // a two-orders-of-magnitude larger relative defect.
    const auto missed = certify_momentum_residual_components(
        {1.e-5, 5.e-4, 1.e-18}, {100., 1., 1.e-19}, 1.e-5);
    if (missed.converged || !(missed.global_relative < 1.e-5) ||
        std::abs(missed.component_relative[1] - 5.e-4) > 1.e-16 ||
        !missed.active[0] || !missed.active[1] || missed.active[2]) return 1;
    const auto passed = certify_momentum_residual_components(
        {1.e-6, 1.e-7, 1.e-16}, {100., 1., 1.e-19}, 1.e-6);
    if (!passed.converged) return 2;
    const auto dormant_failure = certify_momentum_residual_components(
        {1.e-6, 1.e-7, 1.e-1}, {100., 1., 1.e-19}, 1.e-6);
    if (dormant_failure.converged) return 3; // global gate still checks Uz
    bool rejected = false;
    try { (void)certify_momentum_residual_components({0,0,0},{1,1,1},std::numeric_limits<double>::quiet_NaN()); }
    catch(const std::invalid_argument&) { rejected = true; }
    if (!rejected) return 4;
    std::cout << "PASS V11 componentwise and global momentum gates\n";
}
