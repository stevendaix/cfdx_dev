#include "cfdx/physics/pressure_affine_flux_trials.h"
#include "cfdx/physics/pressure_continuity_line_search.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

static void require(bool ok, const char* why) {
    if (!ok) throw std::runtime_error(why);
}

int main() {
    using namespace cfdx::physics;
    // Two internal/boundary cells with a fixed inflow, internal and outflow
    // flux.  The boundary faces MUST participate in the local balance.
    const std::vector<std::size_t> owner{0,0,1,1};
    const std::vector<std::int64_t> neighbour{-1,1,-1,-1};
    const std::vector<double> phi0{-2.0,1.0,3.0,0.0};
    const std::vector<double> phi1{-2.0,2.0,2.0,0.0};
    PressureAffineFluxTrials a(phi0, phi1, owner, neighbour, 2, 0.8);
    require(std::abs(a.continuity_norm(0.0)-3.0)<1e-14, "affine flux test failed");
    require(std::abs(a.continuity_norm(0.8))<1e-14, "affine flux test failed");
    require(std::abs(a.continuity_norm(0.4)-1.5)<1e-14, "affine flux test failed");
    require(a.agrees_with_flux(0.4, {-2.0,1.5,2.5,0.0}), "affine flux test failed");
    require(!a.agrees_with_flux(0.4, {-2.0,1.5,4.5,0.0}), "affine flux test failed");
    require(!a.agrees_with_flux(0.4, {-2.0,1.5,2.5,std::numeric_limits<double>::quiet_NaN()}), "affine flux test failed");
    PressureContinuityLineSearchControls controls;
    controls.enabled = true;
    const auto chosen = select_pressure_continuity_step(0.8, controls,
        [&](double alpha) { return a.continuity_norm(alpha); });
    require(std::abs(chosen.alpha-0.8)<1e-14, "affine flux test failed");
    require(chosen.accepted_norm < chosen.baseline_norm, "affine flux test failed");
    bool rejected = false;
    try { (void)a.continuity_norm(1.2); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "affine flux test failed");
    rejected = false;
    try { PressureAffineFluxTrials bad(phi0, phi1, owner, {0,1,900,-1}, 2, 0.8); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "affine flux test failed");
    std::cout << "PRESSURE_AFFINE_FLUX_TRIALS: PASS\n";
}
