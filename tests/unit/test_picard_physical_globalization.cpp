#include "cfdx/physics/picard_physical_globalization.h"

// The regression assertions must remain active in Release (-DNDEBUG).
#define REQUIRE(expr) do { \
    if (!(expr)) throw std::runtime_error("REQUIRE failed: " #expr); \
} while (false)
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

int main() {
    PicardPhysicalGlobalizationControls c;
    c.enabled = true;
    const PicardPhysicalSample old{
        {1e-3, 2e-4, 0.0}, {1.0, 0.01, 1e-300}, 1e-5, 1.0};
    // Both channels improve and the original, unrelaxed residual is checked.
    const PicardPhysicalSample good{
        {9e-4, 1.5e-4, 0.0}, {1000.0, 1e3, 1e-300}, 8e-6, 5.0};
    const auto old_merit = picard_physical_merit(old, old, 1e-8, 1e-8);
    const auto good_merit = picard_physical_merit(good, old, 1e-8, 1e-8);
    REQUIRE(std::abs(old_merit.momentum - 0.02) < 1e-14);
    REQUIRE(std::abs(good_merit.momentum - 0.015) < 1e-14);
    // Candidate has vastly larger RHS and face throughput; these must NOT
    // make the physical residual appear artificially smaller.
    REQUIRE(std::abs(good_merit.continuity - 8e-6) < 1e-16);
    REQUIRE(picard_physical_accept(old_merit, good_merit, 0.5, c));
    const PicardPhysicalSample bad{{2e-3, 2.1e-4, 0}, {1e8, 1e8, 1e-300},
                                   1.1e-5, 1e8};
    const auto bad_merit = picard_physical_merit(bad, old, 1e-8, 1e-8);
    REQUIRE(!picard_physical_accept(old_merit, bad_merit, 1.0, c));
    PicardPhysicalMerit stagnant{0.1, 0.1, 1e6};
    REQUIRE(!picard_physical_accept(stagnant, stagnant, 0.5, c));
    c.sufficient_decrease = 0;
    REQUIRE(picard_physical_accept(stagnant, stagnant, 0.5, c));
    c.sufficient_decrease = 1e-4;
    c.nonmonotone_window = 6;
    c.maximum_window_growth = 1.05;
    PicardPhysicalMerit worst_window{0.06, 0.02, 6e6};
    PicardPhysicalMerit trial_slight_rise{0.061, 0.0201, 6.10e6};
    REQUIRE(picard_physical_accept(worst_window, trial_slight_rise, 1.0, c));
    PicardPhysicalMerit trial_large_rise{0.10, 0.03, 6.70e6};
    REQUIRE(!picard_physical_accept(worst_window, trial_large_rise, 1.0, c));
    c.nonmonotone_window=1; c.maximum_window_growth=1.0;
    PicardPhysicalMerit near_gate{1e-8,1e-8,1.5};
    REQUIRE(picard_physical_accept(near_gate,near_gate,0.3,c));
    bool invalid = false;
    try { c.channel_growth_cap = 0.5; validate_picard_physical_globalization(c); }
    catch (const std::invalid_argument&) { invalid = true; }
    REQUIRE(invalid);
    std::cout << "PASS physical Picard globalization: fixed-normalization / "
                 "monotone / red-green / guard / unchanged gates\n";
}

