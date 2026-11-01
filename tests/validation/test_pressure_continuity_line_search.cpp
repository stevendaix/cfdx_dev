#include "cfdx/physics/pressure_continuity_line_search.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::physics;

static void require(bool value, const char* msg) {
    if (!value) throw std::runtime_error(msg);
}
int main() {
    try {
        PressureContinuityLineSearchControls c;
        c.enabled = true;
        c.max_trials = 5;
        auto result = select_pressure_continuity_step(1.0, c,
            [](double alpha) { return (alpha-0.25)*(alpha-0.25); });
        require(std::abs(result.alpha-0.25)<1e-14, "missed minimum at alpha=0.25");
        require(result.accepted_norm < result.baseline_norm, "no improvement");
        // A case whose correction always deteriorates must NOT pass.
        bool refused = false;
        try {
            (void)select_pressure_continuity_step(1.0, c,
                [](double alpha) { return 1.0 + alpha; });
        } catch (const std::runtime_error&) { refused = true; }
        require(refused, "accepted nondecreasing continuity");
        refused = false;
        try {
            (void)select_pressure_continuity_step(1.0, c,
                [](double) { return std::numeric_limits<double>::quiet_NaN(); });
        } catch (const std::runtime_error&) { refused = true; }
        require(refused, "accepted NaN");
        refused = false;
        c.reduction = 1.0;
        try { validate_pressure_continuity_line_search(c); }
        catch (const std::invalid_argument&) { refused = true; }
        require(refused, "accepted invalid reduction");
        std::cout << "PRESSURE_CONTINUITY_LINE_SEARCH: PASS\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "PRESSURE_CONTINUITY_LINE_SEARCH: FAIL " << ex.what() << '\n';
        return 1;
    }
}
