// Core temporal integration tests.
#include "cfdx/core/numerics/temporal.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

static void decay_rhs(const Field<double, Location::CELL>& phi,
                      Field<double, Location::CELL>& rhs) {
    for (std::size_t i = 0; i < phi.size(); ++i) {
        rhs.component_data(0)[i] = -2.0 * phi.component_data(0)[i];
    }
}

int main() {
    run_case("core_implicit_euler", []() {
        Field<double, Location::CELL> phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        auto out = advance_time(TimeScheme::EULER_IMPLICIT, phi, 0.1, decay_rhs);
        EXPECT_NEAR(out(0), 1.0 / 1.2, 1e-10);
    });

    run_case("core_crank_nicolson", []() {
        Field<double, Location::CELL> phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        auto out = advance_time(TimeScheme::CRANK_NICOLSON, phi, 0.1, decay_rhs);
        EXPECT_NEAR(out(0), 0.9 / 1.1, 1e-10);
    });

    run_case("core_invalid_dt", []() {
        Field<double, Location::CELL> phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        bool threw = false;
        try {
            (void)advance_time(TimeScheme::EULER_EXPLICIT, phi, 0.0, decay_rhs);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        EXPECT_TRUE(threw);
    });

    return run_all();
}
