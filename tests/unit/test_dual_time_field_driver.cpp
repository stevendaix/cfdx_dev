#include "cfdx/physics/dual_time_field_driver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace cfdx::core;
    using namespace cfdx::physics;

    ScalarCellField state(2, "u", "1", 1);
    state(0) = 1.0; state(1) = 2.0;
    ScalarCellField high = state;
    ScalarCellField low = state;
    high(0) = 1.001;

    DualTimeFieldHistory history(state);
    history.accept(high, 0.1);
    if (!history.has_previous || std::abs(history.dt_previous - 0.1) > 1e-14)
        throw std::runtime_error("dual-time field history did not commit");

    const DualTimeStepAcceptance accepted{true, true, true};
    const auto result = accept_dual_time_field_step(
        high, low, high, 0.1, accepted, 1.0e-6, 1.0e-3);
    if (!result.acceptance.accepted() || result.state(0) != high(0))
        throw std::runtime_error("accepted dual-time field step is incorrect");

    const DualTimeStepAcceptance rejected{true, true, false};
    bool failed = false;
    try {
        (void)accept_dual_time_field_step(
            high, low, high, 0.1, rejected, 1.0e-6, 1.0e-3);
    } catch (const DualTimeConvergenceFailure&) {
        failed = true;
    }
    if (!failed)
        throw std::runtime_error("inadmissible field step was not rejected");

    std::cout << "dual-time field driver tests: PASS\n";
    return 0;
}
