#include "cfdx/physics/continuation.h"
#include "common/test_harness.h"

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("continuation_controls_validation", [] {
        ContinuationControls c;
        validate_continuation_controls(c);
        EXPECT_NEAR(continuation_next_target(0.0, 0.25), 0.25, 1e-12);
        EXPECT_NEAR(continuation_next_target(0.9, 0.25), 1.0, 1e-12);
    });

    run_case("continuation_step_policy", [] {
        ContinuationControls c;
        c.initial_step = 0.25;
        c.maximum_step = 0.5;
        c.minimum_step = 0.01;
        c.step_growth = 1.5;
        c.step_reduction = 0.5;
        EXPECT_NEAR(continuation_step_after_success(0.25, c), 0.375, 1e-12);
        EXPECT_NEAR(continuation_step_after_success(0.5, c), 0.5, 1e-12);
        EXPECT_NEAR(continuation_step_after_failure(0.25, c), 0.125, 1e-12);
    });

    run_case("continuation_invalid_controls", [] {
        ContinuationControls c;
        c.step_reduction = 1.0;
        bool rejected = false;
        try {
            validate_continuation_controls(c);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });

    return run_all();
}
