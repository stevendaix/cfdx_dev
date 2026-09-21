#include "cfdx/physics/adaptive_cfl.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("adaptive_cfl_hits_target_without_limits", [] {
        AdaptiveCflControls c;
        c.target_cfl = 1.0;
        c.growth_limit = 10.0;
        c.shrink_limit = 0.1;
        EXPECT_NEAR(adaptive_time_step(1.0, 4.0, c), 0.5, 1e-14);
    });

    run_case("adaptive_cfl_respects_bounds", [] {
        AdaptiveCflControls c;
        c.min_dt = 0.25;
        c.max_dt = 0.75;
        c.growth_limit = 10.0;
        c.shrink_limit = 0.1;
        EXPECT_NEAR(adaptive_time_step(1.0, 100.0, c), 0.25, 1e-14);
        EXPECT_NEAR(adaptive_time_step(1.0, 0.001, c), 0.75, 1e-14);
    });

    run_case("adaptive_cfl_rejects_invalid_controls", [] {
        AdaptiveCflControls c;
        c.target_cfl = 0.0;
        EXPECT_THROW(adaptive_time_step(1.0, 1.0, c), std::invalid_argument);

        c = {};
        c.min_dt = 2.0;
        c.max_dt = 1.0;
        EXPECT_THROW(adaptive_time_step(1.0, 1.0, c), std::invalid_argument);

        c = {};
        c.growth_limit = 0.9;
        EXPECT_THROW(adaptive_time_step(1.0, 1.0, c), std::invalid_argument);

        c = {};
        c.shrink_limit = 1.1;
        EXPECT_THROW(adaptive_time_step(1.0, 1.0, c), std::invalid_argument);
    });

    run_case("adaptive_cfl_rejects_nonfinite_state", [] {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(adaptive_time_step(nan, 1.0), std::invalid_argument);
        EXPECT_THROW(adaptive_time_step(1.0, nan), std::invalid_argument);
    });

    run_case("pseudo_transient_cfl_is_bounded", [] {
        EXPECT_NEAR(pseudo_transient_cfl(0), 0.5, 1e-14);
        EXPECT_TRUE(pseudo_transient_cfl(1000) <= 100.0);
    });

    run_case("pseudo_transient_cfl_rejects_invalid_controls", [] {
        EXPECT_THROW(pseudo_transient_cfl(0, 0.0), std::invalid_argument);
        EXPECT_THROW(pseudo_transient_cfl(0, 0.5, 1.0), std::invalid_argument);
        EXPECT_THROW(pseudo_transient_cfl(0, 0.5, 1.2, 0.4), std::invalid_argument);
    });

    return run_all();
}
