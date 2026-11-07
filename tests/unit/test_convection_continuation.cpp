#include "cfdx/physics/convection_continuation.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("disabled_continuation_has_full_scheme_from_iteration_one", [] {
        ConvectionContinuationControls c;
        EXPECT_NEAR(convection_high_order_fraction(c,1),1.0,0.0);
        EXPECT_TRUE(convection_full_order_eligible(c,1));
    });
    run_case("ramp_monotone_and_final_gate_after_full_order_hold", [] {
        ConvectionContinuationControls c;
        c.enabled=true; c.upwind_warmup_iterations=2;
        c.high_order_ramp_iterations=4; c.full_order_hold_iterations=3;
        validate_convection_continuation(c,8);
        EXPECT_NEAR(convection_high_order_fraction(c,1),0.0,0.0);
        EXPECT_NEAR(convection_high_order_fraction(c,2),0.0,0.0);
        EXPECT_NEAR(convection_high_order_fraction(c,3),0.25,0.0);
        EXPECT_NEAR(convection_high_order_fraction(c,4),0.5,0.0);
        EXPECT_NEAR(convection_high_order_fraction(c,5),0.75,0.0);
        EXPECT_NEAR(convection_high_order_fraction(c,6),1.0,0.0);
        EXPECT_TRUE(!convection_full_order_eligible(c,6));
        EXPECT_TRUE(!convection_full_order_eligible(c,7));
        EXPECT_TRUE(convection_full_order_eligible(c,8));
    });
    run_case("nonlinear_growth_guard_never_compares_different_operators", [] {
        ConvectionContinuationControls c;
        for (std::size_t it = 1; it <= 10; ++it)
            EXPECT_TRUE(convection_residuals_share_operator(c,it));
        c.enabled=true;
        c.upwind_warmup_iterations=2;
        c.high_order_ramp_iterations=4;
        c.full_order_hold_iterations=3;
        EXPECT_TRUE(convection_residuals_share_operator(c,1));
        EXPECT_TRUE(convection_residuals_share_operator(c,2));
        EXPECT_TRUE(!convection_residuals_share_operator(c,3));
        EXPECT_TRUE(!convection_residuals_share_operator(c,4));
        EXPECT_TRUE(!convection_residuals_share_operator(c,5));
        EXPECT_TRUE(!convection_residuals_share_operator(c,6));
        EXPECT_TRUE(convection_residuals_share_operator(c,7));
        EXPECT_TRUE(convection_residuals_share_operator(c,8));
        EXPECT_THROW(convection_residuals_share_operator(c,0), std::invalid_argument);
    });
    run_case("ramp_refuses_budget_that_can_only_converge_at_partial_order", [] {
        ConvectionContinuationControls c;
        c.enabled=true; c.upwind_warmup_iterations=2;
        EXPECT_TRUE(!convection_full_order_eligible(c,7));
        EXPECT_TRUE(convection_full_order_eligible(c,8));
    });
    run_case("ramp_refuses_budget_that_can_only_converge_at_partial_order", [] {
        ConvectionContinuationControls c;
        c.enabled=true; c.upwind_warmup_iterations=2;
        c.high_order_ramp_iterations=4; c.full_order_hold_iterations=3;
        EXPECT_THROW(validate_convection_continuation(c,7),std::invalid_argument);
        c.high_order_ramp_iterations=0;
        EXPECT_THROW(validate_convection_continuation(c,20),std::invalid_argument);
    });
    run_case("full_order_monitor_uses_normalized_physical_gate_not_pre_ramp_ratio", [] {
        ConvectionContinuationControls c;
        EXPECT_NEAR(convection_monitor_residual_threshold(c,1e-10,1e-8),1e-10,0.0);
        c.enabled=true;
        EXPECT_NEAR(convection_monitor_residual_threshold(c,1e-10,1e-8),1e-8,0.0);
        // A final, independently RHS-normalized nonlinear metric below
        // 1e-8 is acceptable; a pre-ramp reference ratio is not an
        // additional equation nor a valid final-operator comparison.
        const double final_physical_metric = 5.907e-9;
        const double pre_ramp_ratio = 1.80807e-5;
        EXPECT_TRUE(final_physical_metric <
            convection_monitor_residual_threshold(c,1e-10,1e-8));
        EXPECT_TRUE(pre_ramp_ratio >
            convection_monitor_residual_threshold(c,1e-10,1e-8));
        EXPECT_THROW(convection_monitor_residual_threshold(c,0.0,1e-8),
                     std::invalid_argument);
    });
    run_case("ramp_refuses_overflow", [] {
        ConvectionContinuationControls c;
        c.enabled=true;
        c.upwind_warmup_iterations=std::numeric_limits<std::size_t>::max();
        EXPECT_THROW(validate_convection_continuation(c,
            std::numeric_limits<std::size_t>::max()),std::invalid_argument);
    });
    return run_all();
}
