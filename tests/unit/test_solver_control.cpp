#include "cfdx/physics/solver_control.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "common/test_harness.h"

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("convergence_criteria_validation", [] {
        ConvergenceCriteria c;
        validate_convergence_criteria(c);
        EXPECT_TRUE(residual_converged(1.0, 1e-7, c));
        EXPECT_FALSE(residual_converged(1.0, 1e-2, c));
    });

    run_case("conservation_gate", [] {
        ConvergenceCriteria c;
        ConvergenceCriteria c_energy = c;
        c_energy.require_energy = true;
        IterationMetrics m;
        m.continuity_imbalance = 1e-12;
        m.energy_imbalance = 1e-12;
        m.energy_imbalance = 1e-12;
        EXPECT_TRUE(conservation_converged(m, c_energy));
        m.energy_imbalance = 1e-3;
        EXPECT_FALSE(conservation_converged(m, c_energy));
    });

    run_case("adaptive_relaxation", [] {
        AdaptiveRelaxationControls c;
        c.enabled = true;
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 0.8, 0.2, 0.9, c), 0.505, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 1.2, 0.2, 0.9, c), 0.48, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 1.02, 0.2, 0.9, c), 0.5, 1e-12);

        // Regression for the N7 Re=400 limit-cycle mechanism: degradation must
        // not collapse 0.726 directly to the floor; recovery remains gradual.
        double alpha = 0.726;
        alpha = adapt_relaxation_factor(alpha, 1.0, 2.33, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.706, 1e-12);
        alpha = adapt_relaxation_factor(alpha, 2.33, 1.95, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.711, 1e-12);
        alpha = adapt_relaxation_factor(alpha, 1.95, 0.9, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.716, 1e-12);
    });

    run_case("adaptive_relaxation_severe_overshoot_guard", [] {
        AdaptiveRelaxationControls c;
        c.enabled = true;

        // A 2x residual jump is severe only when the corresponding
        // relaxation factor increased between accepted states.
        EXPECT_TRUE(adaptive_relaxation_severe_degradation(1.0, 2.0, 0.70, 0.705, c));
        EXPECT_FALSE(adaptive_relaxation_severe_degradation(1.0, 1.99, 0.70, 0.705, c));
        EXPECT_FALSE(adaptive_relaxation_severe_degradation(1.0, 2.0, 0.70, 0.70, c));
        EXPECT_FALSE(adaptive_relaxation_severe_degradation(1.0, 2.0, 0.705, 0.70, c));

        // The guard is channel-local: a momentum alpha increase does not imply
        // a pressure overshoot and vice versa.
        EXPECT_TRUE(adaptive_relaxation_severe_degradation(1.0, 2.1, 0.30, 0.31, c));
        EXPECT_FALSE(adaptive_relaxation_severe_degradation(1.0, 2.1, 0.31, 0.30, c));

        double alpha = 0.78;
        alpha = std::max(c.min_alpha_u, alpha - c.decrease_step);
        EXPECT_NEAR(alpha, 0.76, 1e-12);

        c.recovery_cooldown_windows = 2;
        EXPECT_TRUE(c.recovery_cooldown_windows == 2);
    });

    run_case("adaptive_relaxation_windowed_hysteresis", [] {
        AdaptiveRelaxationControls c;
        c.enabled = true;
        c.required_trend_windows = 2;

        std::size_t improve = 0;
        std::size_t degrade = 0;
        double alpha = 0.70;

        // One improving window is insufficient to change alpha.
        alpha = adapt_relaxation_factor_windowed(
            alpha, 1.0, 0.90, 0.6, 0.8, c, improve, degrade);
        EXPECT_NEAR(alpha, 0.70, 1e-12);
        EXPECT_TRUE(improve == 1);
        EXPECT_TRUE(degrade == 0);

        // A second sustained improvement changes alpha exactly once.
        alpha = adapt_relaxation_factor_windowed(
            alpha, 1.0, 0.90, 0.6, 0.8, c, improve, degrade);
        EXPECT_NEAR(alpha, 0.705, 1e-12);
        EXPECT_TRUE(improve == 0);

        // A neutral window clears the trend and prevents stale decisions.
        alpha = adapt_relaxation_factor_windowed(
            alpha, 1.0, 1.02, 0.6, 0.8, c, improve, degrade);
        EXPECT_NEAR(alpha, 0.705, 1e-12);
        EXPECT_TRUE(improve == 0);
        EXPECT_TRUE(degrade == 0);

        // Two degrading windows reduce the same factor, without touching the
        // independently tracked pressure controller state.
        alpha = adapt_relaxation_factor_windowed(
            alpha, 1.0, 1.20, 0.6, 0.8, c, improve, degrade);
        EXPECT_NEAR(alpha, 0.705, 1e-12);
        EXPECT_TRUE(degrade == 1);
        alpha = adapt_relaxation_factor_windowed(
            alpha, 1.0, 1.20, 0.6, 0.8, c, improve, degrade);
        EXPECT_NEAR(alpha, 0.685, 1e-12);
        EXPECT_TRUE(degrade == 0);
    });

    run_case("adaptive_relaxation_independent_channels", [] {
        AdaptiveRelaxationControls c;
        c.enabled = true;
        c.required_trend_windows = 1;

        std::size_t u_improve = 0, u_degrade = 0;
        std::size_t p_improve = 0, p_degrade = 0;
        double alpha_u = 0.70;
        double alpha_p = 0.30;

        // Momentum improves while continuity degrades. The two controls must
        // move in opposite directions rather than being driven by one scalar
        // aggregate metric.
        alpha_u = adapt_relaxation_factor_windowed(
            alpha_u, 1.0, 0.90, 0.6, 0.8, c, u_improve, u_degrade);
        alpha_p = adapt_relaxation_factor_windowed(
            alpha_p, 1.0, 1.20, 0.25, 0.35, c, p_improve, p_degrade);

        EXPECT_NEAR(alpha_u, 0.705, 1e-12);
        EXPECT_NEAR(alpha_p, 0.28, 1e-12);
    });

    return run_all();
}
