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
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 0.8, 0.2, 0.9, c), 0.52, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 1.2, 0.2, 0.9, c), 0.48, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 1.02, 0.2, 0.9, c), 0.5, 1e-12);

        // Slow but genuine progress creeps up rather than holding, so a factor
        // that an earlier transient drove down can recover. The creep band is
        // just inside the strong-improvement band, i.e. 0.95 < ratio < 1.
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 0.97, 0.2, 0.9, c), 0.51, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 0.94, 0.2, 0.9, c), 0.52, 1e-12);

        // Regression for the N7 Re=400 limit-cycle mechanism: a large
        // degradation must not collapse 0.726 directly to the floor. The
        // controller should recover gradually after subsequent improvements.
        double alpha = 0.726;
        alpha = adapt_relaxation_factor(alpha, 1.0, 2.33, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.706, 1e-12);
        alpha = adapt_relaxation_factor(alpha, 2.33, 1.95, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.726, 1e-12);
        alpha = adapt_relaxation_factor(alpha, 1.95, 0.9, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.746, 1e-12);

        // Steps are additive and clamped, so no single step can reach a bound
        // from a useful interior value regardless of how extreme the ratio is.
        AdaptiveRelaxationControls narrow = c;
        narrow.min_alpha_u = 0.6;
        narrow.max_alpha_u = 0.8;
        narrow.min_alpha_p = 0.6;
        narrow.max_alpha_p = 0.8;
        EXPECT_NEAR(adapt_relaxation_factor(0.61, 1e-9, 1e9, 0.6, 0.8, narrow), 0.6, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.79, 1e9, 1e-9, 0.6, 0.8, narrow), 0.8, 1e-12);
        // Repeated extreme degradation walks down one step at a time, never jumping.
        double walked = 0.79;
        for (int i = 0; i < 4; ++i)
            walked = adapt_relaxation_factor(walked, 1.0, 100.0, 0.6, 0.8, narrow);
        EXPECT_NEAR(walked, 0.71, 1e-12);
    });

    run_case("adaptive_reference", [] {
        AdaptiveRelaxationControls c;
        c.enabled = true;

        // The first measurement seeds the reference instead of dividing by it.
        EXPECT_NEAR(adaptive_reference_update(0.0, 5.0, c), 5.0, 1e-12);

        // A single step moves the reference by 1/window of the difference, so
        // one noisy measurement cannot swing it.
        const double one_step = adaptive_reference_update(5.0, 1.0, c);
        EXPECT_NEAR(one_step, 4.6, 1e-12);
        AdaptiveRelaxationControls fast = c;
        fast.reference_window = 2;
        EXPECT_NEAR(adaptive_reference_update(5.0, 1.0, fast), 3.0, 1e-12);

        // The reference is a convex combination of the metrics it has seen, so
        // it can never leave the range they occupy. The trend is what moves it.
        double reference = 0.0;
        for (int i = 0; i < 200; ++i) {
            const double metric = (i % 2 == 0) ? 1.0e-3 : 4.0e-3;
            reference = adaptive_reference_update(reference, metric, c);
            EXPECT_TRUE(reference >= 1.0e-3 - 1e-15);
            EXPECT_TRUE(reference <= 4.0e-3 + 1e-15);
        }
        EXPECT_TRUE(reference > 1.0e-3);
        EXPECT_TRUE(reference < 4.0e-3);

        // A sustained decay is read as improvement even though the reference
        // lags above the current metric, so the controller creeps up. Against a
        // geometric decay of factor f the steady-state ratio is
        // (1 - 0.9/f)/0.1, so the decay rate decides which branch is taken: a
        // 2%/iteration decay reads as strong improvement and a 0.33% decay as a
        // creep. A smoothed reference trades noise for lag, and the lag is
        // bounded by the window rather than by a step size.
        const auto settle = [](double factor, AdaptiveRelaxationControls controls) {
            double reference = 0.0;
            double current = 1.0;
            for (int i = 0; i < 4000; ++i) {
                current *= factor;
                reference = adaptive_reference_update(reference, current, controls);
            }
            return std::pair<double, double>(reference, current);
        };
        const auto fast_decay = settle(0.98, c);
        EXPECT_TRUE(fast_decay.first > fast_decay.second);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, fast_decay.first, fast_decay.second, 0.2, 0.9, c), 0.52, 1e-12);
        const auto slow_decay = settle(0.9967, c);
        EXPECT_TRUE(slow_decay.first > slow_decay.second);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, slow_decay.first, slow_decay.second, 0.2, 0.9, c), 0.51, 1e-12);

        // An isolated spike against a settled reference reads as degradation,
        // and costs exactly one bounded additive step.
        const double settled = adaptive_reference_update(0.0, 1.0, c);
        const double after_spike = adaptive_reference_update(settled, 100.0, c);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, after_spike, 100.0, 0.2, 0.9, c), 0.48, 1e-12);
    });

    run_case("adaptive_relaxation_validation", [] {
        AdaptiveRelaxationControls c;
        c.enabled = true;
        c.reference_window = 1;
        EXPECT_THROW(validate_adaptive_relaxation_controls(c), std::invalid_argument);
        c.reference_window = 2;
        validate_adaptive_relaxation_controls(c);
        c.reference_window = 0;
        EXPECT_THROW(validate_adaptive_relaxation_controls(c), std::invalid_argument);
        EXPECT_THROW(adaptive_reference_update(1.0, -1.0, c), std::invalid_argument);
        EXPECT_THROW(adapt_relaxation_factor(0.5, -1.0, 1.0, 0.2, 0.9, c), std::invalid_argument);
    });

    return run_all();
}
