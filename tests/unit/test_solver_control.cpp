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
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 1.2, 0.2, 0.9, c), 0.45, 1e-12);
        EXPECT_NEAR(adapt_relaxation_factor(0.5, 1.0, 1.02, 0.2, 0.9, c), 0.5, 1e-12);

        // Regression for the N7 Re=400 limit-cycle mechanism: degradation must
        // not collapse 0.726 directly to the floor; recovery remains gradual.
        double alpha = 0.726;
        alpha = adapt_relaxation_factor(alpha, 1.0, 2.33, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.676, 1e-12);
        alpha = adapt_relaxation_factor(alpha, 2.33, 1.95, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.696, 1e-12);
        alpha = adapt_relaxation_factor(alpha, 1.95, 0.9, 0.6, 0.8, c);
        EXPECT_NEAR(alpha, 0.716, 1e-12);
    });

    return run_all();
}
