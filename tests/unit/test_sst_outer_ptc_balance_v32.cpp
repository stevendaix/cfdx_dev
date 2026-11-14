#include "cfdx/physics/sst_outer_ptc_balance_v32.h"
#include "cfdx/physics/sst_coupled_pseudo_time_v30.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("both_transport_blocks_are_bounded_relative_to_their_steady_diagonal", [] {
        const auto mass=sst_v30_coupled_pseudo_mass(
            {4.,5.,2.,7.},{5.,8.,2.,14.},{12.,1.,2.,9.},{0,0,1,0},.5,1./32.);
        EXPECT_NEAR(mass[0],4./32.,1e-14);
        EXPECT_NEAR(mass[1],5./32.,1e-14);
        EXPECT_NEAR(mass[2],0.,0.);
        EXPECT_NEAR(mass[3],7./32.,1e-14);
        const auto legacy=sst_v30_coupled_pseudo_mass(
            {4.,5.,2.,7.},{5.,8.,2.,14.},{12.,1.,2.,9.},{0,0,1,0},.5);
        EXPECT_NEAR(legacy[0],12.,0.);
        EXPECT_NEAR(legacy[1],10.,0.);
    });
    run_case("unpublished_retries_monotonically_reduce_both_budgets", [] {
        SSTV32CouplingBudget b;
        const auto original=b;
        EXPECT_NEAR(b.trial_fraction(0),1.e-6,0.);
        EXPECT_NEAR(b.trial_fraction(1),5.e-7,0.);
        EXPECT_NEAR(b.trial_fraction(2),2.5e-7,0.);
        EXPECT_NEAR(b.trial_fraction(100),SSTV32CouplingBudget::minimum_fraction,0.);
        EXPECT_NEAR(b.fraction,original.fraction,0.);
        EXPECT_TRUE(b.committed_outer==0);
    });
    run_case("only_a_certified_outer_contraction_can_adapt_the_budget", [] {
        SSTV32CouplingBudget b;
        EXPECT_THROW(b.committed(1.,1.1),std::invalid_argument);
        EXPECT_NEAR(b.fraction,1.e-6,0.);
        auto trial=b;
        trial.committed(1.,0.7);
        EXPECT_NEAR(trial.fraction,1.2e-6,1e-16);
        EXPECT_NEAR(b.fraction,1.e-6,0.);
        trial.committed(0.7,.69);
        EXPECT_NEAR(trial.fraction,.75*1.2e-6,1e-16);
        EXPECT_TRUE(trial.committed_outer==2);
        EXPECT_THROW(trial.committed(0.,0.),std::invalid_argument);
    });
    run_case("invalid_zero_or_negative_steady_diagonal_or_fraction_fails_closed", [] {
        const double nan=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(sst_v30_coupled_pseudo_mass({0.},{1.},{1.},{},.5,1./32.),std::invalid_argument);
        EXPECT_THROW(sst_v30_coupled_pseudo_mass({4.},{1.},{1.},{},.5,nan),std::invalid_argument);
        SSTV32CouplingBudget c;c.fraction=-.1;
        EXPECT_THROW(c.trial_fraction(0),std::invalid_argument);
    });
    run_case("local_SST_micro_CFL_can_evolve_but_a_rejected_outer_trial_discards_it", [] {
        SSTV30ContinuationState committed;
        const auto baseline=committed;
        auto attempt=committed;
        attempt.frozen_by_outer_v31=false;
        attempt.v32_max_diagonal_fraction=1.e-6;
        attempt.accept_inner(1.0,0.7);
        EXPECT_TRUE(attempt.cfl>baseline.cfl);
        EXPECT_NEAR(committed.cfl,baseline.cfl,0.);
        EXPECT_TRUE(committed.accepted_inner==0);
        // A physical rollback always restores the pretrial controller.
        attempt=baseline;
        EXPECT_NEAR(attempt.cfl,committed.cfl,0.);
        EXPECT_TRUE(attempt.accepted_inner==0);
    });
    run_case("V32_is_opt_in_and_requires_the_full_V31_physical_gate", [] {
        IncompressibleSolverControls c;
        EXPECT_TRUE(!c.outer_ptc_balance_v32_experimental);
        EXPECT_TRUE(c.outer_ptc_balance_momentum_v32 &&
                    !c.outer_ptc_balance_sst_v32);
        c.outer_ptc_balance_v32_experimental=true;
        EXPECT_THROW(validate_incompressible_controls(c,4),std::invalid_argument);
        c.outer_ptc_v31_experimental=true;
        c.outer_ptc_balance_momentum_v32=false;
        c.outer_ptc_balance_sst_v32=false;
        EXPECT_THROW(validate_incompressible_controls(c,4),std::invalid_argument);
    });
    return run_all();
}
