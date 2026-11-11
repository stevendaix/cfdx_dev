#include "cfdx/physics/sst_outer_ptc_v31.h"
#include "cfdx/physics/sst_coupled_pseudo_time_v30.h"
#include "common/test_harness.h"
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("shared_cfl_is_committed_only_on_full_global_physical_contraction", []{
        SSTV31OuterPTC c;
        c.committed({1.,.1,.7,.5,.01},true);
        EXPECT_TRUE(c.committed_outer==1 && c.baseline_valid);
        EXPECT_NEAR(c.cfl,.5,0.);
        const auto snapshot=c;
        auto trial=c;
        EXPECT_THROW(trial.committed({1.1,.1,.6,.4,.01},true),std::invalid_argument);
        EXPECT_NEAR(c.cfl,snapshot.cfl,0.);
        trial=c; trial.committed({.79,.07,.50,.4,.008},true);
        EXPECT_TRUE(trial.committed_outer==2);
        EXPECT_NEAR(trial.cfl,.6,1e-15);
        EXPECT_NEAR(c.cfl,.5,0.); // no mutation before transaction commit
    });
    run_case("component_growth_shrinks_despite_global_max_improvement", []{
        SSTV31OuterPTC c;
        c.committed({1.,.1,.6,.5,.01},true);
        c.committed({.9,.12,.5,.4,.01},true);
        EXPECT_NEAR(c.cfl,.35,1e-15);
    });
    run_case("changed_discrete_operator_resets_comparison_without_acceptance_claim", []{
        SSTV31OuterPTC c;c.committed({1.,.1,.6,.5,.01},true);
        c.committed({2.,.2,1.2,1.,.02},false);
        EXPECT_NEAR(c.cfl,.5,0.);
        EXPECT_TRUE(c.committed_outer==2);
    });
    run_case("frozen_SST_CFL_cannot_be_promoted_by_inner_steps", []{
        SSTV30ContinuationState inner;
        inner.frozen_by_outer_v31=true;
        inner.accept_inner(1.,.1);
        inner.reject_inner();
        EXPECT_NEAR(inner.cfl,.5,0.);
        EXPECT_TRUE(inner.accepted_inner==1 && inner.rejected_inner==1);
        inner.frozen_by_outer_v31=false;
        inner.accept_inner(1.,.5);
        EXPECT_NEAR(inner.cfl,.6,1e-15);
    });
    run_case("globally_converged_plateau_is_admissible_without_CFL_increase", []{
        SSTV31OuterPTC c;
        c.committed({1e-11,1e-12,1e-11,1e-11,1e-12},true);
        const auto original_cfl=c.cfl;
        c.committed({1e-11,1e-12,1e-11,1e-11,1e-12},true,true);
        EXPECT_NEAR(c.cfl,original_cfl,0.);
        EXPECT_TRUE(c.committed_outer==2);
        EXPECT_THROW(c.committed({1e-11,1e-12,1e-11,1e-11,1e-12},true,false),std::invalid_argument);
    });
    run_case("invalid_physical_gate_and_bad_CFL_fail_closed", []{
        SSTV31OuterPTC c;
        EXPECT_THROW(c.committed({1.,-1.,1.,1.,1.},true),std::invalid_argument);
        EXPECT_TRUE(c.committed_outer==0);
        c.cfl=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(c.committed({1.,1.,1.,1.,1.},true),std::invalid_argument);
    });
    return run_all();
}
