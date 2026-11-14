#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_adaptive_trust_v28.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;
namespace {
SSTV25PhysicalTrial physical(const std::vector<double>& x) {
    SSTV25PhysicalTrial t;
    t.fields=x;
    t.k_local=1./x[0];t.omega_local=1./x[1];
    t.merit=std::max(t.k_local,t.omega_local);
    t.smooth_l2=std::hypot(t.k_local,t.omega_local)/std::sqrt(2.);
    t.k_worst=0;t.omega_worst=0;
    return t;
}
SSTV28BlendSearch blend(const SSTV25PhysicalEvaluator& f=physical) {
    return sst_blend_adaptive_trust_v28(
        {1.,1.},{1.5,1.05},{1.05,1.5},1,1e-9,1e-9,1.,1.,0.5,f);
}
}
int main() {
    run_case("v28_log_mixture_beats_both_endpoints_on_physical_linf",[]{
        const auto r=blend();
        EXPECT_TRUE(r.accepted && r.fields.size()==2);
        EXPECT_TRUE(r.merit<1./1.05 && r.smooth_l2<1.);
        EXPECT_TRUE(r.weight_jacobi>0. && r.weight_jacobi<1.);
        EXPECT_TRUE(r.actual_over_predicted>1e-3);
        EXPECT_TRUE(r.actual_reduction>0. && r.predicted_reduction>0.);
        EXPECT_TRUE(r.used_log_radius>0. && r.used_log_radius<=0.5*(1.+1e-9));
        EXPECT_TRUE(r.evaluations>1 && r.admissible_trials>0);
    });
    run_case("v28_does_not_mutate_the_caller_even_after_failed_trials",[]{
        const std::vector<double> x{1.,1.},a{1.5,1.05},b{1.05,1.5};
        auto r=sst_blend_adaptive_trust_v28(x,a,b,1,1e-9,1e-9,
            1.,1.,0.5,[](const auto& z) {
                if (z[0]>1.001) throw std::runtime_error("wall update failed");
                return physical(z);
            });
        EXPECT_TRUE(x==std::vector<double>({1.,1.}));
        EXPECT_TRUE(!r.accepted && r.fields.empty());
    });
    run_case("v28_rejects_l2_improvement_without_strict_linf",[]{
        auto r=blend([](const auto& x) {
            auto t=physical(x);
            t.k_local=1.;t.omega_local=1./x[1];
            t.merit=1.;
            t.smooth_l2=std::hypot(t.k_local,t.omega_local)/std::sqrt(2.);
            return t;
        });
        EXPECT_TRUE(!r.accepted && r.fields.empty());
        EXPECT_TRUE(r.rejected_linf>0);
    });
    run_case("v28_wall_projection_must_obey_trust_region",[]{
        auto r=blend([](const auto& x) {
            auto t=physical(x);
            if (x[0]>1.00001) {
                t.fields[0]=x[0]*4.;
                t.k_local=0.25;t.omega_local=0.25;
                t.merit=0.25;t.smooth_l2=0.25;
            }
            return t;
        });
        EXPECT_TRUE(!r.accepted && r.rejected_projection>0);
    });
    run_case("v28_stale_baseline_or_bad_fv_certificate_fails_closed",[]{
        const auto r=sst_blend_adaptive_trust_v28({1.,1.},{1.5,1.05},
            {1.05,1.5},1,1e-9,1e-9,0.5,1.,0.5,physical);
        EXPECT_TRUE(!r.accepted && r.evaluations==1);
        const auto bad=blend([](const auto& x) {
            auto t=physical(x);
            t.merit=std::numeric_limits<double>::quiet_NaN();return t;
        });
        EXPECT_TRUE(!bad.accepted && bad.evaluations==1);
    });
    run_case("v28_radius_commits_only_explicitly_after_actual_selection",[]{
        SSTV28AdaptiveTrustState state;
        auto r=blend();
        EXPECT_TRUE(r.accepted);
        // Merely constructing a valid proposal never advances the radius.
        EXPECT_TRUE(state.radius==0.5 && state.committed==0);
        state.commit(0.9,0.5);
        EXPECT_TRUE(state.radius==1.0 && state.committed==1);
        state.commit(0.1,0.2);
        EXPECT_TRUE(state.radius==0.5 && state.committed==2);
        state.reject();
        EXPECT_TRUE(state.radius==0.25 && state.rejected==1);
        EXPECT_THROW(state.commit(std::numeric_limits<double>::quiet_NaN(),0.2),
                     std::invalid_argument);
        EXPECT_TRUE(state.committed==2 && state.radius==0.25);
    });
    run_case("v28_rejected_ratio_cannot_expand_state",[]{
        SSTV28AdaptiveTrustState state;
        EXPECT_THROW(state.commit(-0.5,0.1),std::invalid_argument);
        EXPECT_THROW(state.commit(1.0,1.5),std::invalid_argument);
        EXPECT_TRUE(state.radius==0.5 && state.committed==0);
        for(int i=0;i<10;++i) state.reject();
        EXPECT_TRUE(state.radius==SSTV28AdaptiveTrustState::minimum_radius);
    });
    run_case("v28_rejects_nonpositive_fields_and_out_of_range_radius",[]{
        EXPECT_THROW(sst_blend_adaptive_trust_v28({1.,1.},{0.,1.1},
            {1.1,1.2},1,1e-9,1e-9,1.,1.,0.5,physical),std::invalid_argument);
        EXPECT_THROW(sst_blend_adaptive_trust_v28({1.,1.},{1.1,1.1},
            {1.1,1.2},1,1e-9,1e-9,1.,1.,2.5,physical),std::invalid_argument);
    });
    run_case("v28_stall_rescue_is_bounded_and_has_no_default_overhead",[]{
        EXPECT_TRUE(!sst_v28_stagnation_probe(false,22,1));
        EXPECT_TRUE(!sst_v28_stagnation_probe(true,20,1));
        EXPECT_TRUE(!sst_v28_stagnation_probe(true,21,0));
        EXPECT_TRUE(sst_v28_stagnation_probe(true,22,1));
        EXPECT_TRUE(sst_v28_stagnation_probe(true,23,2));
        EXPECT_TRUE(!sst_v28_stagnation_probe(true,24,3));
    });
    run_case("v28_optin_requires_v27_and_is_disabled_by_default",[]{
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_adaptive_trust_v28_experimental);
        c.sst_adaptive_trust_v28_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_subspace_trust_v27_experimental=true;
        c.sst_smooth_globalization_v26_experimental=true;
        c.sst_log_globalization_v25_experimental=true;
        c.sst_spatial_block_preconditioner_v24_experimental=true;
        c.sst_jfnk_v22_experimental=true;
        validate_turbulence_controls(c);
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(0.5,0.50001,0.1));
    });
    return run_all();
}
