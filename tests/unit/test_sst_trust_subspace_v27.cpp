#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_trust_subspace_v27.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
SSTV25PhysicalTrial scalar_trial(const std::vector<double>& x) {
    SSTV25PhysicalTrial t;
    t.fields=x;
    t.k_local=std::exp(-std::log(x[0]));
    t.omega_local=0.5*std::exp(-std::log(x[3]));
    t.merit=std::max(t.k_local,t.omega_local);
    t.smooth_l2=std::hypot(t.k_local,t.omega_local)/std::sqrt(2.);
    t.k_worst=0;t.omega_worst=1;
    return t;
}
SSTV27TrustSearch solve(const SSTV25PhysicalEvaluator& eval=scalar_trial) {
    return sst_cauchy_subspace_v27(
        {1.,1.,1.,1.},{0.,1.,0.,0.},{1.2,1.,1.,1.2},2,
        1e-8,1e-8,1.,std::sqrt(.625),0,1,eval);
}
}
int main() {
    run_case("v27_cauchy_uses_actual_projected_l2_and_strict_linf",[]{
        const auto r=solve();
        EXPECT_TRUE(r.accepted);
        EXPECT_TRUE(r.fields.size()==4 && r.fields[0]>1. && r.fields[3]>1.);
        EXPECT_TRUE(r.merit<1. && r.smooth_l2<std::sqrt(.625));
        EXPECT_TRUE(r.ratio>1e-3 && r.directional_l2_slope<0.);
        EXPECT_TRUE(r.bases==4 && r.evaluations>=11);
        EXPECT_TRUE(r.k_active==0 && r.omega_active==1);
    });
    run_case("v27_stale_baseline_fails_closed",[]{
        auto r=sst_cauchy_subspace_v27({1.,1.,1.,1.},{1.,0.,0.,0.},{},
            2,1e-8,1e-8,0.9,std::sqrt(.625),0,1,scalar_trial);
        EXPECT_TRUE(!r.accepted && r.evaluations==1 && r.fields.empty());
    });
    run_case("v27_never_accepts_l2_only_improvement",[]{
        const auto r=solve([](const auto& x){
            SSTV25PhysicalTrial t;
            t.fields=x;
            t.k_local=1.;t.omega_local=.5/x[3];t.merit=1.;
            t.smooth_l2=std::hypot(1.,t.omega_local)/std::sqrt(2.);
            t.k_worst=0;t.omega_worst=1;
            return t;
        });
        EXPECT_TRUE(!r.accepted && r.fields.empty());
    });
    run_case("v27_rejects_projection_outside_trust_radius",[]{
        auto ev=[](const auto& x){
            auto t=scalar_trial(x);
            if (x[0]>1.0000001) {
                t.fields[0]=x[0]*8.;
                t.k_local=.5;t.omega_local=.4;
                t.merit=.5;t.smooth_l2=.45;
            }
            return t;
        };
        auto r=solve(ev);
        EXPECT_TRUE(!r.accepted && r.fields.empty());
        EXPECT_TRUE(r.rejected_trust>0 || r.reason.find("step")!=std::string::npos);
    });
    run_case("v27_invalid_physical_certificate_fails_closed",[]{
        const auto r=solve([](const auto& x){
            auto t=scalar_trial(x);
            t.merit=std::numeric_limits<double>::quiet_NaN();
            return t;
        });
        EXPECT_TRUE(!r.accepted && r.evaluations==1);
    });
    run_case("v27_caller_state_not_mutated_even_when_projection_fails",[]{
        std::vector<double> x={1.,1.,1.,1.};
        const auto orig=x;
        auto r=sst_cauchy_subspace_v27(x,{0.,1.,0.,0.},{1.2,1.,1.,1.2},
            2,1e-8,1e-8,1.,std::sqrt(.625),0,1,
            [](const auto& z){
                if (z[0]>1.01) throw std::runtime_error("wall update refused");
                return scalar_trial(z);
            });
        EXPECT_TRUE(x==orig);
        EXPECT_TRUE(!r.accepted || r.merit<1.);
    });
    run_case("v27_invalid_floor_and_active_row_are_rejected",[]{
        EXPECT_THROW(sst_cauchy_subspace_v27({1.,1.,1.,1.},{1.,0.,0.,0.},{},
            2,1e-8,1e-8,1.,std::sqrt(.625),2,1,scalar_trial),std::invalid_argument);
        EXPECT_THROW(sst_cauchy_subspace_v27({1.,1.,1.,1.},{1.,0.,0.,0.},{},
            2,1.,1e-8,1.,std::sqrt(.625),0,1,scalar_trial),std::invalid_argument);
    });
    run_case("v27_optin_requires_v26_and_preserves_original_merit_gate",[]{
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_subspace_trust_v27_experimental);
        c.sst_subspace_trust_v27_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_smooth_globalization_v26_experimental=true;
        c.sst_log_globalization_v25_experimental=true;
        c.sst_jfnk_v22_experimental=true;
        c.sst_spatial_block_preconditioner_v24_experimental=true;
        validate_turbulence_controls(c);
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(.5,.50001,.1));
    });
    return run_all();
}
