#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_globalization_v26.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("v26_rms_uses_same_native_and_compensated_FV_rows",[] {
        ScalarEquation eq;
        cfdx::core::SparseMatrix a(2,2);
        a.push_back(0,0,1.);a.push_back(1,1,1.);a.finalize();
        eq.matrix=a;eq.rhs=cfdx::core::Vector(2,0.);
        eq.rhs(1)=1.;eq.diagonal={1.,1.};
        cfdx::core::Vector x(2,2.);
        const auto native=sst_local_equation_defect(eq,x,1.e-8);
        const auto measured=sst_local_residual_precision_audit(eq,x,1.e-8);
        EXPECT_NEAR(native.maximum_relative,1.,1e-14);
        EXPECT_NEAR(native.rms_relative,std::sqrt(.625),1e-14);
        EXPECT_NEAR(measured.compensated.maximum_relative,1.,1e-14);
        EXPECT_NEAR(measured.compensated.rms_relative,std::sqrt(.625),1e-14);
        SSTNonlinearResidualSample sample;
        sample.k_local_l2=.3;sample.omega_local_l2=.4;
        EXPECT_NEAR(sst_nonlinear_smooth_merit_v26(sample),std::sqrt(.125),1e-14);
        sample.k_local_l2=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(sst_nonlinear_smooth_merit_v26(sample),std::invalid_argument);
    });
    run_case("v26_log_jvp_seven_evaluations_and_three_amplitudes",[] {
        const auto r=sst_audit_log_jvp_v26(
            {2.,3.},{1.,-2.},1,1e-8,1e-8,
            [](const auto& x){return std::vector<double>{3.*std::log(x[0])+2.*std::log(x[1]),
                                                           std::log(x[0])-std::log(x[1])};});
        EXPECT_TRUE(r.available && r.evaluations==7);
        EXPECT_TRUE(r.steps[0]>r.steps[1] && r.steps[1]>r.steps[2]);
        EXPECT_TRUE(r.maximum_disagreement<1e-6);
        EXPECT_NEAR(r.derivative_norms[0],std::sqrt(2.5),1e-9);
    });
    run_case("v26_log_jvp_detects_nonsmooth_owner_source",[] {
        const auto r=sst_audit_log_jvp_v26(
            {1.,1.},{1.,1.},1,1e-8,1e-8,
            [](const auto& x){return std::vector<double>{std::abs(std::log(x[0])),
                                                           std::log(x[1])};});
        EXPECT_TRUE(r.available);
        EXPECT_TRUE(r.maximum_disagreement>1.);
    });
    run_case("v26_l2_and_linf_proposals_share_single_set_of_FV_trials",[] {
        std::size_t calls=0;
        const auto x=std::vector<double>{1.,1.};
        const auto r=sst_globalize_log_direction_v26(
            x,{std::log(2.),std::log(2.)},1,1e-8,1e-8,1.,.9,
            [&](const auto& candidate){
                ++calls;
                SSTV25PhysicalTrial t;
                t.fields=candidate;
                if (candidate[0]>1.8) {
                    t.k_local=.5;t.omega_local=.5;
                } else if (candidate[0]>1.3) {
                    t.k_local=.6;t.omega_local=0.;
                } else {
                    t.k_local=.95;t.omega_local=.8;
                }
                t.merit=std::max(t.k_local,t.omega_local);
                t.smooth_l2=std::hypot(t.k_local,t.omega_local)/std::sqrt(2.);
                return t;
            });
        EXPECT_TRUE(calls==13 && r.linf.valid==13);
        EXPECT_TRUE(r.linf.accepted && r.smooth_accepted);
        EXPECT_NEAR(r.linf.alpha,1.,1e-14);
        EXPECT_NEAR(r.smooth_alpha,.5,1e-14);
        EXPECT_NEAR(r.linf.merit,.5,1e-14);
        EXPECT_NEAR(r.selected_linf,.6,1e-14);
        EXPECT_TRUE(r.selected_l2<.5);
        EXPECT_TRUE(x[0]==1. && x[1]==1.);
    });
    run_case("v26_rejects_all_L2_only_improvements_when_physical_Linf_worsens",[] {
        const auto r=sst_globalize_log_direction_v26(
            {1.,1.},{0.1,0.1},1,1e-10,1e-10,1.,.8,
            [](const auto& candidate){
                SSTV25PhysicalTrial t;
                t.fields=candidate;
                t.k_local=1.1;t.omega_local=.1;t.merit=1.1;
                t.smooth_l2=.1;
                return t;
            });
        EXPECT_TRUE(!r.linf.accepted && !r.smooth_accepted);
        EXPECT_TRUE(r.linf.improving==0 && r.smooth_improving==0);
        EXPECT_TRUE(r.linf.fields.empty() && r.smooth_fields.empty());
    });
    run_case("v26_rejects_invalid_projection_and_stale_smooth_metrics",[] {
        const auto r=sst_globalize_log_direction_v26(
            {1.,1.},{-.1,-.1},1,1e-10,1e-10,1.,.8,
            [](const auto& candidate){
                SSTV25PhysicalTrial t;
                t.fields=candidate;
                t.k_local=.2;t.omega_local=.2;t.merit=.2;
                t.smooth_l2=std::numeric_limits<double>::quiet_NaN();
                return t;
            });
        EXPECT_TRUE(!r.linf.accepted && !r.smooth_accepted && r.linf.valid==0);
        EXPECT_THROW(sst_audit_log_jvp_v26({1.,1.},{0.,0.},1,1.,1e-8,
            [](const auto& x){return x;}),std::invalid_argument);
    });
    run_case("v26_opt_in_requires_v25_v24_v22_and_never_changes_threshold",[] {
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_smooth_globalization_v26_experimental);
        c.sst_smooth_globalization_v26_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_jfnk_v22_experimental=true;
        c.sst_spatial_block_preconditioner_v24_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_log_globalization_v25_experimental=true;
        validate_turbulence_controls(c);
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(.5,.51,.5));
        EXPECT_TRUE(sst_nonlinear_step_reduces_merit(.5,.49,.5));
    });
    return run_all();
}
