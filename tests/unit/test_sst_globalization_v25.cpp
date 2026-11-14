#include "cfdx/physics/sst_globalization_v25.h"
#include "cfdx/physics/turbulence_transport.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

static SSTV25PhysicalTrial linear_merit(const std::vector<double>& v) {
    SSTV25PhysicalTrial t;
    t.fields=v;
    t.k_local=std::abs(v[0]-2.0);
    t.omega_local=std::abs(v[1]-3.0);
    t.merit=std::max(t.k_local,t.omega_local);
    t.k_worst=0;t.omega_worst=0;
    return t;
}

int main() {
    run_case("v25_logarithmic_full_newton_step_reaches_coupled_physical_root",[] {
        const auto initial=std::vector<double>{1.,1.};
        const auto r=sst_globalize_log_direction_v25(
            initial,{std::log(2.),std::log(3.)},1,1e-12,1e-12,
            2.,linear_merit);
        EXPECT_TRUE(r.accepted);
        EXPECT_NEAR(r.merit,0.,1e-12);
        EXPECT_NEAR(r.fields[0],2.,1e-12);
        EXPECT_NEAR(r.fields[1],3.,1e-12);
        EXPECT_TRUE(initial[0]==1. && initial[1]==1.);
        EXPECT_TRUE(r.valid==13 && r.improving>0);
    });
    run_case("v25_rejects_every_non_descent_direction_and_never_mutates_base",[] {
        const auto initial=std::vector<double>{1.,1.};
        const auto r=sst_globalize_log_direction_v25(
            initial,{0.5,0.5},1,1e-12,1e-12,1.,
            [](const auto& x){
                SSTV25PhysicalTrial t;
                t.fields=x;t.k_local=x[0];t.omega_local=x[1];
                t.merit=std::max(t.k_local,t.omega_local);
                return t;
            });
        EXPECT_TRUE(!r.accepted);
        EXPECT_TRUE(r.fields.empty());
        EXPECT_TRUE(r.valid==13);
        EXPECT_TRUE(r.best_trial_merit>=1.);
        EXPECT_TRUE(initial[0]==1. && initial[1]==1.);
    });
    run_case("v25_log_trust_region_backtracks_without_clipping_physical_fields",[] {
        const auto r=sst_globalize_log_direction_v25(
            {1.,1.},{100.,0.},1,1e-12,1e-12,0.5,
            [](const auto& x){
                SSTV25PhysicalTrial t;
                t.fields=x;t.k_local=std::abs(x[0]-1.5);
                t.omega_local=0.;t.merit=t.k_local;
                return t;
            });
        EXPECT_TRUE(r.accepted);
        EXPECT_TRUE(r.alpha<=0.015);
        EXPECT_TRUE(r.fields[0]>1.);
        EXPECT_TRUE(r.merit<0.5);
        EXPECT_NEAR(r.maximum_absolute_log_direction,100.,0.);
    });
    run_case("v25_owner_block_uses_true_wall_target_derivative",[] {
        const SSTBlockPreconditionerV23 approx({{8.,0.8,0.,8.}});
        const auto physical=[](const auto& x){
            return std::vector<double>{2.*x[0]+0.1*x[1],
                                       x[1]-3.*std::sqrt(x[0])};
        };
        const auto corrected=sst_correct_wall_blocks_v25(
            approx,{0},{4.,8.},{1.,1.},{1.,1.},physical);
        EXPECT_TRUE(corrected.corrected==1);
        EXPECT_TRUE(corrected.max_relative_error>0.3);
        EXPECT_NEAR(corrected.preconditioner.raw_blocks()[0][0],8.,1e-7);
        EXPECT_NEAR(corrected.preconditioner.raw_blocks()[0][1],0.8,1e-7);
        EXPECT_NEAR(corrected.preconditioner.raw_blocks()[0][2],-3.,1e-7);
        EXPECT_NEAR(corrected.preconditioner.raw_blocks()[0][3],8.,1e-7);
        EXPECT_TRUE(approx.raw_blocks()[0][2]==0.);
        EXPECT_TRUE(corrected.preconditioner.diagonal_fallbacks()==0);
    });
    run_case("v25_singular_wall_derivative_is_not_injected",[] {
        const SSTBlockPreconditionerV23 approx({{2.,0.,0.,3.}});
        const auto result=sst_correct_wall_blocks_v25(
            approx,{0},{2.,3.},{1.,1.},{1.,1.},
            [](const auto& x){return std::vector<double>{x[0]+x[1],2*x[0]+2*x[1]};});
        EXPECT_TRUE(result.corrected==0);
        for(std::size_t j=0;j<4;++j) EXPECT_NEAR(result.preconditioner.raw_blocks()[0][j],
                                        approx.raw_blocks()[0][j],0.);
    });
    run_case("v25_rejects_invalid_projection_and_duplicate_wall_owners",[] {
        const auto x=std::vector<double>{1.,1.};
        const auto r=sst_globalize_log_direction_v25(x,{-0.1,0.},1,1e-12,1e-12,1.,
            [](const auto&){
                SSTV25PhysicalTrial t;
                t.fields={1.,std::numeric_limits<double>::quiet_NaN()};
                t.k_local=0.;t.omega_local=0.;t.merit=0.;return t;
            });
        EXPECT_TRUE(!r.accepted && r.valid==0 && r.fields.empty());
        const SSTBlockPreconditionerV23 b({{1.,0.,0.,1.}});
        EXPECT_THROW(sst_correct_wall_blocks_v25(b,{0,0},{1.,1.},{1.,1.},
            {1.,1.},[](const auto& x){return x;}),std::invalid_argument);
        EXPECT_THROW(sst_globalize_log_direction_v25({1.,1.},{0.,0.},1,1.,1e-12,
            1.,linear_merit),std::invalid_argument);
    });
    run_case("v25_option_is_disabled_and_dependency_checked",[] {
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_log_globalization_v25_experimental);
        c.sst_log_globalization_v25_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_jfnk_v22_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_spatial_block_preconditioner_v24_experimental=true;
        validate_turbulence_controls(c);
    });
    return run_all();
}
