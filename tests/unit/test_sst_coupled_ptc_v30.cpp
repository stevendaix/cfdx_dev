#include "cfdx/physics/sst_coupled_pseudo_time_v30.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("v30_pseudo_mass_uses_transport_scale_and_skips_wall_row",[]{
        const auto mass=sst_v30_coupled_pseudo_mass(
            {4.,5.,1.,7.},{5.,8.,2.,14.},{12.,1.,2.,9.},{0,0,1,0},0.5);
        EXPECT_NEAR(mass[0],12.,1e-14); // max(4,2.5,6)/0.5
        EXPECT_NEAR(mass[1],10.,1e-14);
        EXPECT_NEAR(mass[2],0.,0.);
        EXPECT_NEAR(mass[3],14.,1e-14);
    });
    run_case("v30_shift_preserves_old_physical_residual",[]{
        const std::vector<double> base{2.,-3.},old{1.,2.},mass{4.,5.};
        const auto result=sst_v30_shift_residual(base,old,old,mass);
        EXPECT_NEAR(result[0],2.,0.);
        EXPECT_NEAR(result[1],-3.,0.);
        const auto moved=sst_v30_shift_residual(base,{1.25,1.5},old,mass);
        EXPECT_NEAR(moved[0],3.,1e-14);
        EXPECT_NEAR(moved[1],-5.5,1e-14);
        EXPECT_NEAR(base[0],2.,0.);
    });
    run_case("v30_continuation_updates_after_true_accepted_step_only",[]{
        SSTV30ContinuationState c;
        EXPECT_NEAR(c.cfl,0.5,0.);
        c.accept_inner(1.,0.8);
        EXPECT_NEAR(c.cfl,0.6,1e-14);
        c.accept_inner(1.,0.999);
        EXPECT_NEAR(c.cfl,0.63,1e-14);
        EXPECT_THROW(c.accept_inner(1.,1.),std::invalid_argument);
        EXPECT_NEAR(c.cfl,0.63,1e-14);
        c.reject_inner();
        EXPECT_NEAR(c.cfl,0.315,1e-14);
        EXPECT_TRUE(c.accepted_inner==2 && c.rejected_inner==1);
        c.outer_commit();
        EXPECT_TRUE(c.accepted_outer==1);
    });
    run_case("v30_augmented_row_scaling_matches_shifted_log_jacobian",[]{
        const auto s=sst_v30_augmented_row_scaling(
            {2.,3.,4.},{3.,0.,6.},{2.,4.,1.},{2.,4.,1.});
        EXPECT_NEAR(s[0],5.,1e-14);
        EXPECT_NEAR(s[1],3.,1e-14);
        EXPECT_NEAR(s[2],10.,1e-14);
        EXPECT_THROW(sst_v30_augmented_row_scaling(
            {1.},{-1.},{1.},{1.}),std::invalid_argument);
        EXPECT_THROW(sst_v30_augmented_row_scaling(
            {1.},{1.},{0.},{1.}),std::invalid_argument);
    });
    run_case("v30_invalid_coefficient_is_rejected_without_mutation",[]{
        const double inf=std::numeric_limits<double>::infinity();
        EXPECT_THROW(sst_v30_coupled_pseudo_mass({2.},{2.},{inf},{},1.),
                     std::invalid_argument);
        EXPECT_THROW(sst_v30_coupled_pseudo_mass({2.},{2.},{0.},{},0.),
                     std::invalid_argument);
        EXPECT_THROW(sst_v30_shift_residual({1.},{3.},{1.},{inf}),
                     std::invalid_argument);
        SSTV30ContinuationState c;
        c.cfl=inf;
        EXPECT_THROW(c.reject_inner(),std::invalid_argument);
    });
    run_case("v30_outer_rollback_restores_persistent_controller",[]{
        SSTProductionRuntime runtime;
        runtime.transport.model=TurbulenceModel::SST;
        runtime.transport.sst_jfnk_v22_experimental=true;
        runtime.transport.sst_coupled_ptc_v30_experimental=true;
        runtime.transport.sst_continuation_state_v30=&runtime.coupled_ptc_v30;
        validate_turbulence_controls(runtime.transport);
        const auto snapshot=runtime.coupled_ptc_v30;
        runtime.coupled_ptc_v30.accept_inner(1.,0.8);
        runtime.coupled_ptc_v30=snapshot;
        EXPECT_NEAR(runtime.transport.sst_continuation_state_v30->cfl,0.5,0.);
        EXPECT_TRUE(runtime.coupled_ptc_v30.accepted_inner==0);
    });
    run_case("v30_disabled_by_default_and_requires_jfnk",[]{
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_coupled_ptc_v30_experimental);
        c.sst_coupled_ptc_v30_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_jfnk_v22_experimental=true;
        validate_turbulence_controls(c);
    });
    return run_all();
}
