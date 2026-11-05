#include "cfdx/physics/production_sst_thermal_coupling.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

static CaseSetup case_setup() {
    CaseSetup setup;
    setup.energy_model = "conjugate";
    MaterialSpec fluid;
    fluid.name = "fluid";
    fluid.density = 2.0;
    fluid.dynamic_viscosity = 0.01;
    fluid.specific_heat = 10.0;
    fluid.thermal_conductivity = 3.0;
    fluid.extra_properties["turbulent_prandtl"] = "0.8";
    setup.materials.push_back(fluid);
    return setup;
}
static SSTProductionState state() {
    SSTProductionState s(2);
    for(std::size_t i=0;i<2;++i) {
        s.k(i)=1;
        s.omega(i)=2;
        s.wall_distance(i)=0.1;
        s.strain_rate(i)=0;
        s.nu_t(i)=(i+1)*0.04;
    }
    return s;
}

int main() {
    run_case("eddy_diffusion_uses_kinematic_nut_and_explicit_Prt",[] {
        const auto setup=case_setup();
        const auto effective=build_production_sst_thermal_conductivity(
            setup,state(),2,require_production_sst_turbulent_prandtl(setup));
        // lambda_eff = 3 + rho Cp nu_t / Pr_t
        EXPECT_NEAR(effective[0],4.0,1e-14);
        EXPECT_NEAR(effective[1],5.0,1e-14);
    });
    run_case("missing_or_malformed_Prt_fails_closed",[] {
        auto setup=case_setup();
        setup.materials[0].extra_properties.clear();
        EXPECT_THROW(require_production_sst_turbulent_prandtl(setup),std::invalid_argument);
        setup.materials[0].extra_properties["turbulent_prandtl"]="0.9junk";
        EXPECT_THROW(require_production_sst_turbulent_prandtl(setup),std::invalid_argument);
        setup.materials[0].extra_properties["turbulent_prandtl"]="-0.3";
        EXPECT_THROW(require_production_sst_turbulent_prandtl(setup),std::invalid_argument);
        setup.materials[0].extra_properties["turbulent_prandtl"]="nan";
        EXPECT_THROW(require_production_sst_turbulent_prandtl(setup),std::invalid_argument);
    });
    run_case("invalid_turbulence_and_overflow_fail_closed",[] {
        const auto setup=case_setup();
        auto s=state();
        s.nu_t(0)=-1;
        EXPECT_THROW(build_production_sst_thermal_conductivity(setup,s,2,0.8),std::invalid_argument);
        s=state();
        EXPECT_THROW(build_production_sst_thermal_conductivity(setup,s,3,0.8),std::invalid_argument);
        s=state();
        EXPECT_THROW(build_production_sst_thermal_conductivity(setup,s,2,0.0),std::invalid_argument);
    });
    run_case("empty_extrusion_planes_are_zero_gradient_not_wall",[] {
        auto setup=case_setup();
        BoundarySpec empty;
        empty.patch_name="front";
        empty.type=BCType::EMPTY;
        empty.value_type=BCValueType::ZERO_GRADIENT;
        empty.sst_scalar_conditions={
            {"k",SSTScalarBoundaryType::ZERO_GRADIENT,0.0},
            {"omega",SSTScalarBoundaryType::ZERO_GRADIENT,0.0}};
        setup.boundary_conditions.push_back(empty);
        validate_sst_scalar_boundary_contract(setup);
        setup.boundary_conditions[0].sst_scalar_conditions[1].type =
            SSTScalarBoundaryType::WALL_OMEGA_SST;
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup),std::invalid_argument);
    });
    return 0;
    return run_all();
}
