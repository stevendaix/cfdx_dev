#include "cfdx/physics/production_solid_conduction_contract.h"
#include "common/test_harness.h"

#include <stdexcept>

using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

static CaseSetup solid_case()
{
    CaseSetup setup;
    setup.energy_model = "solid_conduction";
    setup.solver_mode = "steady";
    BoundarySpec wall;
    wall.patch_name = "wall";
    wall.type = BCType::WALL;
    wall.value_type = BCValueType::WALL_THERMAL;
    wall.wall_temperature = 300.0;
    setup.boundary_conditions.push_back(wall);
    return setup;
}

static bool rejected(const CaseSetup& setup)
{
    try {
        validate_production_solid_conduction_selection(setup);
        return false;
    } catch (const std::invalid_argument&) {
        return true;
    }
}

int main()
{
    run_case("solid_conduction_accepts_supported_physics_and_legacy_alias", [] {
        auto setup = solid_case();
        EXPECT_TRUE(!rejected(setup));
        setup.energy_model = "conjugate";
        setup.radiation_model = "off";
        setup.turbulence_model = "laminar";
        EXPECT_TRUE(!rejected(setup));
    });
    run_case("solid_conduction_does_not_silently_drop_radiation", [] {
        auto setup = solid_case();
        for (const auto* model : {"S2S", "P1", "DO", "solar", "???"}) {
            setup.radiation_model = model;
            EXPECT_TRUE(rejected(setup));
        }
    });
    run_case("solid_conduction_does_not_silently_drop_transient_evolution", [] {
        auto setup = solid_case();
        setup.transient = true;
        EXPECT_TRUE(rejected(setup));
        setup = solid_case();
        setup.solver_mode = "unsteady";
        EXPECT_TRUE(rejected(setup));
        setup.solver_mode = "???";
        EXPECT_TRUE(rejected(setup));
        setup = solid_case();
        setup.numerics.transient_scheme = "unsteady_first_order";
        EXPECT_TRUE(rejected(setup));
        setup.numerics.transient_scheme = "???";
        EXPECT_TRUE(rejected(setup));
        setup = solid_case();
        setup.numerical_config.entries.push_back({
            cfdx::core::NumericalMethodFamily::Temporal, "numerics.temporal.bdf2"});
        EXPECT_TRUE(rejected(setup));
        setup.numerical_report = cfdx::core::resolve_case_numerics(
            "solid_transient", setup.numerical_config);
        setup.numerical_config.entries.clear();
        EXPECT_TRUE(setup.numerical_report.valid());
        EXPECT_TRUE(rejected(setup));
    });
    run_case("solid_conduction_rejects_unimplemented_physics", [] {
        auto setup = solid_case();
        setup.turbulence_model = "k_omega_sst";
        EXPECT_TRUE(rejected(setup));
        setup = solid_case();
        setup.multiphase_model = "VOF";
        EXPECT_TRUE(rejected(setup));
        setup = solid_case();
        setup.physics_model = "buoyant";
        EXPECT_TRUE(rejected(setup));
        for (const auto* physics : {"compressible_euler", "compressible_navier_stokes"}) {
            setup = solid_case();
            setup.physics_model = physics;
            EXPECT_TRUE(rejected(setup));
        }
        setup = solid_case();
        setup.energy_model = "isothermal";
        EXPECT_TRUE(rejected(setup));
    });
    run_case("solid_conduction_does_not_replace_fluid_boundaries_with_adiabatic_walls", [] {
        auto setup = solid_case();
        for (const auto type : {BCType::INLET, BCType::OUTLET,
                               BCType::PRESSURE_OUTLET, BCType::INTERFACE,
                               BCType::PERIODIC, BCType::UNKNOWN}) {
            setup.boundary_conditions.front().type = type;
            EXPECT_TRUE(rejected(setup));
        }
    });
    return run_all();
}
