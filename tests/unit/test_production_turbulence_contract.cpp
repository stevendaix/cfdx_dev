#include "cfdx/physics/production_turbulence_contract.h"
#include "common/test_harness.h"

#include <stdexcept>
#include <string>

using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("production_turbulence_accepts_laminar_selection", [] {
        CaseSetup setup;
        setup.turbulence_model = "laminar";
        validate_production_turbulence_selection(setup);

        setup.turbulence_model = "incompressible_laminar";
        validate_production_turbulence_selection(setup);

        setup.turbulence_model = "none";
        validate_production_turbulence_selection(setup);
    });

    run_case("production_turbulence_rejects_implicit_rans_physics", [] {
        CaseSetup setup;
        setup.physics_model = "incompressible_turbulent";
        bool rejected = false;
        try {
            validate_production_turbulence_selection(setup);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find("no turbulence_model") !=
                       std::string::npos;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("production_turbulence_rejects_declared_rans", [] {
        CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        bool rejected = false;
        try {
            validate_production_turbulence_selection(setup);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find("not coupled") !=
                       std::string::npos;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("production_turbulence_rejects_material_rans", [] {
        CaseSetup setup;
        MaterialSpec material;
        material.name = "air";
        material.turbulence_model = "Spalart-Allmaras";
        setup.materials.push_back(material);

        bool rejected = false;
        try {
            validate_production_turbulence_selection(setup);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find("material 'air'") !=
                       std::string::npos;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("production_turbulence_rejects_unknown_model", [] {
        CaseSetup setup;
        setup.turbulence_model = "unrecognized-model";
        bool rejected = false;
        try {
            validate_production_turbulence_selection(setup);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });

    return run_all();
}
