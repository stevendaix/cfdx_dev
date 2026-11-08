#include "cfdx/physics/production_boundary_contract.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("velocity_inlet_and_pressure_outlet_are_supported", [] {
        BoundarySpec inlet;
        inlet.patch_name="in";
        inlet.type=BCType::INLET;
        inlet.value_type=BCValueType::FIXED;
        inlet.velocity_vector={1.0,0.0,0.0};
        validate_production_incompressible_boundary(inlet,PatchType::INLET);
        BoundarySpec outlet;
        outlet.patch_name="out";
        outlet.type=BCType::PRESSURE_OUTLET;
        outlet.value_type=BCValueType::OUTLET_PRESSURE;
        outlet.pressure=123.0;
        validate_production_incompressible_boundary(outlet,PatchType::OUTLET);
    });
    run_case("zero_gradient_inlet_is_supported_as_open_boundary", [] {
        BoundarySpec inlet;
        inlet.patch_name="in";
        inlet.type=BCType::INLET;
        inlet.value_type=BCValueType::ZERO_GRADIENT;
        validate_production_incompressible_boundary(inlet,PatchType::INLET);
    });
    run_case("symmetry_requires_matching_mesh_patch", [] {
        BoundarySpec bc;
        bc.patch_name="symmetry";
        bc.type=BCType::SYMMETRY;
        bc.value_type=BCValueType::ZERO_GRADIENT;
        validate_production_incompressible_boundary(bc,PatchType::SYMMETRY);
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::OUTLET),std::invalid_argument);
    });
    run_case("pressure_outlet_cannot_silently_lose_fixed_p", [] {
        BoundarySpec bc;
        bc.patch_name="out";
        bc.type=BCType::PRESSURE_OUTLET;
        bc.value_type=BCValueType::ZERO_GRADIENT;
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::OUTLET),std::invalid_argument);
        bc.value_type=BCValueType::OUTLET_PRESSURE;
        bc.back_pressure=120.0;
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::OUTLET),std::invalid_argument);
    });
    run_case("explicit_case_conditions_resolve_generic_imported_patches", [] {
        BoundarySpec inlet;
        inlet.patch_name="generic_inlet";
        inlet.type=BCType::INLET;
        inlet.value_type=BCValueType::FIXED;
        inlet.velocity_vector={1.0,0.0,0.0};
        validate_production_incompressible_boundary(inlet,PatchType::UNKNOWN);

        BoundarySpec outlet;
        outlet.patch_name="generic_outlet";
        outlet.type=BCType::OUTLET;
        outlet.value_type=BCValueType::ZERO_GRADIENT;
        validate_production_incompressible_boundary(outlet,PatchType::UNKNOWN);

        BoundarySpec pressure_outlet;
        pressure_outlet.patch_name="generic_pressure_outlet";
        pressure_outlet.type=BCType::PRESSURE_OUTLET;
        pressure_outlet.value_type=BCValueType::OUTLET_PRESSURE;
        pressure_outlet.pressure=0.0;
        validate_production_incompressible_boundary(pressure_outlet,PatchType::UNKNOWN);

        EXPECT_THROW(validate_production_incompressible_boundary(inlet,PatchType::WALL),std::invalid_argument);
    });
    run_case("periodic_and_interface_cannot_fall_through", [] {
        BoundarySpec bc;
        bc.patch_name="paired";
        bc.value_type=BCValueType::ZERO_GRADIENT;
        bc.type=BCType::PERIODIC;
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::PERIODIC),std::invalid_argument);
        bc.type=BCType::INTERFACE;
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::INTERFACE),std::invalid_argument);
    });
    run_case("no_slip_and_slip_are_not_confused", [] {
        BoundarySpec bc;
        bc.patch_name="wall";
        bc.type=BCType::WALL;
        bc.value_type=BCValueType::WALL_NO_SLIP;
        validate_production_incompressible_boundary(bc,PatchType::WALL);
        bc.velocity_vector={1.0,0.0,0.0};
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::WALL),std::invalid_argument);
        bc.value_type=BCValueType::FIXED; // explicit moving wall
        validate_production_incompressible_boundary(bc,PatchType::WALL);
        bc.value_type=BCValueType::WALL_SLIP;
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::WALL),std::invalid_argument);
    });
    run_case("imported_pressure_and_mass_inlets_are_not_velocity_inlets", [] {
        BoundarySpec bc;
        bc.patch_name="in";
        bc.type=BCType::INLET;
        bc.value_type=BCValueType::FIXED;
        bc.source_type_name="pressure-inlet";
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::INLET),std::invalid_argument);
        bc.source_type_name="mass-flow-inlet";
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::INLET),std::invalid_argument);
    });
    run_case("empty_and_outflow_are_explicit", [] {
        BoundarySpec bc;
        bc.patch_name="front";
        bc.type=BCType::EMPTY;
        bc.value_type=BCValueType::ZERO_GRADIENT;
        validate_production_incompressible_boundary(bc,PatchType::EMPTY);
        bc.type=BCType::OUTLET;
        bc.patch_name="open";
        validate_production_incompressible_boundary(bc,PatchType::OUTLET);
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::SYMMETRY),std::invalid_argument);
    });
    run_case("duplicate_case_patch_names_are_rejected", [] {
        CaseSetup setup;
        BoundarySpec b;
        b.patch_name="inlet";
        setup.boundary_conditions={b};
        validate_production_boundary_unique_names(setup);
        setup.boundary_conditions.push_back(b);
        EXPECT_THROW(validate_production_boundary_unique_names(setup),std::invalid_argument);
    });
    run_case("bad_velocity_components_fail_before_assembly", [] {
        BoundarySpec bc;
        bc.patch_name="in";
        bc.type=BCType::INLET;
        bc.value_type=BCValueType::FIXED;
        bc.velocity_vector={1.0,2.0};
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::INLET),std::invalid_argument);
        bc.velocity_vector={1.0,std::numeric_limits<double>::quiet_NaN(),0.0};
        EXPECT_THROW(validate_production_incompressible_boundary(bc,PatchType::INLET),std::invalid_argument);
    });
    return run_all();
}