#include "cfdx/physics/thermal_case.h"
#include "common/test_harness.h"

using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
ThermalMaterial make_material() {
    ThermalMaterial material;
    material.name = "air";
    material.properties.conductivity.reference_value = 0.026;
    material.properties.heat_capacity.reference_value = 1005.0;
    material.properties.density.rho0 = 1.2;
    material.minimum_temperature = 250.0;
    material.maximum_temperature = 500.0;
    return material;
}

ThermalCase make_case() {
    ThermalCase c;
    c.name = "thermal_only_contract";
    c.materials.push_back(make_material());
    ThermalRegion solid;
    solid.name = "wall";
    solid.kind = ThermalRegionKind::SOLID;
    solid.mesh_id = "mesh_wall";
    solid.material = "air";
    c.regions.push_back(solid);
    ThermalRegion fluid;
    fluid.name = "air_domain";
    fluid.kind = ThermalRegionKind::FLUID;
    fluid.mesh_id = "mesh_air";
    fluid.material = "air";
    fluid.fluid_model = ThermalFluidModel::STAGNANT_DIFFUSION;
    c.regions.push_back(fluid);
    ThermalNode node;
    node.name = "air_reservoir";
    node.material = "air";
    node.volume = 0.1;
    c.nodes.push_back(node);
    return c;
}
}

int main() {
    run_case("valid_thermal_only_case", [] {
        auto c = make_case();
        c.validate();
    });

    run_case("thermal_only_fluid_does_not_require_viscosity", [] {
        auto c = make_case();
        c.materials[0].properties.viscosity.reference_value = 0.0;
        c.validate();
    });

    run_case("transient_requires_positive_time_step", [] {
        auto c = make_case();
        c.time_mode = ThermalTimeMode::TRANSIENT;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
        c.time_step = 0.01;
        c.validate();
    });

    run_case("prescribed_advection_requires_authoritative_face_flux", [] {
        auto c = make_case();
        c.regions[1].fluid_model = ThermalFluidModel::PRESCRIBED_MASS_FLUX;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
        c.regions[1].mass_flux_field = "mdot_faces";
        c.validate();
    });

    run_case("solid_rejects_fluid_transport", [] {
        auto c = make_case();
        c.regions[0].fluid_model = ThermalFluidModel::PRESCRIBED_MASS_FLUX;
        c.regions[0].mass_flux_field = "mdot_faces";
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("unknown_material_reference_is_rejected", [] {
        auto c = make_case();
        c.regions[0].material = "missing";
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("lumped_node_requires_positive_volume", [] {
        auto c = make_case();
        c.nodes[0].volume = 0.0;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("conjugate_interface_requires_solid_fluid_pair", [] {
        auto c = make_case();
        ThermalInterface interface;
        interface.name = "wall_air";
        interface.region_a = "wall";
        interface.patch_a = "wall_side";
        interface.region_b = "air_domain";
        interface.patch_b = "air_side";
        c.interfaces.push_back(interface);
        c.validate();
        c.regions[1].kind = ThermalRegionKind::SOLID;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("thermal_link_checks_endpoints_and_conductance", [] {
        auto c = make_case();
        ThermalLink link;
        link.name = "reservoir_to_wall";
        link.source.kind = ThermalEndpointKind::LUMPED_NODE;
        link.source.object = "air_reservoir";
        link.target.kind = ThermalEndpointKind::REGION_PATCH;
        link.target.object = "wall";
        link.target.patch = "inner_wall";
        link.conductance = 4.0;
        c.links.push_back(link);
        c.validate();
        c.links[0].conductance = -1.0;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("duplicate_material_names_are_rejected", [] {
        auto c = make_case();
        c.materials.push_back(make_material());
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("unsupported_enum_values_are_rejected", [] {
        auto c = make_case();
        c.time_mode = static_cast<ThermalTimeMode>(255);
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
        c = make_case();
        c.materials[0].properties.conductivity.model =
            static_cast<ScalarPropertyModel>(255);
        threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("pressure_dependent_fluid_density_has_explicit_pressure", [] {
        auto c = make_case();
        c.materials[0].properties.density.model = DensityModel::IDEAL_GAS;
        c.materials[0].properties.density.gas_constant = 287.05;
        c.regions[1].pressure = 202650.0;
        c.validate();
        c.regions[1].pressure = 0.0;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("lumped_mass_overflow_is_rejected", [] {
        auto c = make_case();
        c.nodes[0].volume = 1.7e308;
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    run_case("anisotropic_conductivity_shape_must_be_spd_and_normalized", [] {
        auto c = make_case();
        c.materials[0].anisotropic_conductivity = true;
        c.materials[0].conductivity_shape = {{1.5, 1.0, 0.5, 0.0, 0.0, 0.0}};
        c.validate();
        c.materials[0].conductivity_shape = {{1.5, 1.0, 0.5, 2.0, 0.0, 0.0}};
        bool threw = false;
        try { c.validate(); } catch (const std::invalid_argument&) { threw = true; }
        EXPECT_TRUE(threw);
    });

    return run_all();
}
