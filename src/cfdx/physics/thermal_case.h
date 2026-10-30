// Thermal-only case data contract (issue #1150, A1).
//
// This is a solver-independent typed model. It deliberately does not depend on
// Navier-Stokes state, CFD runtime state, or HDF5 implementation details.
// HDF5/YAML adapters are added in later integration work.
#pragma once

#include "cfdx/physics/thermophysical_models.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

inline constexpr std::uint32_t THERMAL_CASE_SCHEMA_VERSION = 1;

enum class ThermalTimeMode : std::uint8_t { STEADY, TRANSIENT };
enum class ThermalRegionKind : std::uint8_t { SOLID, FLUID };
enum class ThermalFluidModel : std::uint8_t {
    STAGNANT_DIFFUSION,
    PRESCRIBED_MASS_FLUX
};
enum class ThermalInterfaceModel : std::uint8_t {
    CONJUGATE_CONDUCTION,
    CONTACT_RESISTANCE
};
enum class ThermalEndpointKind : std::uint8_t { REGION_PATCH, LUMPED_NODE };
enum class ThermalLinkModel : std::uint8_t { CONDUCTANCE };

inline const char* to_string(ThermalTimeMode value) {
    return value == ThermalTimeMode::STEADY ? "steady" : "transient";
}
inline const char* to_string(ThermalRegionKind value) {
    return value == ThermalRegionKind::SOLID ? "solid" : "fluid";
}
inline const char* to_string(ThermalFluidModel value) {
    return value == ThermalFluidModel::STAGNANT_DIFFUSION
        ? "stagnant_diffusion" : "prescribed_mass_flux";
}

// All numeric properties use SI units. Temperature values are absolute kelvin.
struct ThermalMaterial {
    std::string name;
    ThermophysicalProperties properties;
    double reference_temperature = 300.0;
    double reference_pressure = 101325.0;
    double minimum_temperature = 1.0;
    double maximum_temperature = 5000.0;

    void validate() const {
        if (name.empty())
            throw std::invalid_argument("thermal material name must not be empty");
        if (!std::isfinite(reference_temperature) || reference_temperature <= 0.0 ||
            !std::isfinite(reference_pressure) || reference_pressure <= 0.0 ||
            !std::isfinite(minimum_temperature) || minimum_temperature <= 0.0 ||
            !std::isfinite(maximum_temperature) ||
            maximum_temperature < minimum_temperature ||
            reference_temperature < minimum_temperature ||
            reference_temperature > maximum_temperature)
            throw std::invalid_argument("thermal material '" + name +
                                        "' has invalid temperature/pressure range");

        // Validate the actual thermal properties used by the thermal-only
        // equations. Viscosity and turbulent Prandtl are intentionally not
        // required: a thermal-only fluid is not a momentum/turbulence model.
        const double k = properties.k(reference_temperature);
        const double cp = properties.cp(reference_temperature);
        const double rho = properties.rho(reference_temperature, reference_pressure);
        if (!(k > 0.0) || !(cp > 0.0) || !(rho > 0.0) ||
            !std::isfinite(k) || !std::isfinite(cp) || !std::isfinite(rho))
            throw std::invalid_argument("thermal material '" + name +
                                        "' must provide finite positive k, cp and rho");
        for (const double temperature : {minimum_temperature, maximum_temperature}) {
            const double edge_k = properties.k(temperature);
            const double edge_cp = properties.cp(temperature);
            const double edge_rho = properties.rho(temperature, reference_pressure);
            if (!(edge_k > 0.0) || !(edge_cp > 0.0) || !(edge_rho > 0.0) ||
                !std::isfinite(edge_k) || !std::isfinite(edge_cp) ||
                !std::isfinite(edge_rho))
                throw std::invalid_argument("thermal material '" + name +
                                            "' has invalid properties in its declared range");
        }
    }
};

struct ThermalRegion {
    std::string name;
    ThermalRegionKind kind = ThermalRegionKind::SOLID;
    std::string mesh_id;
    std::string material;
    double initial_temperature = 300.0;
    // Meaningful only for FLUID. SOLID must retain STAGNANT_DIFFUSION.
    ThermalFluidModel fluid_model = ThermalFluidModel::STAGNANT_DIFFUSION;
    // Required for PRESCRIBED_MASS_FLUX; identifies an authoritative face
    // mass-flux field. It is not reconstructed from a cell velocity field.
    std::string mass_flux_field;

    void validate() const {
        if (name.empty() || mesh_id.empty() || material.empty())
            throw std::invalid_argument("thermal region requires name, mesh_id and material");
        if (!std::isfinite(initial_temperature) || initial_temperature <= 0.0)
            throw std::invalid_argument("thermal region '" + name +
                                        "' has invalid initial temperature");
        if (kind == ThermalRegionKind::SOLID) {
            if (fluid_model != ThermalFluidModel::STAGNANT_DIFFUSION ||
                !mass_flux_field.empty())
                throw std::invalid_argument("solid region '" + name +
                                            "' cannot declare a fluid transport model");
        } else if (fluid_model == ThermalFluidModel::PRESCRIBED_MASS_FLUX) {
            if (mass_flux_field.empty())
                throw std::invalid_argument("fluid region '" + name +
                                            "' requires an explicit face mass-flux field");
        } else if (!mass_flux_field.empty()) {
            throw std::invalid_argument("stagnant fluid region '" + name +
                                        "' must not declare a mass-flux field");
        }
    }
};

struct ThermalNode {
    std::string name;
    std::string material;
    double volume = 0.0; // m^3; mass is derived from rho(T,p) * volume.
    double initial_temperature = 300.0;
    double pressure = 101325.0;

    void validate() const {
        if (name.empty() || material.empty())
            throw std::invalid_argument("thermal node requires name and material");
        if (!std::isfinite(volume) || volume <= 0.0 ||
            !std::isfinite(initial_temperature) || initial_temperature <= 0.0 ||
            !std::isfinite(pressure) || pressure <= 0.0)
            throw std::invalid_argument("thermal node '" + name +
                                        "' has invalid volume, temperature or pressure");
    }
};

struct ThermalInterface {
    std::string name;
    std::string region_a;
    std::string patch_a;
    std::string region_b;
    std::string patch_b;
    ThermalInterfaceModel model = ThermalInterfaceModel::CONJUGATE_CONDUCTION;
    // m^2.K/W; used only by CONTACT_RESISTANCE.
    double contact_resistance = 0.0;

    void validate() const {
        if (name.empty() || region_a.empty() || region_b.empty() ||
            patch_a.empty() || patch_b.empty() || region_a == region_b)
            throw std::invalid_argument("thermal interface requires two distinct regions and patches");
        if (!std::isfinite(contact_resistance) || contact_resistance < 0.0)
            throw std::invalid_argument("thermal interface '" + name +
                                        "' has invalid contact resistance");
        if (model == ThermalInterfaceModel::CONJUGATE_CONDUCTION &&
            contact_resistance != 0.0)
            throw std::invalid_argument("conjugate interface '" + name +
                                        "' cannot also specify contact resistance");
    }
};

struct ThermalEndpoint {
    ThermalEndpointKind kind = ThermalEndpointKind::LUMPED_NODE;
    std::string object;
    // Required for REGION_PATCH and forbidden for LUMPED_NODE.
    std::string patch;
};

struct ThermalLink {
    std::string name;
    ThermalEndpoint source;
    ThermalEndpoint target;
    ThermalLinkModel model = ThermalLinkModel::CONDUCTANCE;
    double conductance = 0.0; // W/K; heat rate is G * (T_source - T_target).

    void validate() const {
        if (name.empty() || source.object.empty() || target.object.empty())
            throw std::invalid_argument("thermal link requires a name and two endpoints");
        if (source.kind == target.kind && source.object == target.object &&
            source.patch == target.patch)
            throw std::invalid_argument("thermal link '" + name + "' is self-connected");
        for (const ThermalEndpoint* endpoint : {&source, &target}) {
            if (endpoint->kind == ThermalEndpointKind::REGION_PATCH) {
                if (endpoint->patch.empty())
                    throw std::invalid_argument("thermal link '" + name +
                                                "' region endpoint requires a patch");
            } else if (!endpoint->patch.empty()) {
                throw std::invalid_argument("thermal link '" + name +
                                            "' lumped-node endpoint cannot specify a patch");
            }
        }
        if (!std::isfinite(conductance) || conductance <= 0.0)
            throw std::invalid_argument("thermal link '" + name +
                                        "' conductance must be finite and positive");
    }
};

struct ThermalCase {
    std::uint32_t schema_version = THERMAL_CASE_SCHEMA_VERSION;
    std::string name;
    ThermalTimeMode time_mode = ThermalTimeMode::STEADY;
    double time_step = 0.0; // s; required positive only for TRANSIENT.
    std::vector<ThermalMaterial> materials;
    std::vector<ThermalRegion> regions;
    std::vector<ThermalNode> nodes;
    std::vector<ThermalInterface> interfaces;
    std::vector<ThermalLink> links;

    void validate() const {
        if (schema_version != THERMAL_CASE_SCHEMA_VERSION)
            throw std::invalid_argument("unsupported thermal case schema version");
        if (name.empty())
            throw std::invalid_argument("thermal case name must not be empty");
        if (!std::isfinite(time_step) || time_step < 0.0 ||
            (time_mode == ThermalTimeMode::TRANSIENT && time_step <= 0.0) ||
            (time_mode == ThermalTimeMode::STEADY && time_step != 0.0))
            throw std::invalid_argument("thermal case time_step is inconsistent with time_mode");
        if (materials.empty() || (regions.empty() && nodes.empty()))
            throw std::invalid_argument("thermal case requires materials and at least one region or node");

        std::map<std::string, const ThermalMaterial*> material_by_name;
        for (const auto& material : materials) {
            material.validate();
            if (!material_by_name.emplace(material.name, &material).second)
                throw std::invalid_argument("duplicate thermal material '" + material.name + "'");
        }

        std::set<std::string> region_names;
        for (const auto& region : regions) {
            region.validate();
            if (!region_names.insert(region.name).second)
                throw std::invalid_argument("duplicate thermal region '" + region.name + "'");
            if (material_by_name.count(region.material) == 0)
                throw std::invalid_argument("thermal region '" + region.name +
                                            "' references unknown material '" + region.material + "'");
            const auto& material = *material_by_name.at(region.material);
            if (region.initial_temperature < material.minimum_temperature ||
                region.initial_temperature > material.maximum_temperature)
                throw std::invalid_argument("thermal region '" + region.name +
                                            "' initial temperature is outside material range");
        }

        std::set<std::string> node_names;
        for (const auto& node : nodes) {
            node.validate();
            if (!node_names.insert(node.name).second)
                throw std::invalid_argument("duplicate thermal node '" + node.name + "'");
            if (material_by_name.count(node.material) == 0)
                throw std::invalid_argument("thermal node '" + node.name +
                                            "' references unknown material '" + node.material + "'");
            const auto& material = *material_by_name.at(node.material);
            if (node.initial_temperature < material.minimum_temperature ||
                node.initial_temperature > material.maximum_temperature)
                throw std::invalid_argument("thermal node '" + node.name +
                                            "' initial temperature is outside material range");
            (void)material.properties.rho(node.initial_temperature, node.pressure);
        }

        std::set<std::string> interface_names;
        for (const auto& interface : interfaces) {
            interface.validate();
            if (!interface_names.insert(interface.name).second)
                throw std::invalid_argument("duplicate thermal interface '" + interface.name + "'");
            const auto a = std::find_if(regions.begin(), regions.end(),
                [&](const ThermalRegion& r) { return r.name == interface.region_a; });
            const auto b = std::find_if(regions.begin(), regions.end(),
                [&](const ThermalRegion& r) { return r.name == interface.region_b; });
            if (a == regions.end() || b == regions.end())
                throw std::invalid_argument("thermal interface '" + interface.name +
                                            "' references an unknown region");
            if (interface.model == ThermalInterfaceModel::CONJUGATE_CONDUCTION &&
                a->kind == b->kind)
                throw std::invalid_argument("conjugate interface '" + interface.name +
                                            "' must connect one solid and one fluid region");
        }

        std::set<std::string> link_names;
        for (const auto& link : links) {
            link.validate();
            if (!link_names.insert(link.name).second)
                throw std::invalid_argument("duplicate thermal link '" + link.name + "'");
            for (const ThermalEndpoint* endpoint : {&link.source, &link.target}) {
                if (endpoint->kind == ThermalEndpointKind::REGION_PATCH) {
                    const auto region = std::find_if(regions.begin(), regions.end(),
                        [&](const ThermalRegion& r) { return r.name == endpoint->object; });
                    if (region == regions.end())
                        throw std::invalid_argument("thermal link '" + link.name +
                                                    "' references an unknown region");
                } else if (node_names.count(endpoint->object) == 0) {
                    throw std::invalid_argument("thermal link '" + link.name +
                                                "' references an unknown lumped node");
                }
            }
        }
    }
};

} // namespace cfdx::physics
