#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/physics/production_turbulence_contract.h"
#include "cfdx/physics/production_time_contract.h"

#include <stdexcept>
#include <string>

namespace cfdx::physics {

// The solid executable solves steady conduction only. A declaration must not
// be silently discarded merely because the conduction subsystem converges.
// "conjugate" is retained as the documented legacy alias for solid-only cases;
// it does not authorize a fluid/solid or radiation coupling in this executable.
inline void validate_production_solid_conduction_selection(
    const cfdx::io::CaseSetup& setup)
{
    if (setup.energy_model != "solid_conduction" &&
        setup.energy_model != "conjugate")
        throw std::invalid_argument(
            "thermal solver: energy_model must be 'solid_conduction' or the "
            "solid-only legacy alias 'conjugate'");

    const auto inactive = [](const std::string& value) {
        const auto key = normalize_turbulence_model_key(value);
        return value.find_first_not_of(" \t\r\n") == std::string::npos ||
            key == "none" || key == "off" || key == "disabled";
    };
    if (!inactive(setup.radiation_model))
        throw std::invalid_argument(
            "thermal solver: declared radiation is not coupled to the solid "
            "conduction executable; use a supported radiation-energy path");
    if (!inactive(setup.multiphase_model))
        throw std::invalid_argument(
            "thermal solver: multiphase physics is not supported by steady solid conduction");
    validate_production_turbulence_selection(setup);

    validate_production_steady_time_selection(setup, "thermal solver");

    const auto physics = normalize_turbulence_model_key(setup.physics_model);
    if (physics == "buoyant" || physics.rfind("compressible", 0) == 0 ||
        physics == "euler")
        throw std::invalid_argument(
            "thermal solver: declared flow physics is not solved by steady solid conduction");

    for (const auto& boundary : setup.boundary_conditions) {
        using cfdx::io::BCType;
        if (boundary.type != BCType::WALL &&
            boundary.type != BCType::SYMMETRY && boundary.type != BCType::EMPTY)
            throw std::invalid_argument(
                "thermal solver: patch '" + boundary.patch_name +
                "' declares flow, coupling, periodic or unknown boundary physics "
                "that steady solid conduction does not implement");
    }
}

}  // namespace cfdx::physics
