#pragma once

// Experimental coupled k-omega SST / scalar-energy seam.  This uses the
// kinematic eddy viscosity of the *accepted* momentum/SST state and a declared
// turbulent Prandtl number, not a fitted effective thermal diffusivity.
// Radiation is not part of this contract.  This is a transport closure, not
// DNS, and does not establish turbulent heat-flux validation by itself.

#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/sst_production_coupling.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

inline double require_production_sst_turbulent_prandtl(
    const cfdx::io::CaseSetup& setup)
{
    if (setup.materials.size() != 1)
        throw std::invalid_argument(
            "SST-energy: exactly one material is required");
    const auto& properties = setup.materials.front().extra_properties;
    const auto it = properties.find("turbulent_prandtl");
    if (it == properties.end())
        throw std::invalid_argument(
            "SST-energy: declare material.extra_properties.turbulent_prandtl "
            "explicitly (dimensionless; no silently assumed default)");
    const double value = parse_case_double(it->second, "turbulent_prandtl");
    if (!(value > 0.0))
        throw std::invalid_argument(
            "SST-energy: turbulent_prandtl must be strictly positive");
    return value;
}

inline std::vector<double> build_production_sst_thermal_conductivity(
    const cfdx::io::CaseSetup& setup,
    const SSTProductionState& state,
    std::size_t expected_cells,
    double turbulent_prandtl)
{
    validate_sst_production_state(state, expected_cells);
    const ProductionEnergyCase energy_case = build_production_energy_case(setup);
    // The fluid energy equation already assembles lambda_eff on each cell,
    // with face diffusion and wall heat flux reconstructed consistently.
    const auto conductivity = build_effective_thermal_conductivity(
        energy_case, state.nu_t, turbulent_prandtl);
    if (conductivity.size() != expected_cells)
        throw std::logic_error("SST-energy: conductivity shape mismatch");
    return conductivity;
}

} // namespace cfdx::physics
