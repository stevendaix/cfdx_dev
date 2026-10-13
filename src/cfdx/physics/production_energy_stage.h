#pragma once

// T2 production energy stage: frozen-flow coupling of the incompressible
// solver's authoritative conservative face mass flux with the steady energy
// equation (capability #885 stage 2, issues #880 / #118).
//
// Contract:
//   * energy_model gate: {"", "none", "isothermal", "solid_conduction",
//     "conjugate"} are accepted declarations; anything else - including
//     "buoyant" - is rejected fail-closed instead of silently degrading to
//     an isothermal run. The stage itself only runs for "conjugate":
//     calling it under any other declared model throws.
//   * Exactly one fluid material with finite positive density, conductivity
//     and specific heat. Volumetric power (W/m^3) comes from the material's
//     extra_properties ("volumetric_power", default 0).
//   * Temperature boundary conditions from the case: WALL_THERMAL with a
//     non-zero wall_temperature -> Dirichlet; INLET or FIXED with a non-zero
//     temperature -> Dirichlet; otherwise raw_params["heat_flux"]
//     (W/m^2, positive into the domain) -> fixed gradient q''/k; otherwise
//     adiabatic. Zero means "not declared" (stage-1 / T1 convention). At
//     least one Dirichlet patch is required: an all-Neumann energy system
//     has a constant null space and is rejected.
//   * The stage consumes the flow result's authoritative conservative face
//     mass flux, never a flux reconstructed from the cell velocity.
//   * Conservation audits are independent of the assembled matrix: the
//     energy flux is rebuilt from the final temperature and the same mass
//     flux (reconstruct_energy_balance), and per-patch fluxes are rebuilt
//     the same way.

#include "cfdx/core/field/field.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/finite_volume_transport.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct ProductionEnergyCase {
    double density = 0.0;          // kg/m^3
    double conductivity = 0.0;     // W/(m.K)
    double specific_heat = 0.0;    // J/(kg.K)
    double volumetric_power = 0.0; // W/m^3
    EnergySolverControls controls;
};

inline bool production_energy_requested(const cfdx::io::CaseSetup& setup)
{
    return setup.energy_model == "conjugate";
}

inline void validate_production_energy_model(const cfdx::io::CaseSetup& setup)
{
    const std::string& model = setup.energy_model;
    if (model.empty() || model == "none" || model == "isothermal" ||
        model == "solid_conduction" || model == "conjugate")
        return;
    throw std::invalid_argument(
        "production energy stage: unsupported energy_model '" + model +
        "' (accepted: none, isothermal, solid_conduction, conjugate; "
        "buoyant coupling is not implemented and is rejected instead of "
        "silently running an isothermal solve)");
}

inline ProductionEnergyCase build_production_energy_case(
    const cfdx::io::CaseSetup& setup)
{
    if (setup.materials.size() != 1)
        throw std::invalid_argument(
            "production energy stage requires exactly one fluid material "
            "(declared: " + std::to_string(setup.materials.size()) + ")");
    const auto& material = setup.materials.front();
    if (!std::isfinite(material.density) || !(material.density > 0.0))
        throw std::invalid_argument(
            "production energy stage: material density must be finite and "
            "positive");
    if (!std::isfinite(material.thermal_conductivity) ||
        !(material.thermal_conductivity > 0.0))
        throw std::invalid_argument(
            "production energy stage: material thermal_conductivity must be "
            "finite and positive");
    if (!std::isfinite(material.specific_heat) ||
        !(material.specific_heat > 0.0))
        throw std::invalid_argument(
            "production energy stage: material specific_heat must be finite "
            "and positive");

    double volumetric_power = 0.0;
    const auto power_it = material.extra_properties.find("volumetric_power");
    if (power_it != material.extra_properties.end()) {
        try {
            volumetric_power = std::stod(power_it->second);
        } catch (const std::exception&) {
            throw std::invalid_argument(
                "production energy stage: volumetric_power must be numeric: '" +
                power_it->second + "'");
        }
        if (!std::isfinite(volumetric_power))
            throw std::invalid_argument(
                "production energy stage: volumetric_power must be finite");
    }

    ProductionEnergyCase energy_case;
    energy_case.density = material.density;
    energy_case.conductivity = material.thermal_conductivity;
    energy_case.specific_heat = material.specific_heat;
    energy_case.volumetric_power = volumetric_power;
    energy_case.controls.density = material.density;
    energy_case.controls.cp = material.specific_heat;
    energy_case.controls.conductivity = material.thermal_conductivity;
    // Frozen-flow steady stage: no transient term, and relaxation stays at 1
    // because the assembled equation is linear for a fixed source, flux and
    // conductivity - a fully converged linear solve is then the exact fixed
    // point, and under-relaxation would only blur that identity.
    energy_case.controls.dt = 0.0;
    energy_case.controls.relaxation = 1.0;
    energy_case.controls.max_iterations = 500;
    energy_case.controls.tolerance = 1e-12;
    energy_case.controls.linear_tolerance = 1e-14;
    energy_case.controls.temperature_tolerance = 1e-12;
    energy_case.controls.energy_balance_tolerance = 1e-12;
    validate_energy_controls(energy_case.controls);
    return energy_case;
}

inline ScalarBoundaryConditions build_production_temperature_bcs(
    const cfdx::core::Mesh& mesh,
    const cfdx::io::CaseSetup& setup,
    const ProductionEnergyCase& energy_case)
{
    if (!std::isfinite(energy_case.conductivity) ||
        !(energy_case.conductivity > 0.0))
        throw std::invalid_argument(
            "build_production_temperature_bcs: invalid conductivity");

    ScalarBoundaryConditions bcs;
    std::size_t dirichlet_patches = 0;
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
        const auto& patch = mesh.boundary().patch(i);
        if (patch.type == cfdx::core::PatchType::EMPTY)
            continue;
        ScalarBoundaryCondition bc;
        bc.type = ScalarBoundaryType::ZERO_GRADIENT;  // adiabatic by default
        const auto* spec = setup.find_boundary(patch.name);
        bool fixed_value = false;
        if (spec != nullptr) {
            // Zero means "not declared" (stage-1 convention); a declared
            // value must be finite.
            if (spec->value_type == cfdx::io::BCValueType::WALL_THERMAL &&
                spec->wall_temperature != 0.0) {
                if (!std::isfinite(spec->wall_temperature))
                    throw std::invalid_argument(
                        "production energy stage: non-finite wall_temperature "
                        "on patch '" + patch.name + "'");
                bc.type = ScalarBoundaryType::FIXED_VALUE;
                bc.value = spec->wall_temperature;
                fixed_value = true;
            } else if ((spec->type == cfdx::io::BCType::INLET ||
                        spec->value_type == cfdx::io::BCValueType::FIXED) &&
                       spec->temperature != 0.0) {
                if (!std::isfinite(spec->temperature))
                    throw std::invalid_argument(
                        "production energy stage: non-finite inlet "
                        "temperature on patch '" + patch.name + "'");
                bc.type = ScalarBoundaryType::FIXED_VALUE;
                bc.value = spec->temperature;
                fixed_value = true;
            }
        }
        if (!fixed_value && spec != nullptr) {
            const auto flux_it = spec->raw_params.find("heat_flux");
            if (flux_it != spec->raw_params.end()) {
                double heat_flux = 0.0;
                try {
                    heat_flux = std::stod(flux_it->second);
                } catch (const std::exception&) {
                    throw std::invalid_argument(
                        "production energy stage: heat_flux must be numeric "
                        "on patch '" + patch.name + "': '" + flux_it->second +
                        "'");
                }
                if (!std::isfinite(heat_flux))
                    throw std::invalid_argument(
                        "production energy stage: heat_flux must be finite on "
                        "patch '" + patch.name + "'");
                // Positive heat_flux enters the domain (stage-1 convention);
                // the FIXED_GRADIENT assembly then contributes +k*A*g to the
                // right-hand side with g = q''/k, an entering flux.
                bc.type = ScalarBoundaryType::FIXED_GRADIENT;
                bc.gradient = heat_flux / energy_case.conductivity;
            }
        }
        if (bc.type == ScalarBoundaryType::FIXED_VALUE)
            ++dirichlet_patches;
        bcs[patch.name] = bc;
    }
    if (dirichlet_patches == 0)
        throw std::invalid_argument(
            "production energy stage: no Dirichlet temperature boundary "
            "(wall_temperature or inlet temperature); an all-Neumann energy "
            "system has a constant null space and is rejected");
    return bcs;
}

struct ProductionEnergyStageResult {
    bool converged = false;
    std::size_t iterations = 0;
    double min_temperature = std::numeric_limits<double>::infinity();
    double max_temperature = -std::numeric_limits<double>::infinity();
    double total_source_power = 0.0;  // W
    double boundary_flux_out = 0.0;    // W, net, positive leaving the domain
    std::map<std::string, double> patch_flux_out;  // W per patch, positive out
    double independent_conservation_residual = 0.0;    // W
    double independent_conservation_normalized = 0.0;  // dimensionless
    std::size_t nonfinite_temperature = 0;
    std::vector<EnergyIteration> history;
    EnergySolverControls controls;
    double volumetric_power = 0.0;  // W/m^3
};

inline ProductionEnergyStageResult run_production_energy_stage(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::io::CaseSetup& setup,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
{
    validate_production_energy_model(setup);
    if (!production_energy_requested(setup))
        throw std::invalid_argument(
            "run_production_energy_stage: the energy stage only runs for "
            "energy_model 'conjugate'; running it under another declared "
            "model is a contract violation, not a silent fallback");
    if (mass_flux.size() != mesh.n_faces() || mass_flux.dimension() != 1)
        throw std::invalid_argument(
            "run_production_energy_stage: mass flux field does not match the "
            "mesh");
    for (std::size_t f = 0; f < mesh.n_faces(); ++f)
        if (!std::isfinite(mass_flux(f)))
            throw std::invalid_argument(
                "run_production_energy_stage: non-finite face mass flux");
    if (temperature.size() != mesh.n_cells() || temperature.dimension() != 1)
        throw std::invalid_argument(
            "run_production_energy_stage: temperature field does not match "
            "the mesh");
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        if (!std::isfinite(temperature(c)))
            throw std::invalid_argument(
                "run_production_energy_stage: non-finite initial temperature");

    const auto energy_case = build_production_energy_case(setup);
    const auto bcs = build_production_temperature_bcs(mesh, setup, energy_case);

    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "volumetric_power", "W/m3", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        source(c) = energy_case.volumetric_power;

    const auto previous_temperature = temperature;
    const auto solve = solve_energy(
        mesh, geometry, mass_flux, temperature, source,
        energy_case.controls, bcs);

    ProductionEnergyStageResult result;
    result.converged = solve.converged;
    result.iterations = solve.iterations;
    result.history = solve.history;
    result.controls = energy_case.controls;
    result.volumetric_power = energy_case.volumetric_power;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const double t = temperature(c);
        if (!std::isfinite(t)) {
            ++result.nonfinite_temperature;
            continue;
        }
        result.min_temperature = std::min(result.min_temperature, t);
        result.max_temperature = std::max(result.max_temperature, t);
    }
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        result.total_source_power += source(c) * geometry.cell_volumes[c];

    // Independent conservation audit: the energy flux is rebuilt from the
    // final temperature and the same authoritative mass flux, never from the
    // assembled matrix.
    const auto audit = reconstruct_energy_balance(
        mesh, geometry, mass_flux, temperature, previous_temperature,
        source, energy_case.controls, bcs);
    result.independent_conservation_residual = audit.residual;
    result.independent_conservation_normalized = audit.normalized_residual;
    result.boundary_flux_out = audit.boundary_flux;

    // Per-patch fluxes use the same independent reconstruction (temperature
    // units first, then scaled back by cp to physical watts).
    auto energy_flux = reconstruct_scalar_transport_flux(
        mesh, geometry, mass_flux, temperature,
        energy_case.conductivity / energy_case.specific_heat,
        bcs, true, ConvectionScheme::UPWIND, nullptr, nullptr, nullptr);
    for (std::size_t f = 0; f < energy_flux.size(); ++f)
        energy_flux(f) *= energy_case.specific_heat;
    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (own.neighbour(f) >= 0)
            continue;
        const std::size_t patch = geometry.face_patch[f];
        if (patch >= mesh.boundary().n_patches())
            continue;
        result.patch_flux_out[mesh.boundary().patch(patch).name] +=
            energy_flux(f);
    }
    return result;
}

} // namespace cfdx::physics
