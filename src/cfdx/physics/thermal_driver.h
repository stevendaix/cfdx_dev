// Autonomous thermal-only driver (issue #1150, A2+A3).
//
// Solves one mesh region of a ThermalCase without any Navier-Stokes state:
// steady or transient solid conduction, stagnant-fluid diffusion and
// prescribed-mass-flux advection, with typed thermal boundary conditions.
// The operator is exactly the production energy assembly (fiches 11/14);
// this driver adds boundary typing, per-cell k(T) outer Picard iteration
// and fail-closed contract enforcement.
// Method sheet: docs/audits/17-thermal-driver.md.
//
// A3: PRESCRIBED_MASS_FLUX fluid regions are advected with the authoritative
// face mass-flux field named by the region (mass_flux_field), supplied by the
// caller. The field must be discretely conservative and every inflow face
// must carry a temperature boundary (fixed_temperature or convection).
//
// Explicitly deferred (fail-closed with explicit diagnostics):
//   - conjugate interfaces              -> A4
//   - thermal links / lumped nodes      -> A5
//   - anisotropic conductivity tensors   (tensor diffusion operator)
#pragma once

#include "cfdx/physics/thermal_case.h"
#include "cfdx/physics/energy_solver.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

enum class ThermalBoundaryKind : std::uint8_t {
    FIXED_TEMPERATURE,  // K, Dirichlet
    FIXED_HEAT_FLUX,    // W/m^2, into the region
    CONVECTION,         // Robin film: h [W/m^2.K], T_ambient [K]
    POWER,              // W total on the patch, into the region
    ADIABATIC           // zero normal gradient
};

inline const char* to_string(ThermalBoundaryKind value) {
    switch (value) {
        case ThermalBoundaryKind::FIXED_TEMPERATURE: return "fixed_temperature";
        case ThermalBoundaryKind::FIXED_HEAT_FLUX: return "fixed_heat_flux";
        case ThermalBoundaryKind::CONVECTION: return "convection";
        case ThermalBoundaryKind::POWER: return "power";
        case ThermalBoundaryKind::ADIABATIC: return "adiabatic";
    }
    return "unknown";
}

// One typed thermal condition for one mesh patch. Validation is fail-closed
// against the region material: a Dirichlet value outside the declared
// material range is rejected before any assembly.
struct ThermalBoundaryCondition {
    ThermalBoundaryKind kind = ThermalBoundaryKind::ADIABATIC;
    // K for FIXED_TEMPERATURE, W/m^2 for FIXED_HEAT_FLUX, W for POWER.
    double value = 0.0;
    // CONVECTION only: film coefficient and ambient temperature.
    double heat_transfer_coefficient = 0.0;
    double ambient_temperature = 0.0;

    void validate(const ThermalMaterial& material) const {
        switch (kind) {
            case ThermalBoundaryKind::FIXED_TEMPERATURE:
                if (!std::isfinite(value))
                    throw std::invalid_argument(
                        "thermal boundary condition: non-finite fixed temperature");
                if (value < material.minimum_temperature ||
                    value > material.maximum_temperature)
                    throw std::invalid_argument(
                        "thermal boundary condition: fixed temperature is outside the material range");
                break;
            case ThermalBoundaryKind::FIXED_HEAT_FLUX:
                if (!std::isfinite(value))
                    throw std::invalid_argument(
                        "thermal boundary condition: non-finite heat flux");
                break;
            case ThermalBoundaryKind::CONVECTION:
                if (!std::isfinite(heat_transfer_coefficient) ||
                    heat_transfer_coefficient <= 0.0)
                    throw std::invalid_argument(
                        "thermal boundary condition: convection requires a positive heat transfer coefficient");
                if (!std::isfinite(ambient_temperature))
                    throw std::invalid_argument(
                        "thermal boundary condition: non-finite ambient temperature");
                break;
            case ThermalBoundaryKind::POWER:
                if (!std::isfinite(value))
                    throw std::invalid_argument(
                        "thermal boundary condition: non-finite patch power");
                break;
            case ThermalBoundaryKind::ADIABATIC:
                break;
            default:
                throw std::invalid_argument(
                    "thermal boundary condition: unsupported boundary kind");
        }
    }
};

using ThermalBoundaryConditions = std::map<std::string, ThermalBoundaryCondition>;

// Authoritative face mass-flux fields by name [kg/s]. The sign convention is
// the production one: a positive flux leaves the owner cell of the face.
using ThermalFaceFields =
    std::map<std::string, cfdx::core::Field<double, cfdx::core::Location::FACE>>;

struct ThermalDriverControls {
    // Outer Picard convergence on the temperature update.
    double outer_temperature_tolerance = 1e-12;
    std::size_t max_outer_iterations = 100;
    // Passed to the energy stage (relaxation of the linear predictor).
    double relaxation = 1.0;
    std::size_t max_inner_iterations = 100;
    // Common energy tolerance (linear, temperature-change and balance gates).
    double tolerance = 1e-12;

    void validate() const {
        if (!std::isfinite(outer_temperature_tolerance) ||
            outer_temperature_tolerance <= 0.0 ||
            max_outer_iterations == 0 ||
            !std::isfinite(relaxation) ||
            !(relaxation > 0.0 && relaxation <= 1.0) ||
            max_inner_iterations == 0 ||
            !std::isfinite(tolerance) || tolerance <= 0.0)
            throw std::invalid_argument("invalid thermal driver controls");
    }
};

struct ThermalDriverResult {
    bool converged = false;
    std::size_t outer_iterations = 0;
    // Magnitude of the last outer Picard temperature update (0 for the
    // single-pass exact path, whose assembled system is state-independent).
    double max_outer_temperature_change = 0.0;
    EnergySolveResult energy;
    IndependentEnergyBalance balance;
    // Per-cell conductivity used by the final pass and the balance audit.
    std::vector<double> cell_conductivity;
};

// Advance one region of a thermal-only case by one solve (steady) or one
// physical time step (transient). temperature is both the input state (T^n,
// validated against the material range) and the output state. Solid and
// stagnant-fluid regions solve with a zero mass flux; a PRESCRIBED_MASS_FLUX
// fluid region advects the face field named by region->mass_flux_field,
// resolved in face_fields (fail-closed: present, mesh-sized, finite and
// discretely conservative; every inflow face needs a temperature boundary).
// Fail-closed: any unsupported coupling, boundary mapping or value is
// rejected before assembly with an explicit diagnostic.
inline ThermalDriverResult solve_thermal_region(
    const ThermalCase& thermal_case,
    const std::string& region_name,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ThermalBoundaryConditions& boundary_conditions,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& source,
    const ThermalDriverControls& controls = {},
    const ThermalFaceFields* face_fields = nullptr)
{
    controls.validate();

    // --- Contract lookups (fail closed) -------------------------------
    const ThermalRegion* region = nullptr;
    for (const auto& candidate : thermal_case.regions)
        if (candidate.name == region_name) { region = &candidate; break; }
    if (!region)
        throw std::invalid_argument(
            "solve_thermal_region: unknown thermal region '" + region_name + "'");

    const ThermalMaterial* material = nullptr;
    for (const auto& candidate : thermal_case.materials)
        if (candidate.name == region->material) { material = &candidate; break; }
    if (!material)
        throw std::invalid_argument(
            "solve_thermal_region: region '" + region_name +
            "' references unknown material '" + region->material + "'");

    const bool prescribed_flux =
        region->kind == ThermalRegionKind::FLUID &&
        region->fluid_model == ThermalFluidModel::PRESCRIBED_MASS_FLUX;

    if (material->anisotropic_conductivity)
        throw std::invalid_argument(
            "solve_thermal_region: anisotropic conductivity of material '" +
            material->name + "' requires the tensor diffusion operator (deferred)");

    for (const auto& interface : thermal_case.interfaces)
        if (interface.region_a == region_name || interface.region_b == region_name)
            throw std::invalid_argument(
                "solve_thermal_region: conjugate interface '" + interface.name +
                "' is deferred to A4");

    for (const auto& link : thermal_case.links)
        for (const ThermalEndpoint* endpoint : {&link.source, &link.target})
            if (endpoint->kind == ThermalEndpointKind::REGION_PATCH &&
                endpoint->object == region_name)
                throw std::invalid_argument(
                    "solve_thermal_region: thermal link '" + link.name +
                    "' is deferred to A5");

    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (temperature.size() != nc || temperature.dimension() != 1 ||
        source.size() != nc || source.dimension() != 1 ||
        geometry.cell_volumes.size() != nc ||
        geometry.cell_centres.size() != nc ||
        geometry.face_area_vectors.size() != nf ||
        geometry.face_patch.size() != nf)
        throw std::invalid_argument(
            "solve_thermal_region: field/geometry shape mismatch");

    for (std::size_t i = 0; i < nc; ++i)
        if (!std::isfinite(temperature(i)) ||
            temperature(i) < material->minimum_temperature ||
            temperature(i) > material->maximum_temperature)
            throw std::invalid_argument(
                "solve_thermal_region: initial temperature is outside the material range of region '" +
                region_name + "'");

    if (thermal_case.time_mode == ThermalTimeMode::TRANSIENT &&
        !(thermal_case.time_step > 0.0))
        throw std::invalid_argument(
            "solve_thermal_region: transient region requires a positive time step");

    // --- Prescribed mass-flux resolution (A3, fail closed) ------------
    cfdx::core::Field<double, cfdx::core::Location::FACE> zero_mass_flux(
        nf, "mass_flux", "kg/s", 1);
    zero_mass_flux.fill(0.0);
    const cfdx::core::Field<double, cfdx::core::Location::FACE>* mass_flux =
        &zero_mass_flux;
    if (prescribed_flux) {
        if (!face_fields)
            throw std::invalid_argument(
                "solve_thermal_region: prescribed-mass-flux region '" +
                region_name + "' requires a face mass-flux field map");
        const auto it = face_fields->find(region->mass_flux_field);
        if (it == face_fields->end())
            throw std::invalid_argument(
                "solve_thermal_region: unknown face mass-flux field '" +
                region->mass_flux_field + "' for region '" + region_name + "'");
        const auto& field = it->second;
        if (field.size() != nf || field.dimension() != 1)
            throw std::invalid_argument(
                "solve_thermal_region: face mass-flux field '" +
                region->mass_flux_field + "' does not match the mesh of region '" +
                region_name + "'");
        for (std::size_t f = 0; f < nf; ++f)
            if (!std::isfinite(field(f)))
                throw std::invalid_argument(
                    "solve_thermal_region: face mass-flux field '" +
                    region->mass_flux_field + "' has a non-finite value");
        // Discrete conservation: with frozen density a net mass source has
        // no physical counterpart in this model, so an unbalanced field is
        // rejected instead of silently creating phantom energy.
        std::vector<double> net(nc, 0.0);
        std::vector<double> through(nc, 0.0);
        for (std::size_t f = 0; f < nf; ++f) {
            const std::size_t owner = mesh.ownership().owner(f);
            const double flux = field(f);
            if (mesh.ownership().neighbour(f) >= 0) {
                const std::size_t neighbour =
                    static_cast<std::size_t>(mesh.ownership().neighbour(f));
                net[owner] += flux;
                net[neighbour] -= flux;
                through[owner] += std::abs(flux);
                through[neighbour] += std::abs(flux);
            } else {
                net[owner] += flux;
                through[owner] += std::abs(flux);
            }
        }
        for (std::size_t i = 0; i < nc; ++i)
            if (std::abs(net[i]) > 1e-9 * std::max(1.0, through[i]))
                throw std::invalid_argument(
                    "solve_thermal_region: face mass-flux field '" +
                    region->mass_flux_field +
                    "' is not discretely conservative in region '" + region_name + "'");
        mass_flux = &field;
    }

    // --- Boundary mapping (fail closed) -------------------------------
    // Every non-EMPTY mesh patch must have exactly one typed condition and
    // every condition must reference an existing patch. EMPTY patches have
    // no area-integrated thermal meaning for this driver and are rejected.
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        const auto& patch = mesh.boundary().patch(p);
        if (patch.type == cfdx::core::PatchType::EMPTY)
            throw std::invalid_argument(
                "solve_thermal_region: EMPTY patch '" + patch.name +
                "' is not supported by the thermal driver");
        if (boundary_conditions.find(patch.name) == boundary_conditions.end())
            throw std::invalid_argument(
                "solve_thermal_region: missing thermal boundary condition for patch '" +
                patch.name + "'");
    }
    for (const auto& entry : boundary_conditions) {
        bool found = false;
        for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p)
            if (mesh.boundary().patch(p).name == entry.first) { found = true; break; }
        if (!found)
            throw std::invalid_argument(
                "solve_thermal_region: thermal boundary condition references unknown patch '" +
                entry.first + "'");
        entry.second.validate(*material);
    }

    // Integrated patch areas, required to spread a POWER patch into a flux.
    std::map<std::string, double> patch_area;
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        const auto& patch = mesh.boundary().patch(p);
        double area = 0.0;
        for (const std::size_t face : patch.face_ids) {
            if (face >= nf)
                throw std::out_of_range(
                    "solve_thermal_region: patch '" + patch.name +
                    "' references an invalid face");
            area += geometry.face_area_vectors[face].mag();
        }
        patch_area.emplace(patch.name, area);
    }

    // --- Frozen material state ----------------------------------------
    // rho and cp are evaluated once at the declared initial state of the
    // region. Only conductivity participates in the outer Picard loop.
    const double rho0 =
        material->properties.rho(region->initial_temperature, region->pressure);
    const double cp0 = material->properties.cp(region->initial_temperature);
    const double dt =
        thermal_case.time_mode == ThermalTimeMode::TRANSIENT
            ? thermal_case.time_step : 0.0;

    EnergySolverControls energy_controls;
    energy_controls.density = rho0;
    energy_controls.cp = cp0;
    energy_controls.conductivity =
        material->properties.k(region->initial_temperature);
    energy_controls.dt = dt;
    energy_controls.relaxation = controls.relaxation;
    energy_controls.max_iterations = controls.max_inner_iterations;
    energy_controls.tolerance = controls.tolerance;

    // With a CONSTANT conductivity law and no CONVECTION boundary the
    // assembled system does not depend on the temperature state: a single
    // energy pass is exact and the outer loop degenerates to one pass.
    bool single_pass_exact =
        material->properties.conductivity.model == ScalarPropertyModel::CONSTANT;
    for (const auto& entry : boundary_conditions)
        if (entry.second.kind == ThermalBoundaryKind::CONVECTION)
            single_pass_exact = false;

    const cfdx::core::Field<double, cfdx::core::Location::CELL> previous_time =
        temperature;  // T^n for transient; identity for steady.

    ThermalDriverResult result;
    std::vector<double> conductivity(nc, 0.0);
    ScalarBoundaryFaceConditions face_conditions(nf);

    const std::size_t max_outer =
        single_pass_exact ? static_cast<std::size_t>(1) : controls.max_outer_iterations;
    for (std::size_t outer = 1; outer <= max_outer; ++outer) {
        // Per-cell conductivity at the current Picard state. k(T) throws on
        // non-positive values; no silent clamping is performed.
        for (std::size_t i = 0; i < nc; ++i)
            conductivity[i] = material->properties.k(temperature(i));

        // Face-resolved boundary conditions, rebuilt at every outer
        // iteration because FIXED_HEAT_FLUX, POWER and CONVECTION all
        // depend on the owner conductivity or the owner temperature.
        for (std::size_t f = 0; f < nf; ++f) {
            if (mesh.ownership().neighbour(f) >= 0)
                continue;  // internal face: harmonic cell conductance
            const std::size_t owner = mesh.ownership().owner(f);
            if (owner >= nc)
                throw std::runtime_error(
                    "solve_thermal_region: boundary face without an owner");
            const std::size_t patch_id = geometry.face_patch[f];
            if (patch_id >= mesh.boundary().n_patches())
                throw std::invalid_argument(
                    "solve_thermal_region: boundary face without a patch");
            const auto& patch = mesh.boundary().patch(patch_id);
            const ThermalBoundaryCondition& bc = boundary_conditions.at(patch.name);
            // A3: an inflow face needs an external temperature to advect
            // in; a flux or adiabatic condition would silently recirculate
            // the owner state instead.
            if (prescribed_flux && (*mass_flux)(f) < 0.0 &&
                bc.kind != ThermalBoundaryKind::FIXED_TEMPERATURE &&
                bc.kind != ThermalBoundaryKind::CONVECTION)
                throw std::invalid_argument(
                    "solve_thermal_region: inflow boundary patch '" + patch.name +
                    "' of prescribed-mass-flux region '" + region_name +
                    "' requires a fixed_temperature or convection condition");
            const double k_owner = conductivity[owner];
            ScalarBoundaryFaceCondition resolved;
            switch (bc.kind) {
                case ThermalBoundaryKind::FIXED_TEMPERATURE:
                    resolved.type = ScalarBoundaryFaceCondition::Type::FIXED_VALUE;
                    resolved.value = bc.value;
                    break;
                case ThermalBoundaryKind::ADIABATIC:
                    resolved.type = ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT;
                    resolved.value = 0.0;
                    break;
                case ThermalBoundaryKind::FIXED_HEAT_FLUX: {
                    // q into the region => outward normal gradient q/k_owner,
                    // integrated by assembly as an outward flux -k*A*g.
                    resolved.type = ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT;
                    resolved.value = bc.value / k_owner;
                    break;
                }
                case ThermalBoundaryKind::POWER: {
                    const double area = patch_area.at(patch.name);
                    if (!(area > 0.0) || !std::isfinite(area))
                        throw std::invalid_argument(
                            "solve_thermal_region: patch '" + patch.name +
                            "' has a degenerate area for a POWER boundary");
                    resolved.type = ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT;
                    resolved.value = (bc.value / area) / k_owner;
                    break;
                }
                case ThermalBoundaryKind::CONVECTION: {
                    // Discrete Robin fixed point: equating the conductive
                    // half-cell flux k*(T_b - T_P)/d with the film flux
                    // h*(T_amb - T_b) gives the face value below. The same
                    // face value is used by assembly and by the balance
                    // audit, so the reconstructed film flux is consistent.
                    const double distance =
                        boundary_normal_distance(geometry, f, owner);
                    const double h = bc.heat_transfer_coefficient;
                    resolved.type = ScalarBoundaryFaceCondition::Type::FIXED_VALUE;
                    resolved.value =
                        (k_owner * temperature(owner) + h * distance * bc.ambient_temperature) /
                        (k_owner + h * distance);
                    break;
                }
                default:
                    throw std::invalid_argument(
                        "solve_thermal_region: unsupported thermal boundary kind for patch '" +
                        patch.name + "'");
            }
            face_conditions.conditions[f] = resolved;
            face_conditions.specified[f] = static_cast<unsigned char>(1);
        }

        const cfdx::core::Field<double, cfdx::core::Location::CELL> state_before =
            temperature;
        result.energy = solve_energy(
            mesh, geometry, *mass_flux, temperature, source, energy_controls,
            {}, nullptr, nullptr, &previous_time, &conductivity, &face_conditions);

        double max_change = 0.0;
        for (std::size_t i = 0; i < nc; ++i)
            max_change = std::max(
                max_change, std::abs(temperature(i) - state_before(i)));
        result.outer_iterations = outer;
        result.max_outer_temperature_change = max_change;

        if (single_pass_exact) {
            result.converged = result.energy.converged;
            break;
        }
        // A failed linear solve cannot improve by rebuilding the identical
        // system; preserve the diagnostics and stop (converged stays false).
        if (!result.energy.linear_solver_converged)
            break;
        if (result.energy.converged &&
            max_change <= controls.outer_temperature_tolerance) {
            result.converged = true;
            break;
        }
    }

    if (result.converged) {
        for (std::size_t i = 0; i < nc; ++i)
            if (!std::isfinite(temperature(i)) ||
                temperature(i) < material->minimum_temperature ||
                temperature(i) > material->maximum_temperature)
                throw std::runtime_error(
                    "solve_thermal_region: converged temperature is outside the material range of region '" +
                    region_name + "'");
    }

    // Independent post-solve audit: physical watts reconstructed from the
    // final state with the same cell conductivities and face conditions as
    // the final assembly pass.
    result.balance = reconstruct_energy_balance(
        mesh, geometry, *mass_flux, temperature, previous_time, source,
        energy_controls, {}, nullptr, nullptr, &conductivity, &face_conditions);
    result.cell_conductivity = conductivity;
    return result;
}

} // namespace cfdx::physics
