// Autonomous thermal-only driver (issue #1150, A2+A3+A4).
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
// A4: conjugate and contact-resistance interfaces between two meshed
// regions of the case. Each paired interface face is coupled with the
// same discrete Robin fixed point as CONVECTION: the supplied film
// coefficient is the inverse of the remaining series resistance
// (neighbour half-cell plus contact resistance) and the ambient state is
// the paired neighbour cell temperature. solve_thermal_regions_coupled
// iterates the injected ambient states to their fixed point, so the
// converged solution is the exact discrete series-coupled system and the
// interface fluxes are equal and opposite.
//
// Explicitly deferred (fail-closed with explicit diagnostics):
//   - thermal links / lumped nodes      -> A5
//   - anisotropic conductivity tensors   (tensor diffusion operator)
#pragma once

#include "cfdx/physics/thermal_case.h"
#include "cfdx/physics/energy_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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

// Per-face conjugate interface condition (A4). The driver resolves it with
// the same discrete Robin fixed point as CONVECTION, which already accounts
// for the owner half-cell conductance k/d; the supplied coefficient must
// therefore be the inverse of the REMAINING series resistance (neighbour
// half-cell plus contact resistance) and the ambient state is the paired
// neighbour cell temperature. Keyed by patch name, one entry per patch face
// in patch face order.
struct ThermalInterfaceFaceCondition {
    double heat_transfer_coefficient = 0.0;  // W/(m^2.K), strictly positive
    double ambient_temperature = 0.0;         // K, finite
};

using ThermalInterfaceFaceConditions =
    std::map<std::string, std::vector<ThermalInterfaceFaceCondition>>;

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
// A conjugate or contact-resistance interface patch of this region is solved
// with the per-face conditions supplied in interface_conditions (A4); every
// interface patch of the region must be covered, fail-closed.
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
    const ThermalFaceFields* face_fields = nullptr,
    const ThermalInterfaceFaceConditions* interface_conditions = nullptr)
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

    for (const auto& interface : thermal_case.interfaces) {
        if (interface.region_a != region_name && interface.region_b != region_name)
            continue;
        const std::string& patch_name =
            interface.region_a == region_name ? interface.patch_a : interface.patch_b;
        if (!interface_conditions ||
            interface_conditions->find(patch_name) == interface_conditions->end())
            throw std::invalid_argument(
                "solve_thermal_region: conjugate interface '" + interface.name +
                "' couples patch '" + patch_name + "' of region '" + region_name +
                "'; supply interface conditions or use the coupled solver");
    }

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
        const bool interface_patch =
            interface_conditions != nullptr &&
            interface_conditions->find(patch.name) != interface_conditions->end();
        if (interface_patch &&
            boundary_conditions.find(patch.name) != boundary_conditions.end())
            throw std::invalid_argument(
                "solve_thermal_region: patch '" + patch.name + "' of region '" +
                region_name +
                "' has both an external boundary condition and an interface condition");
        if (!interface_patch &&
            boundary_conditions.find(patch.name) == boundary_conditions.end())
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
    // The linear stage stops on a RELATIVE residual scaled by ||rhs||,
    // but the independent audit closes the balance in absolute watts
    // against a gross scale floored at 1 W. On a balanced steady state
    // the signed boundary fluxes cancel (in = out), so that floor becomes
    // the audit scale and a relative criterion alone can leave an absolute
    // cell residual of ||rhs||*tolerance that the audit then rejects. Cap
    // the linear tolerance so the assembled system is solved to an
    // absolute residual the audit can certify; the temperature-change and
    // balance gates keep the configured common tolerance.
    energy_controls.linear_tolerance = std::min(controls.tolerance, 1e-14);

    // With a CONSTANT conductivity law and no CONVECTION boundary the
    // assembled system does not depend on the temperature state: a single
    // energy pass is exact and the outer loop degenerates to one pass.
    bool single_pass_exact =
        material->properties.conductivity.model == ScalarPropertyModel::CONSTANT;
    for (const auto& entry : boundary_conditions)
        if (entry.second.kind == ThermalBoundaryKind::CONVECTION)
            single_pass_exact = false;
    // A4: an interface face is resolved through the state-dependent Robin
    // fixed point, exactly like CONVECTION.
    if (interface_conditions && !interface_conditions->empty())
        single_pass_exact = false;

    // Interface conditions validated once and indexed by face (A4).
    std::map<std::size_t, ThermalInterfaceFaceCondition> interface_by_face;
    if (interface_conditions) {
        for (const auto& entry : *interface_conditions) {
            bool patch_found = false;
            for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
                const auto& patch = mesh.boundary().patch(p);
                if (patch.name != entry.first)
                    continue;
                patch_found = true;
                if (entry.second.size() != patch.face_ids.size())
                    throw std::invalid_argument(
                        "solve_thermal_region: interface condition vector for patch '" +
                        patch.name + "' of region '" + region_name +
                        "' does not match the patch");
                for (std::size_t k = 0; k < patch.face_ids.size(); ++k) {
                    const ThermalInterfaceFaceCondition& ic = entry.second[k];
                    if (!std::isfinite(ic.heat_transfer_coefficient) ||
                        !(ic.heat_transfer_coefficient > 0.0) ||
                        !std::isfinite(ic.ambient_temperature))
                        throw std::invalid_argument(
                            "solve_thermal_region: interface condition of patch '" +
                            patch.name + "' requires a positive finite film "
                            "coefficient and a finite ambient temperature");
                    const std::size_t f = patch.face_ids[k];
                    if (!interface_by_face.emplace(f, ic).second)
                        throw std::invalid_argument(
                            "solve_thermal_region: face " + std::to_string(f) +
                            " carries more than one interface condition");
                }
                break;
            }
            if (!patch_found)
                throw std::invalid_argument(
                    "solve_thermal_region: interface condition references unknown patch '" +
                    entry.first + "' of region '" + region_name + "'");
        }
    }

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
            const double k_owner = conductivity[owner];
            ScalarBoundaryFaceCondition resolved;
            const auto iface_it = interface_by_face.find(f);
            if (iface_it != interface_by_face.end()) {
                // A4: paired interface face. The neighbour half-cell plus
                // contact resistance act as the film; the same discrete
                // Robin fixed point as CONVECTION resolves the face value,
                // so assembly and the balance audit stay consistent. A
                // prescribed mass flux crossing the interface is rejected:
                // the advected state of a paired face is not defined yet.
                if (prescribed_flux && (*mass_flux)(f) != 0.0)
                    throw std::invalid_argument(
                        "solve_thermal_region: interface face of patch '" +
                        patch.name + "' of region '" + region_name +
                        "' carries a prescribed mass flux, which is not supported");
                const ThermalInterfaceFaceCondition& ic = iface_it->second;
                const double distance =
                    boundary_normal_distance(geometry, f, owner);
                resolved.type = ScalarBoundaryFaceCondition::Type::FIXED_VALUE;
                resolved.value =
                    (k_owner * temperature(owner) +
                     ic.heat_transfer_coefficient * distance * ic.ambient_temperature) /
                    (k_owner + ic.heat_transfer_coefficient * distance);
            } else {
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

// ---------------------------------------------------------------------------
// A4: coupled solve of every meshed region of a case across its conjugate
// and contact-resistance interfaces.
// ---------------------------------------------------------------------------
struct ThermalCoupledControls {
    // Convergence of the injected interface ambient temperatures [K].
    double interface_temperature_tolerance = 1e-10;
    std::size_t max_interface_iterations = 200;

    void validate() const {
        if (!std::isfinite(interface_temperature_tolerance) ||
            !(interface_temperature_tolerance > 0.0) ||
            max_interface_iterations == 0)
            throw std::invalid_argument(
                "invalid thermal coupled controls");
    }
};

// Per-region solve inputs. The temperature field is both the initial state
// and the output state of the coupled solve.
struct ThermalRegionInputs {
    const cfdx::core::Mesh* mesh = nullptr;
    const FvGeometry* geometry = nullptr;
    ThermalBoundaryConditions boundary_conditions = {};
    cfdx::core::Field<double, cfdx::core::Location::CELL>* temperature = nullptr;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* source = nullptr;
    const ThermalFaceFields* face_fields = nullptr;
};

struct ThermalCoupledResult {
    bool converged = false;
    std::size_t interface_iterations = 0;
    // Largest ambient-temperature change of the final sweep [K].
    double interface_max_change = 0.0;
    // Largest violation of the equal/opposite interface flux audit [W],
    // reconstructed from the final states with the conditions the final
    // passes actually used.
    double interface_flux_imbalance = 0.0;
    std::map<std::string, ThermalDriverResult> regions;
};

inline ThermalCoupledResult solve_thermal_regions_coupled(
    const ThermalCase& thermal_case,
    const std::map<std::string, ThermalRegionInputs>& region_inputs,
    const ThermalDriverControls& controls = {},
    const ThermalCoupledControls& coupled_controls = {})
{
    coupled_controls.validate();

    // --- Region inputs (fail closed) ----------------------------------
    if (region_inputs.empty())
        throw std::invalid_argument(
            "solve_thermal_regions_coupled: no region inputs");
    for (const auto& entry : region_inputs) {
        bool known = false;
        for (const auto& region : thermal_case.regions)
            if (region.name == entry.first) { known = true; break; }
        if (!known)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: unknown region input '" +
                entry.first + "'");
    }
    for (const auto& region : thermal_case.regions) {
        const auto it = region_inputs.find(region.name);
        if (it == region_inputs.end())
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: missing region input '" +
                region.name + "'");
        const ThermalRegionInputs& in = it->second;
        if (!in.mesh || !in.geometry || !in.temperature || !in.source)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: region '" + region.name +
                "' requires a mesh, a geometry, a temperature field and a source field");
    }

    // --- Interface pairing (fail closed) ------------------------------
    struct InterfacePair {
        std::size_t face_a = 0, face_b = 0;
        std::size_t cell_a = 0, cell_b = 0;
        double area = 0.0;
        double distance_a = 0.0, distance_b = 0.0;
    };
    struct InterfaceData {
        const ThermalInterface* interface = nullptr;
        std::vector<InterfacePair> pairs;
    };
    std::vector<InterfaceData> interfaces;
    for (const auto& interface : thermal_case.interfaces) {
        InterfaceData data;
        data.interface = &interface;
        const ThermalRegionInputs& in_a = region_inputs.at(interface.region_a);
        const ThermalRegionInputs& in_b = region_inputs.at(interface.region_b);

        const cfdx::core::Mesh& mesh_a = *in_a.mesh;
        const cfdx::core::Mesh& mesh_b = *in_b.mesh;
        const FvGeometry& geom_a = *in_a.geometry;
        const FvGeometry& geom_b = *in_b.geometry;

        const std::size_t no_patch = static_cast<std::size_t>(-1);
        std::size_t index_a = no_patch;
        std::size_t index_b = no_patch;
        for (std::size_t p = 0; p < mesh_a.boundary().n_patches(); ++p)
            if (mesh_a.boundary().patch(p).name == interface.patch_a) {
                index_a = p;
                break;
            }
        for (std::size_t p = 0; p < mesh_b.boundary().n_patches(); ++p)
            if (mesh_b.boundary().patch(p).name == interface.patch_b) {
                index_b = p;
                break;
            }
        if (index_a == no_patch)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: interface '" + interface.name +
                "' references unknown patch '" + interface.patch_a +
                "' of region '" + interface.region_a + "'");
        if (index_b == no_patch)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: interface '" + interface.name +
                "' references unknown patch '" + interface.patch_b +
                "' of region '" + interface.region_b + "'");
        const auto& patch_a = mesh_a.boundary().patch(index_a);
        const auto& patch_b = mesh_b.boundary().patch(index_b);
        if (patch_a.type == cfdx::core::PatchType::EMPTY ||
            patch_b.type == cfdx::core::PatchType::EMPTY)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: interface '" + interface.name +
                "' references an EMPTY patch");
        if (in_a.boundary_conditions.count(interface.patch_a) != 0 ||
            in_b.boundary_conditions.count(interface.patch_b) != 0)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: interface patch of '" +
                interface.name + "' also has an external boundary condition");
        if (patch_a.face_ids.size() != patch_b.face_ids.size())
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: interface '" + interface.name +
                "' patches are not conformal (different face counts)");

        std::vector<bool> taken(patch_b.face_ids.size(), false);
        for (const std::size_t fa : patch_a.face_ids) {
            std::size_t best_k = patch_b.face_ids.size();
            double best_distance = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < patch_b.face_ids.size(); ++k) {
                const double distance =
                    (geom_a.face_centres[fa] -
                     geom_b.face_centres[patch_b.face_ids[k]]).mag();
                if (distance < best_distance) {
                    best_distance = distance;
                    best_k = k;
                }
            }
            const std::size_t fb = patch_b.face_ids[best_k];
            // Mutual nearest neighbour.
            std::size_t reverse_k = patch_a.face_ids.size();
            double reverse_distance = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < patch_a.face_ids.size(); ++k) {
                const double distance =
                    (geom_b.face_centres[fb] -
                     geom_a.face_centres[patch_a.face_ids[k]]).mag();
                if (distance < reverse_distance) {
                    reverse_distance = distance;
                    reverse_k = k;
                }
            }
            const double area_a = geom_a.face_area_vectors[fa].mag();
            const double area_b = geom_b.face_area_vectors[fb].mag();
            if (taken[best_k] ||
                patch_a.face_ids[reverse_k] != fa ||
                !(best_distance <= 1e-9 * std::sqrt(std::max(area_a, area_b))))
                throw std::invalid_argument(
                    "solve_thermal_regions_coupled: interface '" + interface.name +
                    "' patches are not conformal (unpaired faces)");
            if (std::abs(area_a - area_b) > 1e-9 * std::max(area_a, area_b))
                throw std::invalid_argument(
                    "solve_thermal_regions_coupled: interface '" + interface.name +
                    "' paired faces have mismatched areas");
            if (geom_a.face_area_vectors[fa].dot(geom_b.face_area_vectors[fb]) >= 0.0)
                throw std::invalid_argument(
                    "solve_thermal_regions_coupled: interface '" + interface.name +
                    "' paired faces do not face each other");
            taken[best_k] = true;

            InterfacePair pair;
            pair.face_a = fa;
            pair.face_b = fb;
            pair.area = area_a;
            pair.cell_a = mesh_a.ownership().owner(fa);
            pair.cell_b = mesh_b.ownership().owner(fb);
            if (pair.cell_a >= mesh_a.n_cells() ||
                pair.cell_b >= mesh_b.n_cells())
                throw std::invalid_argument(
                    "solve_thermal_regions_coupled: interface '" + interface.name +
                    "' paired face without an owner cell");
            pair.distance_a = boundary_normal_distance(geom_a, fa, pair.cell_a);
            pair.distance_b = boundary_normal_distance(geom_b, fb, pair.cell_b);
            data.pairs.push_back(pair);
        }
        interfaces.push_back(data);
    }

    // Materials for the interface film coefficients.
    std::map<std::string, const ThermalMaterial*> materials;
    for (const auto& region : thermal_case.regions) {
        const ThermalMaterial* material = nullptr;
        for (const auto& candidate : thermal_case.materials)
            if (candidate.name == region.material) { material = &candidate; break; }
        if (!material)
            throw std::invalid_argument(
                "solve_thermal_regions_coupled: region '" + region.name +
                "' references unknown material '" + region.material + "'");
        materials.emplace(region.name, material);
    }

    // --- Interface coupling iteration ---------------------------------
    // For transient cases every sweep restarts from the same T^n, so one
    // coupled sweep always advances the physical state by exactly one
    // time step regardless of the number of coupling sweeps.
    std::map<std::string, cfdx::core::Field<double, cfdx::core::Location::CELL>>
        start_state;
    for (const auto& region : thermal_case.regions)
        start_state.emplace(region.name, *region_inputs.at(region.name).temperature);

    // Ambients injected during the current sweep, for convergence and for
    // the final equal/opposite flux audit.
    std::vector<std::vector<double>> ambient_a(interfaces.size());
    std::vector<std::vector<double>> ambient_b(interfaces.size());
    std::vector<std::vector<double>> h_a(interfaces.size());
    std::vector<std::vector<double>> h_b(interfaces.size());
    for (std::size_t ii = 0; ii < interfaces.size(); ++ii) {
        const std::size_t np = interfaces[ii].pairs.size();
        ambient_a[ii].assign(np, std::numeric_limits<double>::quiet_NaN());
        ambient_b[ii].assign(np, std::numeric_limits<double>::quiet_NaN());
        h_a[ii].assign(np, 0.0);
        h_b[ii].assign(np, 0.0);
    }

    ThermalCoupledResult result;
    const double tolerance = coupled_controls.interface_temperature_tolerance;
    for (std::size_t sweep = 1; sweep <= coupled_controls.max_interface_iterations;
         ++sweep) {
        result.interface_iterations = sweep;
        for (const auto& region : thermal_case.regions) {
            const ThermalRegionInputs& in = region_inputs.at(region.name);
            *in.temperature = start_state.at(region.name);

            // Interface conditions of this region from the freshest
            // neighbour states (Gauss-Seidel within the sweep).
            ThermalInterfaceFaceConditions conditions;
            for (std::size_t ii = 0; ii < interfaces.size(); ++ii) {
                const ThermalInterface& iface = *interfaces[ii].interface;
                const bool side_a = iface.region_a == region.name;
                if (!side_a && iface.region_b != region.name)
                    continue;
                const std::string& other_name =
                    side_a ? iface.region_b : iface.region_a;
                const ThermalRegionInputs& other = region_inputs.at(other_name);
                const ThermalMaterial* other_material = materials.at(other_name);
                std::vector<ThermalInterfaceFaceCondition> per_face;
                for (std::size_t pi = 0; pi < interfaces[ii].pairs.size(); ++pi) {
                    const InterfacePair& pair = interfaces[ii].pairs[pi];
                    const std::size_t my_cell = side_a ? pair.cell_a : pair.cell_b;
                    const std::size_t other_cell =
                        side_a ? pair.cell_b : pair.cell_a;
                    const double other_distance =
                        side_a ? pair.distance_b : pair.distance_a;
                    (void)my_cell;
                    // Neighbour conductivity: final Picard conductivity of
                    // the neighbour's latest solve when available, else the
                    // material law at the current neighbour state.
                    double other_k = std::numeric_limits<double>::quiet_NaN();
                    const auto solved = result.regions.find(other_name);
                    if (solved != result.regions.end() &&
                        solved->second.cell_conductivity.size() == other.mesh->n_cells())
                        other_k = solved->second.cell_conductivity[other_cell];
                    else
                        other_k = other_material->properties.k(
                            (*other.temperature)(other_cell));
                    const double remaining =
                        other_distance / other_k + iface.contact_resistance;
                    ThermalInterfaceFaceCondition ic;
                    ic.heat_transfer_coefficient = 1.0 / remaining;
                    ic.ambient_temperature = (*other.temperature)(other_cell);
                    per_face.push_back(ic);
                    if (side_a) {
                        ambient_a[ii][pi] = ic.ambient_temperature;
                        h_a[ii][pi] = ic.heat_transfer_coefficient;
                    } else {
                        ambient_b[ii][pi] = ic.ambient_temperature;
                        h_b[ii][pi] = ic.heat_transfer_coefficient;
                    }
                }
                conditions.emplace(side_a ? iface.patch_a : iface.patch_b,
                                   std::move(per_face));
            }

            result.regions[region.name] = solve_thermal_region(
                thermal_case, region.name, *in.mesh, *in.geometry,
                in.boundary_conditions, *in.temperature, *in.source, controls,
                in.face_fields,
                conditions.empty() ? nullptr : &conditions);
        }

        // Convergence: the ambients the converged sweeps would inject next
        // must match the ambients actually injected by this sweep.
        double max_change = 0.0;
        for (std::size_t ii = 0; ii < interfaces.size(); ++ii) {
            const ThermalInterface& iface = *interfaces[ii].interface;
            const auto& in_a = region_inputs.at(iface.region_a);
            const auto& in_b = region_inputs.at(iface.region_b);
            for (std::size_t pi = 0; pi < interfaces[ii].pairs.size(); ++pi) {
                const InterfacePair& pair = interfaces[ii].pairs[pi];
                const double next_a = (*in_b.temperature)(pair.cell_b);
                const double next_b = (*in_a.temperature)(pair.cell_a);
                if (std::isfinite(ambient_a[ii][pi]))
                    max_change = std::max(max_change,
                        std::abs(next_a - ambient_a[ii][pi]));
                if (std::isfinite(ambient_b[ii][pi]))
                    max_change = std::max(max_change,
                        std::abs(next_b - ambient_b[ii][pi]));
            }
        }
        result.interface_max_change = max_change;
        if (max_change <= tolerance) {
            result.converged = true;
            break;
        }
    }

    // --- Equal/opposite interface flux audit ---------------------------
    // Reconstruct the boundary flux each final pass embedded at every
    // paired face: the Robin fixed point resolved with the final Picard
    // conductivity and the actually used film coefficient and ambient.
    double imbalance = 0.0;
    for (std::size_t ii = 0; ii < interfaces.size(); ++ii) {
        const ThermalInterface& iface = *interfaces[ii].interface;
        const auto& in_a = region_inputs.at(iface.region_a);
        const auto& in_b = region_inputs.at(iface.region_b);
        const auto& res_a = result.regions.at(iface.region_a);
        const auto& res_b = result.regions.at(iface.region_b);
        for (std::size_t pi = 0; pi < interfaces[ii].pairs.size(); ++pi) {
            const InterfacePair& pair = interfaces[ii].pairs[pi];
            const double k_a = res_a.cell_conductivity[pair.cell_a];
            const double k_b = res_b.cell_conductivity[pair.cell_b];
            const double t_a = (*in_a.temperature)(pair.cell_a);
            const double t_b = (*in_b.temperature)(pair.cell_b);
            const double face_a =
                (k_a * t_a + h_a[ii][pi] * pair.distance_a * ambient_a[ii][pi]) /
                (k_a + h_a[ii][pi] * pair.distance_a);
            const double face_b =
                (k_b * t_b + h_b[ii][pi] * pair.distance_b * ambient_b[ii][pi]) /
                (k_b + h_b[ii][pi] * pair.distance_b);
            const double q_out_a =
                k_a * pair.area * (face_a - t_a) / pair.distance_a;
            const double q_out_b =
                k_b * pair.area * (face_b - t_b) / pair.distance_b;
            imbalance = std::max(imbalance, std::abs(q_out_a + q_out_b));
        }
    }
    result.interface_flux_imbalance = imbalance;

    return result;
}

} // namespace cfdx::physics

