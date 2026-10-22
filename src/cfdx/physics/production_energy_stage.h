#pragma once

// T2 production energy stage: frozen-flow coupling of the incompressible
// solver's authoritative conservative face mass flux with the steady energy
// equation (capability #885 stage 2, issues #880 / #118).
//
// Contract:
//   * energy_model gate: {"", "none", "isothermal", "solid_conduction",
//     "conjugate", "buoyant", "compressible"} are accepted declarations;
//     anything else is rejected fail-closed instead of silently
//     degrading to an isothermal run. "compressible" and "solid_conduction"
//     are valid SCHEMA declarations (e.g. "compressible" is produced by the
//     SU2 conversion of Euler / Navier-Stokes compressible cases,
//     adapters/su2.py), and they never trigger this stage (only "conjugate"
//     and "buoyant" do) - but the production solver APPLICATION rejects them
//     with an
//     explicit routing error: this executable implements the incompressible
//     carrier, and running it on a case that declares another physics would
//     silently compute something else than the case asks for. The stage
//     itself only runs for "conjugate" (frozen-flow coupled model) and
//     "buoyant" (re-solved per outer iteration of the flow-energy buoyancy
//     loop, production_buoyancy_stage.h, issue #885 step 4; the stage
//     kernel, gates and conservation audits are identical under both
//     models): calling it under any other declared
//     model throws.
//   * Exactly one fluid material with finite positive density, dynamic
//     viscosity, conductivity and specific heat. The carrier flow solves
//     the same declared fluid, so an incomplete material is rejected whole
//     instead of letting solver defaults silently substitute for the
//     declared properties. Volumetric power (W/m^3) comes from the
//     material's extra_properties ("volumetric_power", default 0) and must
//     be a strict finite numeric literal: a numeric prefix such as
//     "500abc" is rejected whole, not silently truncated to 500.
//   * Temperature boundary conditions from the case: WALL_THERMAL with a
//     non-zero wall_temperature -> Dirichlet; INLET or FIXED with a non-zero
//     temperature -> Dirichlet; otherwise raw_params["heat_flux"]
//     (W/m^2, positive into the domain) plus the optional
//     radiative_heat_flux_into_domain patch map -> fixed gradient q''/k;
//     otherwise adiabatic. Both flux inputs are rejected on non-WALL
//     patches, and a patch declaring both a retained Dirichlet temperature
//     and either flux input is rejected as contradictory: exactly one thermal
//     boundary condition per patch. Zero means "not declared" (stage-1 /
//     T1 convention), so a Dirichlet at exactly 0 K is not representable
//     in this API; any declared temperature must be finite and strictly
//     positive (kelvin). At least one Dirichlet patch is required: an
//     all-Neumann energy system has a constant null space and is rejected.
//   * The stage consumes the flow result's authoritative conservative face
//     mass flux, never a flux reconstructed from the cell velocity, and the
//     finite-volume geometry must match the mesh (sizes, positive finite
//     volumes, finite centres and area vectors, valid patch indices).
//   * Conservation audits are independent of the assembled matrix: the
//     energy flux is rebuilt from the final temperature and the same mass
//     flux (reconstruct_energy_balance), and per-patch fluxes are rebuilt
//     the same way. These audits are BLOCKING acceptance criteria, not
//     diagnostics: result.converged requires the linear solve to converge
//     AND the independent global energy balance AND the independent local
//     (per-cell) balance to close within energy_balance_tolerance, the
//     local one in two normalizations - the gross physical power scale
//     (total boundary throughput + source) and each cell's OWN gross
//     power scale (its face throughputs + source + accumulation, floor
//     1 W), so a local defect cannot hide behind a large domain scale -
//     AND every temperature/flux/cell residual to be finite. A solution
//     whose independent audit fails cannot be declared converged.

#include "cfdx/core/field/field.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/finite_volume_transport.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
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
    // "conjugate" = frozen-flow coupled model (one solve after the
    // carrier). "buoyant" (issue #885 step 4) re-solves the same
    // frozen-flow equation per outer iteration of the flow-energy
    // buoyancy loop; the coupled fixed-point verdict belongs to that
    // loop.
    return setup.energy_model == "conjugate" ||
           setup.energy_model == "buoyant";
}

inline void validate_production_energy_model(const cfdx::io::CaseSetup& setup)
{
    const std::string& model = setup.energy_model;
    if (model.empty() || model == "none" || model == "isothermal" ||
        model == "solid_conduction" || model == "conjugate" ||
        model == "compressible" || model == "buoyant")
        return;
    throw std::invalid_argument(
        "production energy stage: unsupported energy_model '" + model +
        "' (accepted: none, isothermal, solid_conduction, conjugate, "
        "buoyant, compressible; anything else is rejected instead of "
        "silently running an isothermal solve)");
}

// Strict numeric parsing for case-declared parameters. std::stod alone
// accepts a numeric prefix without consuming the whole string ("500abc"
// silently becomes 500); here the entire text must parse as one finite
// number, and callers still own the physical sign/quotient constraints.
inline double parse_case_double(const std::string& text, const std::string& name)
{
    if (text.empty())
        throw std::invalid_argument(
            "production energy stage: " + name + " must be numeric: ''");
    std::size_t parsed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &parsed);
    } catch (const std::exception&) {
        throw std::invalid_argument(
            "production energy stage: " + name + " must be numeric: '" +
            text + "'");
    }
    if (parsed != text.size())
        throw std::invalid_argument(
            "production energy stage: " + name + " must be numeric: '" +
            text + "' (trailing characters)");
    if (!std::isfinite(value))
        throw std::invalid_argument(
            "production energy stage: " + name + " must be finite: '" +
            text + "'");
    return value;
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
    // The carrier flow is solved with the same declared material, so its
    // dynamic viscosity is part of the physical contract too: an incomplete
    // fluid material must be rejected here instead of letting the flow run
    // with silently substituted default properties (a velocity field and
    // mass flux computed with a different fluid than the one the energy
    // balance declares).
    if (!std::isfinite(material.dynamic_viscosity) ||
        !(material.dynamic_viscosity > 0.0))
        throw std::invalid_argument(
            "production energy stage: material dynamic_viscosity must be "
            "finite and positive (the carrier flow solves the same "
            "declared fluid)");

    double volumetric_power = 0.0;
    const auto power_it = material.extra_properties.find("volumetric_power");
    if (power_it != material.extra_properties.end())
        volumetric_power =
            parse_case_double(power_it->second, "volumetric_power");

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
    // Calibrated audit floor (diag run 39633631822): the local and
    // per-cell normalized energy balances bottom out at the
    // floating-point summation floor of the independent flux
    // reconstruction (~1.4e-12 measured at audit_power_scale ~8.4e3
    // W) while the global balance stays at ~1.2e-14 - a summation
    // artifact, not a physical defect. The gate is 1e-11: a margin
    // of more than 7x above the measured floor, still far below any
    // locally significant defect. All three forms of the gate
    // (global, local, per-cell normalized) share this control.
    energy_case.controls.energy_balance_tolerance = 1e-11;
    validate_energy_controls(energy_case.controls);
    return energy_case;
}

inline ScalarBoundaryConditions build_production_temperature_bcs(
    const cfdx::core::Mesh& mesh,
    const cfdx::io::CaseSetup& setup,
    const ProductionEnergyCase& energy_case,
    const std::map<std::string, double>& radiative_heat_flux_into_domain = {})
{
    if (!std::isfinite(energy_case.conductivity) ||
        !(energy_case.conductivity > 0.0))
        throw std::invalid_argument(
            "build_production_temperature_bcs: invalid conductivity");

    // Coupling seam for a radiation outer iteration: values are patch-
    // averaged net heat fluxes in W/m^2, positive into the fluid. The caller
    // owns radiosity/temperature iteration; this routine only maps the
    // already-computed flux into the conservative energy boundary operator.
    // Reject misspelled/non-wall patch names and invalid values fail-closed.
    for (const auto& entry : radiative_heat_flux_into_domain) {
        if (entry.first.empty() || !std::isfinite(entry.second))
            throw std::invalid_argument(
                "production energy stage: radiative wall heat flux must have "
                "a non-empty patch name and finite value");
    }

    ScalarBoundaryConditions bcs;
    std::set<std::string> applied_radiative_fluxes;
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
            // value must be finite and physically admissible. Temperatures
            // are absolute (kelvin): a declared value must be > 0. 0 K is
            // not representable in this API, exactly like "not declared".
            if (spec->value_type == cfdx::io::BCValueType::WALL_THERMAL &&
                spec->wall_temperature != 0.0) {
                if (!std::isfinite(spec->wall_temperature) ||
                    !(spec->wall_temperature > 0.0))
                    throw std::invalid_argument(
                        "production energy stage: wall_temperature on patch "
                        "'" + patch.name + "' must be finite and positive "
                        "(kelvin)");
                bc.type = ScalarBoundaryType::FIXED_VALUE;
                bc.value = spec->wall_temperature;
                fixed_value = true;
            } else if ((spec->type == cfdx::io::BCType::INLET ||
                        spec->value_type == cfdx::io::BCValueType::FIXED) &&
                       spec->temperature != 0.0) {
                if (!std::isfinite(spec->temperature) ||
                    !(spec->temperature > 0.0))
                    throw std::invalid_argument(
                        "production energy stage: inlet temperature on "
                        "patch '" + patch.name + "' must be finite and "
                        "positive (kelvin)");
                bc.type = ScalarBoundaryType::FIXED_VALUE;
                bc.value = spec->temperature;
                fixed_value = true;
            }
        }
        const auto radiative_it =
            radiative_heat_flux_into_domain.find(patch.name);
        const bool has_radiative_flux =
            radiative_it != radiative_heat_flux_into_domain.end();
        if (has_radiative_flux) {
            if (spec == nullptr || spec->type != cfdx::io::BCType::WALL)
                throw std::invalid_argument(
                    "production energy stage: radiative heat flux is only "
                    "supported on declared WALL patches, not patch '" +
                    patch.name + "'");
            if (fixed_value)
                throw std::invalid_argument(
                    "production energy stage: patch '" + patch.name +
                    "' declares a fixed temperature and a radiative heat "
                    "flux; a radiating wall must be flux-coupled, not "
                    "temperature-clamped");
            applied_radiative_fluxes.insert(patch.name);
        }

        if (spec != nullptr) {
            const auto flux_it = spec->raw_params.find("heat_flux");
            if (flux_it != spec->raw_params.end() || has_radiative_flux) {
                // heat_flux (W/m^2, positive into the domain) is a WALL
                // boundary condition: fiche 14 defines it as an imposed
                // wall flux. Declaring it on an inlet or outlet patch is a
                // case-contract violation, rejected fail-closed instead of
                // silently imposing a temperature gradient on a flow
                // boundary.
                if (spec->type != cfdx::io::BCType::WALL)
                    throw std::invalid_argument(
                        "production energy stage: heat_flux is only "
                        "supported on WALL patches, not on patch '" +
                        patch.name + "'");
                // A retained Dirichlet temperature (wall_temperature or
                // inlet temperature) contradicts a declared flux: both
                // cannot hold on the same patch. Reject the ambiguous
                // declaration instead of silently letting one win.
                if (fixed_value)
                    throw std::invalid_argument(
                        "production energy stage: patch '" + patch.name +
                        "' declares both a fixed temperature and "
                        "heat_flux; exactly one thermal boundary condition "
                        "per patch is required");
                const double imposed_heat_flux =
                    flux_it == spec->raw_params.end()
                        ? 0.0
                        : parse_case_double(
                              flux_it->second,
                              "heat_flux on patch '" + patch.name + "'");
                const double radiative_heat_flux =
                    has_radiative_flux ? radiative_it->second : 0.0;
                const double heat_flux =
                    imposed_heat_flux + radiative_heat_flux;
                if (!std::isfinite(heat_flux))
                    throw std::invalid_argument(
                        "production energy stage: combined wall heat flux "
                        "on patch '" + patch.name + "' must be finite");
                // The quotient can overflow even when both operands are
                // finite; the assembled gradient must itself be finite. Both
                // flux contributions use the same positive-into-domain sign.
                const double gradient = heat_flux / energy_case.conductivity;
                if (!std::isfinite(gradient))
                    throw std::invalid_argument(
                        "production energy stage: heat_flux / conductivity "
                        "on patch '" + patch.name + "' must be finite");
                // Positive heat_flux enters the domain (stage-1
                // convention); the FIXED_GRADIENT assembly then contributes
                // +k*A*g to the right-hand side with g = q''/k, an
                // entering flux.
                bc.type = ScalarBoundaryType::FIXED_GRADIENT;
                bc.gradient = gradient;
            }
        }
        if (bc.type == ScalarBoundaryType::FIXED_VALUE)
            ++dirichlet_patches;
        bcs[patch.name] = bc;
    }
    for (const auto& entry : radiative_heat_flux_into_domain)
        if (applied_radiative_fluxes.count(entry.first) == 0)
            throw std::invalid_argument(
                "production energy stage: radiative heat flux references "
                "unknown boundary patch '" + entry.first + "'");
    if (dirichlet_patches == 0)
        throw std::invalid_argument(
            "production energy stage: no Dirichlet temperature boundary "
            "(wall_temperature or inlet temperature); an all-Neumann energy "
            "system has a constant null space and is rejected");
    return bcs;
}

struct ProductionEnergyStageResult {
    // Blocking verdict, see run_production_energy_stage: linear convergence
    // AND independent global balance AND independent local balance AND all
    // results finite. A red independent audit can never yield converged.
    bool converged = false;
    bool linear_solver_converged = false;
    int linear_solver_status = -1;
    double linear_residual_relative = std::numeric_limits<double>::infinity();
    bool global_energy_balance_passed = false;
    bool local_energy_balance_passed = false;
    std::size_t iterations = 0;
    double min_temperature = std::numeric_limits<double>::infinity();
    double max_temperature = -std::numeric_limits<double>::infinity();
    double total_source_power = 0.0;  // W
    double boundary_flux_out = 0.0;   // W, net, positive leaving the domain
    std::map<std::string, double> patch_flux_out;  // W per patch, positive out
    double independent_conservation_residual = 0.0;    // W
    double independent_conservation_normalized = 0.0;  // legacy reconstruct
                                                      // scale (floor 1 W)
    // Gross physical power scale the blocking audits are normalized by:
    // total boundary throughput |q_f| plus the source power. Unlike a net
    // scale, it never collapses on cases whose boundary fluxes cancel
    // (steady cases without volumetric source).
    double audit_power_scale = 0.0;          // W
    double global_balance_normalized = 0.0;  // |residual| / audit_power_scale
    double max_cell_residual = 0.0;           // W, independent local audit
    std::size_t worst_cell = 0;
    double local_balance_normalized = 0.0;  // max_cell_residual / scale
    // The same per-cell residual normalized by that cell's OWN gross
    // power scale (its face throughputs + its source + its accumulation,
    // floor 1 W): a locally significant defect cannot hide behind a
    // large domain-wide scale. Blocking gate, see local_energy_balance.
    double max_cell_residual_normalized = 0.0;
    std::size_t worst_cell_normalized = 0;
    std::size_t nonfinite_temperature = 0;
    std::size_t nonfinite_faces = 0;
    std::size_t nonfinite_cells = 0;
    std::vector<EnergyIteration> history;
    EnergySolverControls controls;
    double volumetric_power = 0.0;  // W/m^3
};

inline ProductionEnergyStageResult run_production_energy_stage(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::io::CaseSetup& setup,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    const std::map<std::string, double>& radiative_heat_flux_into_domain = {})
{
    validate_production_energy_model(setup);
    if (!production_energy_requested(setup))
        throw std::invalid_argument(
            "run_production_energy_stage: the energy stage only runs for "
            "energy_model 'conjugate' (frozen-flow coupled model) or "
            "'buoyant' (per-iteration model of the flow-energy buoyancy "
            "outer loop); running it under another declared model is a "
            "contract violation, not a silent fallback");
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
        if (!std::isfinite(temperature(c)) || !(temperature(c) > 0.0))
            throw std::invalid_argument(
                "run_production_energy_stage: initial temperature must be "
                "finite and positive (kelvin)");
    // Entry-point geometry contract: the kernel is a public header API and
    // can be called with a geometry built for another mesh; the assembly
    // would then index foreign arrays. Validate the contract here instead
    // of trusting the caller.
    if (geometry.cell_volumes.size() != mesh.n_cells() ||
        geometry.cell_centres.size() != mesh.n_cells() ||
        geometry.face_centres.size() != mesh.n_faces() ||
        geometry.face_area_vectors.size() != mesh.n_faces() ||
        geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument(
            "run_production_energy_stage: finite-volume geometry does not "
            "match the mesh");
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        if (!(geometry.cell_volumes[c] > 0.0) ||
            !std::isfinite(geometry.cell_volumes[c]))
            throw std::invalid_argument(
                "run_production_energy_stage: cell volumes must be finite "
                "and positive");
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto& centre = geometry.cell_centres[c];
        if (!std::isfinite(centre.x) || !std::isfinite(centre.y) ||
            !std::isfinite(centre.z))
            throw std::invalid_argument(
                "run_production_energy_stage: cell centres must be finite");
    }
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto& centre = geometry.face_centres[f];
        const auto& area = geometry.face_area_vectors[f];
        if (!std::isfinite(centre.x) || !std::isfinite(centre.y) ||
            !std::isfinite(centre.z) ||
            !std::isfinite(area.x) || !std::isfinite(area.y) ||
            !std::isfinite(area.z))
            throw std::invalid_argument(
                "run_production_energy_stage: face centres and area "
                "vectors must be finite");
    }
    // A boundary face must map to a declared patch: a geometry built for
    // another mesh (or a degenerate one) can carry patch indices that are
    // silently skipped by the per-patch flux audit, hiding flux that then
    // misses from the balance. Reject instead of skipping.
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0)
            continue;
        if (geometry.face_patch[f] >= mesh.boundary().n_patches())
            throw std::invalid_argument(
                "run_production_energy_stage: boundary face maps to an "
                "invalid patch index (geometry does not match the mesh)");
    }

    const auto energy_case = build_production_energy_case(setup);
    const auto bcs = build_production_temperature_bcs(
        mesh, setup, energy_case, radiative_heat_flux_into_domain);

    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "volumetric_power", "W/m3", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        source(c) = energy_case.volumetric_power;

    const auto previous_temperature = temperature;
    const auto solve = solve_energy(
        mesh, geometry, mass_flux, temperature, source,
        energy_case.controls, bcs);

    ProductionEnergyStageResult result;
    result.linear_solver_converged = solve.linear_solver_converged;
    result.linear_solver_status = solve.linear_solver_status;
    result.linear_residual_relative = solve.linear_residual_relative;
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
    // final temperature and the same authoritative mass flux, never from
    // the assembled matrix.
    const auto audit = reconstruct_energy_balance(
        mesh, geometry, mass_flux, temperature, previous_temperature,
        source, energy_case.controls, bcs);
    result.independent_conservation_residual = audit.residual;
    result.independent_conservation_normalized = audit.normalized_residual;
    result.boundary_flux_out = audit.boundary_flux;
    result.max_cell_residual = audit.max_cell_residual;
    result.worst_cell = audit.worst_cell;
    result.max_cell_residual_normalized = audit.max_cell_residual_normalized;
    result.worst_cell_normalized = audit.worst_cell_normalized;
    result.nonfinite_faces = audit.nonfinite_faces;
    result.nonfinite_cells = audit.nonfinite_cells;

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

    // Blocking acceptance criteria (fiche 14): the independent audit is a
    // verdict, not a diagnostic. Both the global and the local balance are
    // normalized by the gross physical power scale (total boundary
    // throughput + source), so a case whose boundary fluxes cancel (no
    // volumetric source) still has a meaningful, non-collapsing scale. The
    // local criterion is strictly stronger than the global one: it cannot
    // be satisfied by compensating positive and negative cell residuals.
    double audit_power_scale = std::abs(result.total_source_power);
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (own.neighbour(f) >= 0)
            continue;
        audit_power_scale += std::abs(energy_flux(f));
    }
    audit_power_scale = std::max(audit_power_scale, 1.0);
    const double balance_tolerance = energy_case.controls.energy_balance_tolerance;
    result.audit_power_scale = audit_power_scale;
    result.global_balance_normalized =
        std::abs(result.independent_conservation_residual) / audit_power_scale;
    result.local_balance_normalized = audit.max_cell_residual / audit_power_scale;
    result.global_energy_balance_passed =
        std::isfinite(result.independent_conservation_residual) &&
        std::isfinite(audit_power_scale) &&
        result.global_balance_normalized <= balance_tolerance;
    result.local_energy_balance_passed =
        std::isfinite(audit.max_cell_residual) &&
        result.local_balance_normalized <= balance_tolerance &&
        // Locally normalized per-cell balance: each cell's residual is
        // measured against its own gross power scale, so a locally
        // significant defect cannot hide behind a large domain-wide
        // scale. Both forms must pass.
        std::isfinite(audit.max_cell_residual_normalized) &&
        audit.max_cell_residual_normalized <= balance_tolerance;

    result.converged =
        result.linear_solver_converged &&
        result.global_energy_balance_passed &&
        result.local_energy_balance_passed &&
        result.nonfinite_temperature == 0 &&
        result.nonfinite_faces == 0 &&
        result.nonfinite_cells == 0 &&
        std::isfinite(result.min_temperature) &&
        std::isfinite(result.max_temperature);
    return result;
}

} // namespace cfdx::physics
