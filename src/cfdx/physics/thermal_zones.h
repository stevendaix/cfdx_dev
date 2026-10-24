// Production solid-conduction thermal path: multi-material cell zones,
// volumetric power sources and contact resistances at material interfaces.
//
// Stage 1 of the multiphysics solver-level orchestration capability
// (issue #885). Qualification case T1 (issue #880): steady conduction
// through composite solids with per-zone conductivity, per-zone volumetric
// power (W/m^3) and per-interface contact resistance (m^2.K/W).
//
// Design rules (fail-closed):
//   * A declared zone must have a finite positive conductivity.
//   * Volumetric power must be finite (zero is allowed).
//   * Contact resistance must be finite, non-negative, and applies at
//     internal faces separating two *different* declared zones. When both
//     sides declare a contact resistance the two values must be identical.
//   * Every cell must map onto a declared zone id.
//   * A conduction system with no fixed-temperature boundary anywhere is
//     singular (its constant null space leaves the temperature level
//     unpinned) and is rejected at assembly time.
//
// Discretization: standard two-point finite volumes. Internal faces use the
// exact series-resistance conductance
//     G = A / ( d/(2*k_owner) + R_c + d/(2*k_neighbour) )
// which reduces to the harmonic mean when R_c = 0 and is exact for linear
// temperature profiles in uniform cells. Boundary Dirichlet faces use
// G_b = k_owner * A / d_bc. Boundary fixed-gradient faces inject
// rhs += k_owner * A * g (positive g = heat flowing into the domain).
#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// ---------------------------------------------------------------------------
// Zone specification (user-declared, one per material)
// ---------------------------------------------------------------------------
struct ThermalZoneSpec {
    std::string name;                  // material name
    int zone_id = 0;                    // integer id referenced by the cell zone map
    double conductivity = 0.0;          // W/(m.K), finite and > 0
    double volumetric_power = 0.0;      // W/m^3, finite (uniform per zone)
    double contact_resistance = 0.0;   // m^2.K/W at interfaces with other zones
    bool solid = true;                  // reserved for CHT stage (issue #885)
};

// ---------------------------------------------------------------------------
// Assembled per-cell/per-face thermal state
// ---------------------------------------------------------------------------
struct ThermalZoneDiagnostics {
    std::size_t n_cells = 0;
    std::size_t n_zones = 0;
    std::size_t n_interface_faces = 0;  // internal faces between two zones
    double min_conductivity = std::numeric_limits<double>::infinity();
    double max_conductivity = 0.0;
    double max_contact_resistance = 0.0;
    std::map<std::string, std::size_t> cells_per_zone;
};

struct ThermalZoneModel {
    std::vector<double> cell_conductivity;         // per cell, W/(m.K)
    std::vector<double> cell_power;                // per cell, W/m^3
    std::vector<double> face_contact_resistance;   // per face, m^2.K/W (0 on
                                                    // same-zone and boundary faces)
    std::vector<int> cell_zone;                    // zone_id per cell
    ThermalZoneDiagnostics diagnostics;
    std::vector<ThermalZoneSpec> specs;
};

// ---------------------------------------------------------------------------
// Zone construction
// ---------------------------------------------------------------------------
// Builds the per-cell conductivity/power fields from a cell -> zone-id map.
// Throws std::invalid_argument on any inconsistent declaration.
inline ThermalZoneModel build_thermal_zones(
    const cfdx::core::Mesh& mesh,
    const std::vector<int>& cell_zone,
    const std::vector<ThermalZoneSpec>& specs)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    if (cell_zone.size() != nc)
        throw std::invalid_argument(
            "build_thermal_zones: cell zone map size does not match mesh cell count");
    if (specs.empty())
        throw std::invalid_argument("build_thermal_zones: no thermal zone declared");

    std::map<int, std::size_t> zone_index;
    for (std::size_t i = 0; i < specs.size(); ++i) {
        const auto& spec = specs[i];
        if (spec.name.empty())
            throw std::invalid_argument("build_thermal_zones: zone name must not be empty");
        if (!(spec.conductivity > 0.0) || !std::isfinite(spec.conductivity))
            throw std::invalid_argument(
                "build_thermal_zones: zone '" + spec.name +
                "' must declare a finite positive conductivity");
        if (!std::isfinite(spec.volumetric_power))
            throw std::invalid_argument(
                "build_thermal_zones: zone '" + spec.name +
                "' must declare a finite volumetric power");
        if (!(spec.contact_resistance >= 0.0) ||
            !std::isfinite(spec.contact_resistance))
            throw std::invalid_argument(
                "build_thermal_zones: zone '" + spec.name +
                "' must declare a non-negative finite contact resistance");
        if (!zone_index.emplace(spec.zone_id, i).second)
            throw std::invalid_argument(
                "build_thermal_zones: duplicate zone id " +
                std::to_string(spec.zone_id));
    }

    ThermalZoneModel zones;
    zones.specs = specs;
    zones.cell_conductivity.assign(nc, 0.0);
    zones.cell_power.assign(nc, 0.0);
    zones.cell_zone.assign(nc, 0);
    zones.face_contact_resistance.assign(mesh.n_faces(), 0.0);

    auto& diag = zones.diagnostics;
    diag.n_cells = nc;
    diag.n_zones = specs.size();

    for (std::size_t c = 0; c < nc; ++c) {
        const auto it = zone_index.find(cell_zone[c]);
        if (it == zone_index.end())
            throw std::invalid_argument(
                "build_thermal_zones: cell " + std::to_string(c) +
                " references undeclared zone id " + std::to_string(cell_zone[c]));
        const auto& spec = specs[it->second];
        zones.cell_conductivity[c] = spec.conductivity;
        zones.cell_power[c] = spec.volumetric_power;
        zones.cell_zone[c] = spec.zone_id;
        diag.cells_per_zone[spec.name] += 1;
    }

    // Contact resistance lives at internal faces between two *different*
    // declared zones. A same-zone face carries no contact resistance.
    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto nraw = own.neighbour(f);
        if (nraw < 0) continue;
        const std::size_t o = own.owner(f);
        const std::size_t n = static_cast<std::size_t>(nraw);
        if (zones.cell_zone[o] == zones.cell_zone[n]) continue;
        const double r_o =
            specs[zone_index.at(zones.cell_zone[o])].contact_resistance;
        const double r_n =
            specs[zone_index.at(zones.cell_zone[n])].contact_resistance;
        if (r_o > 0.0 && r_n > 0.0 && r_o != r_n)
            throw std::invalid_argument(
                "build_thermal_zones: inconsistent contact resistance at face " +
                std::to_string(f) + " between zones " +
                std::to_string(zones.cell_zone[o]) + " and " +
                std::to_string(zones.cell_zone[n]) + " (" +
                std::to_string(r_o) + " vs " + std::to_string(r_n) + ")");
        zones.face_contact_resistance[f] = r_o > 0.0 ? r_o : r_n;
        diag.n_interface_faces += 1;
        diag.max_contact_resistance =
            std::max(diag.max_contact_resistance, zones.face_contact_resistance[f]);
    }

    for (std::size_t c = 0; c < nc; ++c) {
        diag.min_conductivity = std::min(diag.min_conductivity,
                                        zones.cell_conductivity[c]);
        diag.max_conductivity = std::max(diag.max_conductivity,
                                        zones.cell_conductivity[c]);
    }
    return zones;
}

// ---------------------------------------------------------------------------
// Conduction assembly
// ---------------------------------------------------------------------------
// Assembles the steady conduction system
//     div( k grad(T) ) + q''' = 0
// with per-zone k, per-zone volumetric power and per-interface contact
// resistance. Boundary conditions follow ScalarBoundaryConditions semantics
// keyed by patch name.
inline ScalarEquation assemble_thermal_conduction(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ThermalZoneModel& zones,
    const ScalarBoundaryConditions& boundary_conditions = {})
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (zones.cell_conductivity.size() != nc ||
        zones.cell_power.size() != nc ||
        zones.face_contact_resistance.size() != nf ||
        zones.cell_zone.size() != nc)
        throw std::invalid_argument(
            "assemble_thermal_conduction: thermal zone model does not match mesh");

    std::vector<std::map<std::size_t, double>> rows(nc);
    std::vector<double> rhs(nc, 0.0);
    std::vector<double> diag(nc, 0.0);
    double min_internal_distance = std::numeric_limits<double>::infinity();
    double max_internal_distance = 0.0;
    double min_face_area = std::numeric_limits<double>::infinity();
    double max_face_area = 0.0;
    bool has_fixed_temperature = false;

    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t o = own.owner(f);
        if (o >= nc)
            throw std::runtime_error("assemble_thermal_conduction: invalid face owner");
        const auto nraw = own.neighbour(f);

        if (nraw >= 0) {
            const std::size_t n = static_cast<std::size_t>(nraw);
            if (n >= nc)
                throw std::runtime_error(
                    "assemble_thermal_conduction: invalid face neighbour");
            const double distance =
                (geometry.cell_centres[n] - geometry.cell_centres[o]).mag();
            const double area = geometry.face_area_vectors[f].mag();
            min_internal_distance = std::min(min_internal_distance, distance);
            max_internal_distance = std::max(max_internal_distance, distance);
            min_face_area = std::min(min_face_area, area);
            max_face_area = std::max(max_face_area, area);
            if (!(distance > 0.0) || !(area > 0.0) || !std::isfinite(distance) ||
                !std::isfinite(area))
                throw std::runtime_error(
                    "assemble_thermal_conduction: degenerate internal face");

            const double k_o = zones.cell_conductivity[o];
            const double k_n = zones.cell_conductivity[n];
            const double r_c = zones.face_contact_resistance[f];
            const double g = area /
                (distance / (2.0 * k_o) + r_c + distance / (2.0 * k_n));
            if (!(g > 0.0) || !std::isfinite(g))
                throw std::runtime_error(
                    "assemble_thermal_conduction: non-positive face conductance");

            diag[o] += g;
            diag[n] += g;
            rows[o][n] -= g;
            rows[n][o] -= g;
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;

            ScalarBoundaryCondition bc;  // default: adiabatic
            if (patch < mesh.boundary().n_patches()) {
                const auto& patch_name = mesh.boundary().patch(patch).name;
                const auto it = boundary_conditions.find(patch_name);
                if (it != boundary_conditions.end()) bc = it->second;
            }

            const double area = geometry.face_area_vectors[f].mag();
            const double distance = boundary_normal_distance(geometry, f, o);
            min_face_area = std::min(min_face_area, area);
            max_face_area = std::max(max_face_area, area);
            if (!(area > 0.0) || !(distance > 0.0))
                throw std::runtime_error(
                    "assemble_thermal_conduction: degenerate boundary face");

            const double k_o = zones.cell_conductivity[o];
            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                if (!std::isfinite(bc.value))
                    throw std::invalid_argument(
                        "assemble_thermal_conduction: non-finite boundary temperature");
                const double g_b = k_o * area / distance;
                diag[o] += g_b;
                rhs[o] += g_b * bc.value;
                has_fixed_temperature = true;
            } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                if (!std::isfinite(bc.gradient))
                    throw std::invalid_argument(
                        "assemble_thermal_conduction: non-finite boundary gradient");
                // Positive gradient = heat flowing into the domain.
                rhs[o] += k_o * area * bc.gradient;
            }
            // ZERO_GRADIENT: adiabatic, no contribution.
        }
    }

    // A purely Neumann conduction system (no fixed temperature anywhere)
    // has a constant null space: every diagonal entry stays positive on
    // a connected mesh, so the per-cell diagonal check below cannot
    // detect the singular case. Reject it explicitly.
    if (!has_fixed_temperature)
        throw std::runtime_error(
            "assemble_thermal_conduction: no fixed-temperature boundary "
            "found in the system; the pure conduction operator is singular "
            "(constant null space, temperature level unpinned). Declare at "
            "least one fixed-temperature patch");

    ScalarEquation eq{
        cfdx::core::SparseMatrix(nc, nc), cfdx::core::Vector(nc, 0.0), {}, 0.0};
    eq.diagonal.assign(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        rhs[c] += zones.cell_power[c] * geometry.cell_volumes[c];

        if (!(diag[c] > 0.0) || !std::isfinite(diag[c]))
            throw std::runtime_error(
                "assemble_thermal_conduction: cell " + std::to_string(c) +
                " has no fixed-temperature boundary in its stencil; the pure "
                "conduction system is singular (declare at least one "
                "fixed-temperature patch or a fixed-gradient patch)");

        std::vector<std::pair<std::size_t, double>> entries;
        entries.reserve(rows[c].size() + 1);
        for (const auto& [col, value] : rows[c]) entries.push_back({col, value});
        entries.push_back({c, diag[c]});
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        std::vector<std::pair<std::size_t, double>> merged;
        for (const auto& e : entries) {
            if (!merged.empty() && merged.back().first == e.first)
                merged.back().second += e.second;
            else
                merged.push_back(e);
        }
        for (const auto& [col, value] : merged) eq.matrix.push_back(c, col, value);
        eq.rhs(c) = rhs[c];
        eq.diagonal[c] = diag[c];
    }
    eq.matrix.finalize();
    eq.max_imbalance = 0.0;  // no convection in the conduction system
    eq.min_internal_distance = min_internal_distance;
    eq.max_internal_distance = max_internal_distance;
    eq.min_face_area = min_face_area;
    eq.max_face_area = max_face_area;
    return eq;
}

// ---------------------------------------------------------------------------
// Boundary flux audit
// ---------------------------------------------------------------------------
struct ThermalBoundaryFluxAudit {
    double boundary_flux_out = 0.0;               // W, positive leaving domain
    std::map<std::string, double> patch_flux_out;  // W, positive leaving domain
};

// Independent recomputation of the conductive boundary flux directly from
// the temperature field. This deliberately does not reuse the assembled
// matrix: it is the energy-balance audit for qualification (T1, issue #880).
inline ThermalBoundaryFluxAudit audit_thermal_boundary_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ThermalZoneModel& zones,
    const ScalarBoundaryConditions& boundary_conditions,
    const cfdx::core::Vector& temperature)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (temperature.size() != nc)
        throw std::invalid_argument(
            "audit_thermal_boundary_flux: temperature size does not match mesh");
    if (zones.cell_conductivity.size() != nc ||
        zones.face_contact_resistance.size() != nf)
        throw std::invalid_argument(
            "audit_thermal_boundary_flux: thermal zone model does not match mesh");

    ThermalBoundaryFluxAudit audit;
    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < nf; ++f) {
        const auto nraw = own.neighbour(f);
        if (nraw >= 0) continue;
        const std::size_t o = own.owner(f);
        const std::size_t patch = geometry.face_patch[f];
        if (patch >= mesh.boundary().n_patches()) continue;
        const auto& boundary_patch = mesh.boundary().patch(patch);
        if (boundary_patch.type == PatchType::EMPTY) continue;

        ScalarBoundaryCondition bc;  // default: adiabatic
        const auto it = boundary_conditions.find(boundary_patch.name);
        if (it != boundary_conditions.end()) bc = it->second;

        const double area = geometry.face_area_vectors[f].mag();
        const double distance = boundary_normal_distance(geometry, f, o);
        const double k_o = zones.cell_conductivity[o];
        double flux_out = 0.0;  // W, positive leaving the domain
        if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
            flux_out = k_o * area * (temperature(o) - bc.value) / distance;
        } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
            flux_out = -k_o * area * bc.gradient;
        }
        audit.boundary_flux_out += flux_out;
        audit.patch_flux_out[boundary_patch.name] += flux_out;
    }
    return audit;
}

// ---------------------------------------------------------------------------
// Solve + qualification result
// ---------------------------------------------------------------------------
struct ThermalConductionResult {
    bool converged = false;
    cfdx::core::SolverStatus status = cfdx::core::SolverStatus::MAX_ITER_REACHED;
    std::size_t iterations = 0;
    // Status and candidate_converged describe the linear candidate;
    // converged describes the actual returned temperature field.
    double residual_inf = 0.0;
    double residual_relative_inf = 0.0;
    bool candidate_converged = false;
    bool returned_state_converged = false;
    double min_temperature = 0.0;             // K
    double max_temperature = 0.0;             // K
    double total_source_power = 0.0;         // W
    double boundary_flux_out = 0.0;           // W
    std::map<std::string, double> patch_flux_out;  // W
    double energy_imbalance_relative = 0.0;  // |produced - net out| / scale
};

inline ThermalConductionResult solve_thermal_conduction(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ThermalZoneModel& zones,
    const ScalarBoundaryConditions& boundary_conditions,
    cfdx::core::Vector& temperature,
    const ScalarSolveControls& controls = {})
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    if (temperature.size() != nc)
        throw std::invalid_argument(
            "solve_thermal_conduction: temperature size does not match mesh");

    const ScalarEquation equation =
        assemble_thermal_conduction(mesh, geometry, zones, boundary_conditions);

    validate_scalar_controls(controls);
    const auto solve = solve_scalar_equation(equation, temperature, controls);

    ThermalConductionResult result;
    result.status = solve.status;
    result.iterations = solve.iterations;
    result.candidate_converged = solve.status == SolverStatus::CONVERGED;
    result.residual_inf = scalar_equation_residual_inf(equation, temperature);
    double rhs_norm_inf = 0.0;
    for (std::size_t i = 0; i < equation.rhs.size(); ++i)
        rhs_norm_inf = std::max(rhs_norm_inf, std::abs(equation.rhs(i)));
    result.residual_relative_inf =
        result.residual_inf / std::max(rhs_norm_inf, 1e-300);
    result.returned_state_converged =
        std::isfinite(result.residual_inf) &&
        std::isfinite(result.residual_relative_inf) &&
        result.residual_relative_inf <= controls.tolerance;
    result.converged = result.candidate_converged && result.returned_state_converged;

    for (std::size_t c = 0; c < nc; ++c) {
        if (c == 0) {
            result.min_temperature = temperature(c);
            result.max_temperature = temperature(c);
        } else {
            result.min_temperature = std::min(result.min_temperature, temperature(c));
            result.max_temperature = std::max(result.max_temperature, temperature(c));
        }
        result.total_source_power +=
            zones.cell_power[c] * geometry.cell_volumes[c];
    }

    const auto audit = audit_thermal_boundary_flux(
        mesh, geometry, zones, boundary_conditions, temperature);
    result.boundary_flux_out = audit.boundary_flux_out;
    result.patch_flux_out = audit.patch_flux_out;

    const double imbalance = result.total_source_power - result.boundary_flux_out;
    // Normalize by the gross energy throughput (sources plus per-patch flux
    // magnitudes), never by the net boundary flux: balanced cases (flux in ==
    // flux out, no source) have a net flux that is zero up to the linear
    // solver residual, and a net-flux normalization would turn that rounding
    // noise into an O(1) "imbalance".
    double gross_throughput = std::abs(result.total_source_power);
    for (const auto& entry : result.patch_flux_out)
        gross_throughput += std::abs(entry.second);
    result.energy_imbalance_relative =
        gross_throughput > 0.0 ? std::abs(imbalance) / gross_throughput : 0.0;
    return result;
}

}  // namespace cfdx::physics
