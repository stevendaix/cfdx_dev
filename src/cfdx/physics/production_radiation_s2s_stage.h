#pragma once

// Production S2S radiation stage (capability #885 stage 2, issue #885):
// couple surface-to-surface radiosity to the production energy stage on a
// frozen flow, one delayed-Picard pass.
//
// Contract (fiche methode docs/audits/11-rayonnement.md, APPROVED):
//   * One S2S surface = one boundary patch. A participating patch is a
//     mesh WALL patch (non-EMPTY) whose case BoundarySpec declares
//     raw_params["emissivity"]; the emissivity is parsed with the strict
//     case-double parser (whole text, finite) and must lie in (0, 1].
//     Declaring an emissivity on a non-WALL spec (or an EMPTY patch) is
//     rejected fail-closed instead of silently ignoring the declaration.
//   * At least two participating surfaces are required; a participant
//     cannot carry a retained Dirichlet temperature (wall_temperature or
//     inlet/fixed temperature): a radiosity surface whose temperature is
//     imposed by the case contradicts the wall-temperature feedback the
//     stage computes from the energy field.
//   * View factors come from the case: raw_params
//     ["s2s_view_factor_to_<target_patch>"], parsed with the strict
//     case-double parser, each in [0, 1]. A missing pair means 0. A view
//     factor toward the patch itself, toward a non-participant patch, or
//     declared on a patch without emissivity is rejected fail-closed.
//     Row closure sum_j F_ij = 1 is REQUIRED within 1e-10 (fiche 11, B5:
//     an enclosure the case does not close is a case-contract violation,
//     never silently patched), and reciprocity A_i F_ij = A_j F_ji is
//     validated with the production validator at machine tolerance.
//   * Surface areas are the sum of |face_area_vectors| over the patch's
//     boundary faces and must be finite and positive.
//   * The wall temperature of each surface is the |Sf|-weighted average of
//     the adjacent (owner) cell temperatures of the first-pass energy
//     field: this is the link between the radiation and the resolved
//     thermal field. The radiosity system is the production
//     solve_s2s_radiosity with machine controls (tolerance 1e-10, direct
//     dense solve).
//   * Sign convention (fiche 11, B3): the net radiative flux per surface is
//     q_i = J_i - G_i = eps_i (E_bi - G_i) in W/m^2, POSITIVE when the
//     surface loses radiant power. The flux injected into the energy
//     equation is -q_i per patch (W/m^2, positive INTO the domain), the
//     same convention as the declared heat_flux.
//   * Machine gates, BLOCKING: per-surface identity |J_i - G_i -
//     eps_i (E_bi - G_i)| <= 1e-10 * max(eps_i E_bi + G_i, 1 W/m^2), and
//     global enclosure balance |sum_i A_i q_i| <= 1e-10 *
//     max(emitted + absorbed totals, 1 W). With closure + reciprocity the
//     global balance closes exactly; a red gate can never yield converged.
//   * Orchestrated stage (run_production_s2s_energy_stage): energy pass 1
//     WITHOUT radiation -> T0; S2S solve on the T0 wall temperatures; the
//     net radiative fluxes are re-injected as extra wall heat fluxes
//     (-q_i per patch, summed with any declared heat_flux) -> energy pass
//     2 -> T1. This is ONE delayed-Picard pass: no claim of a converged
//     coupled fixed point is made (the external energy-radiation
//     iteration is step 3 of the integration plan, fiche 14 section 2
//     delegates it explicitly). max_wall_temperature_lag is reported as a
//     DIAGNOSTIC, never a gate. If pass 1 does not converge, or the
//     radiation gates are red, the stage is converged=false and no second
//     energy solve runs.

#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/core/field/field.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/radiation_s2s.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct ProductionS2SRadiationCase {
    std::vector<std::string> surface_names;  // patch names, patch order
    std::vector<double> areas;               // m^2, one per surface
    std::vector<double> emissivities;        // (0, 1], one per surface
    std::vector<double> view_factors;        // n*n row-major
};

namespace detail_s2s_stage {

inline const std::string kEmissivityKey = "emissivity";
inline const std::string kViewFactorPrefix = "s2s_view_factor_to_";

// Mirrors the Dirichlet detection of build_production_temperature_bcs: a
// retained fixed temperature (WALL_THERMAL wall_temperature or
// inlet/fixed temperature, zero meaning "not declared") contradicts the
// radiative wall-temperature feedback.
inline bool s2s_spec_has_retained_dirichlet(
    const cfdx::io::BoundarySpec& spec)
{
    if (spec.value_type == cfdx::io::BCValueType::WALL_THERMAL &&
        spec.wall_temperature != 0.0)
        return true;
    if ((spec.type == cfdx::io::BCType::INLET ||
         spec.value_type == cfdx::io::BCValueType::FIXED) &&
        spec.temperature != 0.0)
        return true;
    return false;
}

}  // namespace detail_s2s_stage

// Builds the S2S surface system from the case declaration. Fail-closed on
// every contract violation listed in the header comment; never invents an
// emissivity or a view factor.
inline ProductionS2SRadiationCase build_production_s2s_radiation_case(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::io::CaseSetup& setup)
{
    if (geometry.face_area_vectors.size() != mesh.n_faces() ||
        geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument(
            "production S2S radiation stage: finite-volume geometry does "
            "not match the mesh");

    using detail_s2s_stage::kEmissivityKey;
    using detail_s2s_stage::kViewFactorPrefix;

    // Classification pass: participants, emissivity violations, view
    // factors declared without emissivity.
    std::vector<std::size_t> participant_patches;
    std::vector<std::string> names;
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        const auto& patch = mesh.boundary().patch(p);
        const auto* spec = setup.find_boundary(patch.name);
        bool declares_emissivity = false;
        bool declares_view_factor = false;
        if (spec != nullptr) {
            declares_emissivity =
                spec->raw_params.find(kEmissivityKey) !=
                spec->raw_params.end();
            for (const auto& param : spec->raw_params) {
                if (param.first.rfind(kViewFactorPrefix, 0) == 0) {
                    declares_view_factor = true;
                    break;
                }
            }
        }
        if (!declares_emissivity && !declares_view_factor)
            continue;
        if (patch.type == cfdx::core::PatchType::EMPTY ||
            spec == nullptr || spec->type != cfdx::io::BCType::WALL) {
            if (declares_emissivity)
                throw std::invalid_argument(
                    "production S2S radiation stage: emissivity is only "
                    "supported on WALL patches, not on patch '" +
                    patch.name + "'");
        }
        if (declares_view_factor && !declares_emissivity)
            throw std::invalid_argument(
                "production S2S radiation stage: patch '" + patch.name +
                "' declares S2S view factors without declaring an "
                "emissivity; a non-participating patch cannot carry view "
                "factors");
        if (!declares_emissivity)
            continue;
        if (patch.type != cfdx::core::PatchType::WALL)
            throw std::invalid_argument(
                "production S2S radiation stage: participating patch '" +
                patch.name + "' must be a mesh WALL patch");
        if (detail_s2s_stage::s2s_spec_has_retained_dirichlet(*spec))
            throw std::invalid_argument(
                "production S2S radiation stage: participating patch '" +
                patch.name + "' declares a fixed temperature; an S2S "
                "surface cannot both participate in radiosity and have "
                "its wall temperature imposed by the case");
        participant_patches.push_back(p);
        names.push_back(patch.name);
    }
    if (participant_patches.size() < 2)
        throw std::invalid_argument(
            "production S2S radiation stage requires at least two "
            "participating WALL patches declaring an emissivity "
            "(declared: " +
            std::to_string(participant_patches.size()) + ")");

    const std::size_t n = participant_patches.size();

    // Areas: sum of |face area vectors| over the patch faces.
    std::vector<double> areas(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& patch = mesh.boundary().patch(participant_patches[i]);
        for (const auto face : patch.face_ids) {
            if (face >= mesh.n_faces())
                throw std::invalid_argument(
                    "production S2S radiation stage: patch '" +
                    patch.name + "' references an invalid face id");
            const auto& area = geometry.face_area_vectors[face];
            const double magnitude = std::sqrt(
                area.x * area.x + area.y * area.y + area.z * area.z);
            if (!std::isfinite(magnitude))
                throw std::invalid_argument(
                    "production S2S radiation stage: face area vector on "
                    "patch '" + patch.name + "' must be finite");
            areas[i] += magnitude;
        }
        if (!std::isfinite(areas[i]) || !(areas[i] > 0.0))
            throw std::invalid_argument(
                "production S2S radiation stage: surface '" + names[i] +
                "' must have a finite positive area");
    }

    // Emissivities: strict case-double parse, (0, 1].
    std::vector<double> emissivities(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto* spec =
            setup.find_boundary(mesh.boundary().patch(participant_patches[i]).name);
        const double emissivity = parse_case_double(
            spec->raw_params.at(kEmissivityKey),
            "emissivity on patch '" + names[i] + "'");
        if (!(emissivity > 0.0) || emissivity > 1.0)
            throw std::invalid_argument(
                "production S2S radiation stage: emissivity on patch '" +
                names[i] + "' must lie in (0, 1]");
        emissivities[i] = emissivity;
    }

    // View factors: strict parse, [0, 1], missing pair = 0, closure
    // required, reciprocity validated at machine tolerance (fiche 11 B5).
    std::vector<double> view_factors(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto* spec =
            setup.find_boundary(mesh.boundary().patch(participant_patches[i]).name);
        for (const auto& param : spec->raw_params) {
            if (param.first.rfind(kViewFactorPrefix, 0) != 0)
                continue;
            const std::string target = param.first.substr(
                kViewFactorPrefix.size());
            if (target == names[i])
                throw std::invalid_argument(
                    "production S2S radiation stage: patch '" + names[i] +
                    "' declares a view factor toward itself; an S2S "
                    "surface cannot see itself");
            bool target_is_participant = false;
            for (std::size_t j = 0; j < n; ++j)
                target_is_participant =
                    target_is_participant || target == names[j];
            if (!target_is_participant)
                throw std::invalid_argument(
                    "production S2S radiation stage: patch '" + names[i] +
                    "' declares a view factor toward '" + target +
                    "', which is not a participating S2S surface");
        }
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j)
                continue;
            const auto entry = spec->raw_params.find(
                kViewFactorPrefix + names[j]);
            if (entry == spec->raw_params.end())
                continue;  // missing pair = 0
            const double factor = parse_case_double(
                entry->second,
                "view factor '" + entry->first + "' on patch '" +
                names[i] + "'");
            if (!(factor >= 0.0) || factor > 1.0)
                throw std::invalid_argument(
                    "production S2S radiation stage: view factor '" +
                    entry->first + "' on patch '" + names[i] +
                    "' must lie in [0, 1]");
            view_factors[i * n + j] = factor;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            row_sum += view_factors[i * n + j];
        if (std::abs(row_sum - 1.0) > 1e-10)
            throw std::invalid_argument(
                "production S2S radiation stage: view factors of surface '" +
                names[i] + "' do not close (sum " +
                std::to_string(row_sum) +
                " != 1); the case must declare a closed enclosure");
    }
    validate_s2s_view_factors(view_factors, areas, 1e-10);

    ProductionS2SRadiationCase radiation_case;
    radiation_case.surface_names = names;
    radiation_case.areas = areas;
    radiation_case.emissivities = emissivities;
    radiation_case.view_factors = view_factors;
    return radiation_case;
}

// |Sf|-weighted average of the adjacent (owner) cell temperatures over
// each participating patch: the wall temperature the radiosity system
// sees. This is the only coupling from the resolved thermal field to the
// radiation.
inline std::vector<double> production_s2s_wall_temperatures(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ProductionS2SRadiationCase& radiation_case,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
{
    if (temperature.size() != mesh.n_cells() || temperature.dimension() != 1)
        throw std::invalid_argument(
            "production S2S radiation stage: temperature field does not "
            "match the mesh");
    const std::size_t n = radiation_case.surface_names.size();
    std::vector<double> wall_temperatures(n, 0.0);
    const auto& own = mesh.ownership();
    for (std::size_t i = 0; i < n; ++i) {
        const auto* patch = [&]() {
            for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p)
                if (mesh.boundary().patch(p).name ==
                    radiation_case.surface_names[i])
                    return &mesh.boundary().patch(p);
            return static_cast<const cfdx::core::Patch*>(nullptr);
        }();
        if (patch == nullptr)
            throw std::invalid_argument(
                "production S2S radiation stage: surface '" +
                radiation_case.surface_names[i] +
                "' does not match a mesh patch");
        double weighted = 0.0;
        double total_weight = 0.0;
        for (const auto face : patch->face_ids) {
            if (face >= mesh.n_faces() || own.neighbour(face) >= 0)
                throw std::invalid_argument(
                    "production S2S radiation stage: surface '" +
                    radiation_case.surface_names[i] +
                    "' must be built from boundary faces");
            const auto& area = geometry.face_area_vectors[face];
            const double weight = std::sqrt(
                area.x * area.x + area.y * area.y + area.z * area.z);
            const std::size_t cell = static_cast<std::size_t>(own.owner(face));
            if (cell >= mesh.n_cells())
                throw std::invalid_argument(
                    "production S2S radiation stage: boundary face owner "
                    "is invalid");
            weighted += weight * temperature(cell);
            total_weight += weight;
        }
        if (!(total_weight > 0.0) || !std::isfinite(weighted))
            throw std::invalid_argument(
                "production S2S radiation stage: surface '" +
                radiation_case.surface_names[i] +
                "' has no boundary faces to average");
        const double wall_temperature = weighted / total_weight;
        if (!std::isfinite(wall_temperature) || !(wall_temperature > 0.0))
            throw std::invalid_argument(
                "production S2S radiation stage: wall temperature of "
                "surface '" + radiation_case.surface_names[i] +
                "' must be finite and positive (kelvin)");
        wall_temperatures[i] = wall_temperature;
    }
    return wall_temperatures;
}

struct ProductionS2SRadiationResult {
    // Blocking verdict: radiosity solve AND per-surface identity gate AND
    // global enclosure balance gate. A red machine gate can never yield
    // converged.
    bool converged = false;
    bool balance_passed = false;
    bool surface_identity_passed = false;
    std::vector<std::string> surface_names;
    std::vector<double> areas;              // m^2
    std::vector<double> emissivities;
    std::vector<double> wall_temperatures;   // K (T0-based)
    std::vector<double> radiosity;           // W/m^2, J_i
    std::vector<double> irradiation;         // W/m^2, G_i
    std::vector<double> net_flux;            // W/m^2, q_i = J_i - G_i,
                                             // positive = surface loses
    std::vector<double> net_power;           // W, A_i * q_i
    // Flux injected into the energy equation per patch (W/m^2, positive
    // into the domain): -q_i, the same convention as declared heat_flux.
    std::map<std::string, double> injected_heat_flux;
    double total_net_power = 0.0;     // W, sum A_i q_i (closes to ~0)
    double balance_scale = 0.0;       // W, emitted + absorbed, floor 1 W
    double balance_normalized = 0.0;  // |total| / scale
    double max_surface_identity_normalized = 0.0;
    std::size_t surface_count = 0;
};

// Solves the production S2S radiosity system on the wall temperatures
// derived from the given temperature field, with blocking machine gates.
inline ProductionS2SRadiationResult solve_production_s2s_radiation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::io::CaseSetup& setup,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
{
    if (!production_s2s_radiation_requested(setup))
        throw std::invalid_argument(
            "solve_production_s2s_radiation: radiation_model 'S2S' is not "
            "declared; the production S2S radiation stage only runs on a "
            "declared S2S case");
    if (setup.energy_model != "conjugate")
        throw std::invalid_argument(
            "solve_production_s2s_radiation: S2S radiation is coupled into "
            "the production energy stage and requires energy_model "
            "'conjugate' (declared: '" + setup.energy_model + "')");
    if (temperature.size() != mesh.n_cells() || temperature.dimension() != 1)
        throw std::invalid_argument(
            "solve_production_s2s_radiation: temperature field does not "
            "match the mesh");
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        if (!std::isfinite(temperature(c)) || !(temperature(c) > 0.0))
            throw std::invalid_argument(
                "solve_production_s2s_radiation: temperature must be "
                "finite and positive (kelvin)");

    const auto radiation_case =
        build_production_s2s_radiation_case(mesh, geometry, setup);
    const auto wall_temperatures = production_s2s_wall_temperatures(
        mesh, geometry, radiation_case, temperature);

    S2SControls controls;
    controls.tolerance = 1e-10;
    controls.max_iterations = 1;
    controls.ambient_irradiation = 0.0;
    const auto solution = solve_s2s_radiosity(
        radiation_case.areas, radiation_case.emissivities,
        wall_temperatures, radiation_case.view_factors, controls);

    ProductionS2SRadiationResult result;
    result.surface_names = radiation_case.surface_names;
    result.areas = radiation_case.areas;
    result.emissivities = radiation_case.emissivities;
    result.wall_temperatures = wall_temperatures;
    result.radiosity = solution.radiosity;
    result.irradiation = solution.irradiation;
    result.net_flux.resize(radiation_case.surface_names.size(), 0.0);
    result.net_power.resize(radiation_case.surface_names.size(), 0.0);
    result.surface_count = radiation_case.surface_names.size();

    // Machine gates (fiche 11): per-surface identity q_i = eps_i (E_bi -
    // G_i) holds by construction of the radiosity system; it is verified
    // here as an independent blocking gate, normalized by that surface's
    // own gross radiant scale (emitted + irradiation, floor 1 W/m^2).
    double emitted_total = 0.0;
    double absorbed_total = 0.0;
    for (std::size_t i = 0; i < result.surface_count; ++i) {
        const double q = solution.radiosity[i] - solution.irradiation[i];
        result.net_flux[i] = q;
        result.net_power[i] = radiation_case.areas[i] * q;
        result.total_net_power += result.net_power[i];
        const double black_body_power =
            STEFAN_BOLTZMANN *
            std::pow(wall_temperatures[i], 4);
        const double identity =
            std::abs(q - radiation_case.emissivities[i] *
                              (black_body_power - solution.irradiation[i]));
        const double surface_scale =
            std::max(radiation_case.emissivities[i] * black_body_power +
                         solution.irradiation[i],
                     1.0);
        result.max_surface_identity_normalized =
            std::max(result.max_surface_identity_normalized,
                     identity / surface_scale);
        emitted_total += radiation_case.areas[i] *
                         radiation_case.emissivities[i] * black_body_power;
        absorbed_total += radiation_case.areas[i] *
                          radiation_case.emissivities[i] *
                          solution.irradiation[i];
        // Injected flux: positive into the domain, same convention as the
        // declared heat_flux (fiche 11 B3 + fiche 14).
        result.injected_heat_flux[radiation_case.surface_names[i]] = -q;
    }

    // Global enclosure balance: with closure + reciprocity the enclosure
    // loses nothing, so sum A_i q_i closes to machine zero. Blocking.
    result.balance_scale =
        std::max(emitted_total + absorbed_total, 1.0);
    result.balance_normalized =
        std::abs(result.total_net_power) / result.balance_scale;
    const double machine_tolerance = 1e-10;
    result.surface_identity_passed =
        std::isfinite(result.max_surface_identity_normalized) &&
        result.max_surface_identity_normalized <= machine_tolerance;
    result.balance_passed = std::isfinite(result.balance_normalized) &&
                             result.balance_normalized <= machine_tolerance;
    result.converged =
        solution.converged && result.balance_passed &&
        result.surface_identity_passed;
    return result;
}

struct ProductionS2SEnergyStageResult {
    // Blocking verdict of the orchestrated delayed-Picard pass: energy
    // pass 1 converged AND S2S radiation gates green AND energy pass 2
    // (with the radiative fluxes re-injected) converged. No claim is made
    // about the coupled fixed point; max_wall_temperature_lag is a
    // diagnostic only (the external iteration is integration step 3).
    bool converged = false;
    bool energy_pass1_converged = false;
    bool s2s_radiation_converged = false;
    bool energy_pass2_converged = false;
    ProductionEnergyStageResult energy;  // final (pass 2) energy result
    ProductionS2SRadiationResult radiation;
    double max_wall_temperature_lag = 0.0;  // K, diagnostic, not gated
};

// Orchestrates one delayed-Picard pass of the production energy stage
// coupled with S2S radiation on a frozen flow: (1) energy pass 1 without
// radiation -> T0; (2) S2S radiosity on the T0 wall temperatures with
// blocking machine gates; (3) re-injection of the net radiative fluxes as
// extra wall heat fluxes -> energy pass 2 -> T1. If pass 1 does not
// converge or the radiation gates are red, the stage returns
// converged=false and no second energy solve runs.
inline ProductionS2SEnergyStageResult run_production_s2s_energy_stage(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::io::CaseSetup& setup,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
{
    if (!production_s2s_radiation_requested(setup))
        throw std::invalid_argument(
            "run_production_s2s_energy_stage: radiation_model 'S2S' is not "
            "declared; the production S2S energy stage only runs on a "
            "declared S2S case");
    if (setup.energy_model != "conjugate")
        throw std::invalid_argument(
            "run_production_s2s_energy_stage: S2S radiation is coupled into "
            "the production energy stage and requires energy_model "
            "'conjugate' (declared: '" + setup.energy_model + "')");

    ProductionS2SEnergyStageResult result;
    const auto pass1 =
        run_production_energy_stage(mesh, geometry, setup, mass_flux,
                                    temperature);
    result.energy = pass1;
    result.energy_pass1_converged = pass1.converged;
    if (!pass1.converged)
        return result;

    const auto radiation =
        solve_production_s2s_radiation(mesh, geometry, setup, temperature);
    result.radiation = radiation;
    result.s2s_radiation_converged = radiation.converged;
    if (!radiation.converged)
        return result;

    // Re-injection: the net radiative flux of each surface enters the
    // domain as -q_i (W/m^2, positive into the domain), summed with any
    // declared heat_flux by the energy BC builder.
    const auto extra_flux = radiation.injected_heat_flux;
    const auto pass2 = run_production_energy_stage(
        mesh, geometry, setup, mass_flux, temperature, extra_flux);
    result.energy = pass2;
    result.energy_pass2_converged = pass2.converged;
    result.converged = pass2.converged;

    // Delayed-Picard diagnostic: how much the wall temperatures moved
    // between the radiation input (T0) and the re-solved field (T1). Not
    // a gate; the external energy-radiation iteration is step 3.
    const auto radiation_case =
        build_production_s2s_radiation_case(mesh, geometry, setup);
    const auto wall_t1 = production_s2s_wall_temperatures(
        mesh, geometry, radiation_case, temperature);
    for (std::size_t i = 0; i < radiation.wall_temperatures.size(); ++i)
        result.max_wall_temperature_lag =
            std::max(result.max_wall_temperature_lag,
                     std::abs(wall_t1[i] - radiation.wall_temperatures[i]));
    return result;
}

}  // namespace cfdx::physics
