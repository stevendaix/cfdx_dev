#pragma once

// T2 production buoyancy stage: inverse energy->flow coupling through the
// Boussinesq momentum source (capability #885 step 4, issue #885).
//
// Contract:
//   * energy_model gate: the stage only runs for "buoyant". A buoyant case
//     does NOT solve a frozen carrier: the flow solve moves inside this
//     outer loop, so the momentum equation sees the temperature-dependent
//     Boussinesq source and the energy equation sees the responding mass
//     flux.
//   * Physics (fiche methode docs/audits/10-thermique-cht.md, APPROVED
//     2026-10-12 PR #861): rho = rho0 [1 - beta (T - T0)] in the gravity
//     term ONLY, valid for beta dT << 1. The stage builds the per-cell
//     acceleration
//         a_i = beta (T_ref - T_i) g_vec                 [m/s^2]
//     which is the Boussinesq perturbation of the gravity force per unit
//     mass (hot pushes against gravity). The constant hydrostatic part
//     rho0 g_vec is deliberately NOT included: the coupled solve is
//     written in modified pressure (the constant part is absorbed by the
//     pressure field), the same convention as the existing BoussinesqModel
//     (src/cfdx/physics/boussinesq.h). T == T_ref gives exactly zero.
//   * Required declarations, validated fail-closed BEFORE any solve: a
//     gravity_vector (exactly three strict finite numbers, non-zero), the
//     Boussinesq thermal expansion coefficient beta > 0 in the single
//     material's extra_properties["thermal_expansion_coefficient"] (same
//     strict whole-string parse as volumetric_power, fiche 14), a finite
//     positive ref_temperature, exactly one material, and NO case
//     body_force (the outer loop owns the momentum source; a uniform force
//     cannot coexist with the cellwise buoyancy source). The energy
//     material contract (fiche 14) is pre-flighted here too so an invalid
//     case fails at load instead of after the first carrier solve.
//   * Outer iteration: k = 1..200 fail-closed. Pass k builds the cellwise
//     source from the relaxed temperature seed, solves the carrier flow
//     with the declared initialization on k = 1 and a warm start (Provided)
//     afterwards, runs the production energy stage (fiche 14) on the new
//     authoritative face mass flux - through the nested S2S
//     energy-radiation iteration (fiche 11, step 3) when radiation is
//     declared - then relaxes the seed T += 0.25 (T_out - T). Both inner
//     stages are reused UNCHANGED.
//   * Blocking criteria on the UNRELAXED response, never the relaxed path:
//     stationarity of T, U and p between successive responses, each
//     normalized by max(1, max|value|), all <= 1e-10. k = 1 can never
//     converge (no previous response); two identical responses are the
//     proof of the coupled fixed point. Response stationarity plus the
//     relaxation update imply the seed equals the response at the passing
//     iteration (a lagging seed would move the next source and break
//     pressure stationarity), so the verdict is the exact coupled fixed
//     point, never an intermediate pass. INCONCLUSIVE is never PASS.
//   * A red inner stage (carrier flow, energy, S2S) fails the whole loop
//     closed; an exhausted budget is fail-closed. A buoyant case whose
//     initial temperature sits exactly at T_ref starts from a zero source
//     and may honestly stagnate on the first carrier solve (lecon #663):
//     declare an initial temperature away from T_ref unless that is the
//     intent.

#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

namespace detail_buoyancy_stage {

inline constexpr double kOuterRelaxation = 0.25;
inline constexpr std::size_t kOuterMaxIterations = 200;
inline constexpr double kOuterFixedPointTolerance = 1e-10;

// Normalized max change between successive scalar responses:
// max|current - previous| / max(1, max|current|). The floor keeps the
// criterion absolute near zero; the normalization keeps it scale-free on
// large fields. A non-finite value fails the gate closed.
inline double normalized_max_change(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& current,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& previous)
{
    double max_abs_current = 0.0;
    double max_abs_change = 0.0;
    for (std::size_t cell = 0; cell < current.size(); ++cell) {
        const double value = current(cell);
        const double change = value - previous(cell);
        if (!std::isfinite(value) || !std::isfinite(change))
            return std::numeric_limits<double>::infinity();
        max_abs_current = std::max(max_abs_current, std::abs(value));
        max_abs_change = std::max(max_abs_change, std::abs(change));
    }
    return max_abs_change / std::max(max_abs_current, 1.0);
}

// Same normalized max change over all three velocity components.
inline double normalized_max_change_vector(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& current,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& previous)
{
    double max_abs_current = 0.0;
    double max_abs_change = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const auto& current_component = current.component_data(d);
        const auto& previous_component = previous.component_data(d);
        for (std::size_t cell = 0; cell < current_component.size(); ++cell) {
            const double value = current_component[cell];
            const double change = value - previous_component[cell];
            if (!std::isfinite(value) || !std::isfinite(change))
                return std::numeric_limits<double>::infinity();
            max_abs_current = std::max(max_abs_current, std::abs(value));
            max_abs_change = std::max(max_abs_change, std::abs(change));
        }
    }
    return max_abs_change / std::max(max_abs_current, 1.0);
}

}  // namespace detail_buoyancy_stage

// Routing predicate of the production driver (issue #885 step 4): true when
// the case declares the buoyant energy model, whatever the radiation model
// is.
inline bool production_buoyancy_requested(const cfdx::io::CaseSetup& setup)
{
    return setup.energy_model == "buoyant";
}

// Strict gravity declaration: exactly three whitespace-separated finite
// numbers. A numeric prefix without consuming the whole string, a wrong
// token count, nan/inf and a zero vector are all rejected (same
// whole-string contract as parse_case_double; a zero gravity would
// silently run the buoyant case uncoupled).
inline cfdx::core::Vec3 parse_gravity_vector(const std::string& declared)
{
    std::istringstream stream(declared);
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!(stream >> x >> y >> z))
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must be three "
            "numbers, e.g. \"0 0 -9.81\" (declared: '" + declared + "')");
    std::string trailing;
    if (stream >> trailing)
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must have exactly "
            "three components (declared: '" + declared + "')");
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector components must be "
            "finite (declared: '" + declared + "')");
    if (x == 0.0 && y == 0.0 && z == 0.0)
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must be non-zero; a "
            "buoyant case without gravity would silently run uncoupled");
    return cfdx::core::Vec3{x, y, z};
}

// Boussinesq thermal expansion coefficient from the single material's
// extra_properties, with the strict whole-string case parse (fiche 14
// volumetric_power convention). beta must be strictly positive and finite.
inline double production_buoyancy_beta(const cfdx::io::CaseSetup& setup)
{
    if (setup.materials.size() != 1)
        throw std::invalid_argument(
            "production buoyancy stage: energy_model 'buoyant' requires "
            "exactly one material declaring the Boussinesq thermal "
            "expansion coefficient (declared materials: " +
            std::to_string(setup.materials.size()) + ")");
    const auto& material = setup.materials.front();
    const auto beta_it =
        material.extra_properties.find("thermal_expansion_coefficient");
    if (beta_it == material.extra_properties.end())
        throw std::invalid_argument(
            "production buoyancy stage: energy_model 'buoyant' requires the "
            "material extra_properties 'thermal_expansion_coefficient'");
    const double beta = parse_case_double(
        beta_it->second,
        "material extra_properties 'thermal_expansion_coefficient'");
    if (!(beta > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: material extra_properties "
            "'thermal_expansion_coefficient' must be strictly positive "
            "(declared: '" + beta_it->second + "')");
    return beta;
}

// Fail-closed validation of the full buoyant declaration. Runs at case
// load (before any expensive computation) and again at the entry of the
// outer loop. Non-buoyant cases pass through untouched.
inline void validate_production_buoyancy_model(
    const cfdx::io::CaseSetup& setup)
{
    if (!production_buoyancy_requested(setup))
        return;
    if (setup.gravity_vector.find_first_not_of(" \t\r\n") ==
        std::string::npos)
        throw std::invalid_argument(
            "production buoyancy stage: energy_model 'buoyant' requires a "
            "declared gravity_vector (e.g. \"0 0 -9.81\")");
    (void)parse_gravity_vector(setup.gravity_vector);
    (void)production_buoyancy_beta(setup);
    if (!std::isfinite(setup.ref_temperature) ||
        !(setup.ref_temperature > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: ref_temperature must be finite and "
            "positive (kelvin) for energy_model 'buoyant'");
    if (!setup.body_force.empty())
        throw std::invalid_argument(
            "production buoyancy stage: energy_model 'buoyant' declares a "
            "case body_force; the buoyancy outer loop owns the momentum "
            "source, and a uniform body force cannot coexist with the "
            "cellwise Boussinesq source");
    // Pre-flight the energy material contract (fiche 14) so an invalid
    // material fails the run at case load, before the first carrier solve,
    // instead of after it.
    (void)build_production_energy_case(setup);
}

// Per-cell Boussinesq buoyancy acceleration a_i = beta (T_ref - T_i) g_vec.
// T == T_ref gives exactly zero; hot cells push against gravity. The
// temperature must be finite and positive (kelvin).
inline std::vector<cfdx::core::Vec3> build_cell_buoyancy_acceleration(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    double beta, double reference_temperature,
    const cfdx::core::Vec3& gravity)
{
    if (temperature.dimension() != 1)
        throw std::invalid_argument(
            "production buoyancy stage: temperature field must be scalar");
    if (!std::isfinite(beta) || !(beta > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: thermal expansion coefficient must "
            "be finite and strictly positive");
    if (!std::isfinite(reference_temperature) ||
        !(reference_temperature > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: reference temperature must be "
            "finite and positive (kelvin)");
    std::vector<cfdx::core::Vec3> acceleration(temperature.size());
    for (std::size_t cell = 0; cell < temperature.size(); ++cell) {
        const double t = temperature(cell);
        if (!std::isfinite(t) || !(t > 0.0))
            throw std::invalid_argument(
                "production buoyancy stage: temperature must be finite and "
                "positive (kelvin)");
        const double scaled = beta * (reference_temperature - t);
        acceleration[cell] = cfdx::core::Vec3{
            scaled * gravity.x, scaled * gravity.y, scaled * gravity.z};
    }
    return acceleration;
}

struct ProductionBuoyancyIteration {
    std::size_t iteration = 0;
    bool flow_converged = false;
    bool energy_converged = false;
    // Inner S2S energy-radiation loop verdict on this pass; false when S2S
    // is not declared.
    bool s2s_converged = false;
    std::size_t flow_iterations = 0;
    std::size_t energy_iterations = 0;
    // Stationarity of the UNRELAXED response against the previous
    // iteration's response; infinity while no previous response exists.
    double max_temperature_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_velocity_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_pressure_change_normalized =
        std::numeric_limits<double>::infinity();
};

struct ProductionBuoyancyOuterLoopResult {
    // Blocking verdict, see the contract above: green inner stages on the
    // final pass AND T/U/p response stationarity within 1e-10. A red inner
    // stage or an exhausted budget can never yield converged.
    bool converged = false;
    bool energy_ran = false;
    bool s2s_ran = false;
    std::size_t outer_iterations = 0;
    double max_temperature_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_velocity_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_pressure_change_normalized =
        std::numeric_limits<double>::infinity();
    std::vector<ProductionBuoyancyIteration> history;
    IncompressibleSolveResult flow;               // final carrier solve
    ProductionEnergyStageResult energy;           // final energy pass
    ProductionS2SOuterLoopResult s2s_outer_result;  // nested loop when S2S
    // Final coupled temperature (the unrelaxed response of the final
    // pass; the diagnostics artifact even when the loop fails closed
    // after at least one green pass).
    cfdx::core::Field<double, cfdx::core::Location::CELL> temperature;
};

inline ProductionBuoyancyOuterLoopResult run_production_buoyant_outer_loop(
    const cfdx::core::Mesh& mesh,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& base_controls,
    const cfdx::io::CaseSetup& setup,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const std::string& restart_path = {})
{
    if (!production_buoyancy_requested(setup))
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: the buoyancy outer loop "
            "only runs for energy_model 'buoyant'; running it under "
            "another declared model is a contract violation, not a silent "
            "fallback");
    validate_production_buoyancy_model(setup);
    if (U.dimension() != 3 || U.size() != mesh.n_cells() ||
        p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: invalid fields");

    ProductionBuoyancyOuterLoopResult result;
    const bool s2s_requested = production_s2s_radiation_requested(setup);
    const cfdx::core::Vec3 gravity = parse_gravity_vector(
        setup.gravity_vector);
    const double beta = production_buoyancy_beta(setup);
    const double reference_temperature = setup.ref_temperature;

    const FvGeometry geometry = build_fv_geometry(mesh);

    // Relaxed temperature seed (driver convention, stage 2): the declared
    // initial temperature when present, else the case reference
    // temperature. Zero means "not declared".
    double seed_temperature = setup.ref_temperature;
    if (setup.has_initial_condition &&
        setup.initial_condition.temperature != 0.0)
        seed_temperature = setup.initial_condition.temperature;
    if (!std::isfinite(seed_temperature) || !(seed_temperature > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: initial temperature must be finite "
            "and positive (kelvin)");
    cfdx::core::Field<double, cfdx::core::Location::CELL> relaxed_temperature(
        mesh.n_cells(), "temperature", "K", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
        relaxed_temperature(cell) = seed_temperature;

    // Trackers of the previous iteration's UNRELAXED response.
    bool have_previous_response = false;
    cfdx::core::Field<double, cfdx::core::Location::CELL> previous_temperature(
        mesh.n_cells(), "temperature_previous", "K", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> previous_U(
        mesh.n_cells(), "U_previous", "m/s", 3);
    cfdx::core::Field<double, cfdx::core::Location::CELL> previous_p(
        mesh.n_cells(), "p_previous", "Pa", 1);

    for (std::size_t iteration = 1;
         iteration <= detail_buoyancy_stage::kOuterMaxIterations; ++iteration) {
        ProductionBuoyancyIteration record;
        record.iteration = iteration;
        result.outer_iterations = iteration;

        // Momentum source of this pass, built from the relaxed seed. The
        // k = 1 pass keeps the declared initialization; every later pass
        // warm-starts the carrier from the current fields (Provided keeps
        // the caller state as-is).
        IncompressibleSolverControls controls = base_controls;
        controls.cell_body_force = build_cell_buoyancy_acceleration(
            relaxed_temperature, beta, reference_temperature, gravity);
        if (iteration > 1)
            controls.initialization.mode = InitializationMode::Provided;

        const auto flow = solve_steady_incompressible(
            mesh, U, p, velocity_bcs, pressure_bcs, controls, restart_path);
        result.flow = flow;
        record.flow_converged = flow.converged;
        record.flow_iterations = flow.iterations;
        if (!flow.converged) {
            result.history.push_back(record);
            return result;  // fail-closed: no energy pass on a diverged
                            // carrier
        }

        // Energy response on the new authoritative mass flux. Both stages
        // write their solution into relaxed_temperature, so the seed is
        // copied first for the relaxation update below.
        const auto seed = relaxed_temperature;
        if (s2s_requested) {
            const auto s2s = run_production_s2s_energy_outer_loop(
                mesh, geometry, setup, flow.authoritative_mass_flux,
                relaxed_temperature);
            result.s2s_outer_result = s2s;
            result.energy = s2s.energy;
            result.s2s_ran = true;
            record.energy_converged = s2s.energy.converged;
            record.s2s_converged = s2s.converged;
            record.energy_iterations = s2s.energy.iterations;
            if (!s2s.energy.converged || !s2s.converged) {
                result.history.push_back(record);
                return result;  // fail-closed: red inner energy/S2S stage
            }
        } else {
            const auto energy = run_production_energy_stage(
                mesh, geometry, setup, flow.authoritative_mass_flux,
                relaxed_temperature);
            result.energy = energy;
            record.energy_converged = energy.converged;
            record.energy_iterations = energy.iterations;
            if (!energy.converged) {
                result.history.push_back(record);
                return result;  // fail-closed: red inner energy stage
            }
        }
        result.energy_ran = true;
        // The unrelaxed response of this pass is the artifact a
        // diagnostics consumer sees, even when the loop later fails
        // closed.
        result.temperature = relaxed_temperature;

        // Blocking criteria on the UNRELAXED response against the previous
        // response. k = 1 has no previous response and can never converge;
        // two identical responses are the proof of the coupled fixed
        // point.
        if (have_previous_response) {
            record.max_temperature_change_normalized =
                detail_buoyancy_stage::normalized_max_change(
                    relaxed_temperature, previous_temperature);
            record.max_velocity_change_normalized =
                detail_buoyancy_stage::normalized_max_change_vector(
                    U, previous_U);
            record.max_pressure_change_normalized =
                detail_buoyancy_stage::normalized_max_change(
                    p, previous_p);
            result.max_temperature_change_normalized =
                record.max_temperature_change_normalized;
            result.max_velocity_change_normalized =
                record.max_velocity_change_normalized;
            result.max_pressure_change_normalized =
                record.max_pressure_change_normalized;
            const double tolerance =
                detail_buoyancy_stage::kOuterFixedPointTolerance;
            if (std::isfinite(record.max_temperature_change_normalized) &&
                record.max_temperature_change_normalized <= tolerance &&
                std::isfinite(record.max_velocity_change_normalized) &&
                record.max_velocity_change_normalized <= tolerance &&
                std::isfinite(record.max_pressure_change_normalized) &&
                record.max_pressure_change_normalized <= tolerance) {
                result.converged = true;  // exact coupled fixed point
                result.history.push_back(record);
                return result;
            }
        }

        previous_temperature = relaxed_temperature;
        previous_U = U;
        previous_p = p;
        have_previous_response = true;

        // Under-relaxed seed update for the next pass:
        // T(k+1) = T(k) + omega (T_out(k) - T(k)).
        const double omega = detail_buoyancy_stage::kOuterRelaxation;
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
            relaxed_temperature(cell) =
                seed(cell) +
                omega * (relaxed_temperature(cell) - seed(cell));

        result.history.push_back(record);
    }

    // Budget exhausted without a stationary response: fail-closed. The
    // verdict can never be an intermediate pass.
    return result;
}

}  // namespace cfdx::physics
