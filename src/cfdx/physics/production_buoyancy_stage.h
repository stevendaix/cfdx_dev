#pragma once

// Production buoyancy stage (capability #885 step 4): couple the energy
// field back into the incompressible carrier through the Boussinesq
// buoyancy source (fiche methode docs/audits/10-thermique-cht.md,
// APPROVED; carrier fiche 04; energy stage fiche 14; S2S stage fiche 11
// when radiation is declared).
//
// Contract (fiche methode 10, Boussinesq approximation):
//   * The buoyancy source is a MOMENTUM source only:
//     a_i = -beta (T_i - T_ref) g_vec per unit mass (m/s^2), with g_vec
//     the declared gravity vector (case gravity_vector, e.g. "0 0 -9.81")
//     and beta the material thermal expansion coefficient. Density
//     variations appear in the gravity term ONLY
//     (rho = rho0 [1 - beta (T - T_ref)]); mass continuity, inertia and
//     the energy equation keep the constant declared rho0, so the
//     approximation is valid for beta (T - T_ref) << 1. Full
//     variable-density coupling is out of scope of the approved contract.
//   * The constant rho0 g_vec part is absorbed in the modified pressure
//     (the carrier carries no hydrostatic term): at T = T_ref the source
//     is exactly zero and the buoyant run reproduces the buoyancy-free
//     carrier. This matches the scalar magnitude convention of
//     boussinesq.h (buoyancy_acceleration = gravity * beta * (T - T_ref),
//     positive upward for a hot cell).
//   * Activation: energy_model = "buoyant" (case schema). Fail-closed
//     declarations required up front: a non-zero gravity_vector, a
//     material thermal expansion coefficient (extra_properties
//     ["thermal_expansion_coefficient"], strict case-double parse,
//     finite and strictly positive), a finite positive ref_temperature
//     and exactly one declared material. A buoyant case must NOT declare
//     a uniform body_force: the buoyancy source is the only body force
//     and the cellwise solver channel replaces the uniform one.
//   * Orchestration: run_production_buoyant_outer_loop iterates the
//     primitives of the approved stages UNCHANGED - the incompressible
//     carrier solve (fiche 04) with the cellwise body force built from
//     the under-relaxed temperature, then the production energy stage
//     (fiche 14), or the stage-3 energy-radiation outer loop (fiche 11)
//     when S2S is declared, on the new authoritative mass flux. Every
//     outer iteration therefore costs one full carrier solve (warm
//     started from the previous outer state after the first iteration).
//   * Blocking criteria (verdict = the exact coupled fixed point;
//     INCONCLUSIVE is never PASS): every inner stage green (carrier
//     converged, energy gates, S2S fixed point when declared) AND
//     (a) the non-relaxed energy response reproduces the relaxed input
//     that built the buoyancy source (normalized source-input lag
//     <= 1e-10: the flow was solved with the source of the input, so a
//     response equal to its input certifies the coupled fixed point
//     exactly, the direct analog of the stage-3 radiative fixed-point
//     residual), AND (b) the non-relaxed responses are stationary to
//     1e-10 normalized (temperature in K, velocity in m/s, pressure in
//     Pa; gross scale = the max magnitude of the response field, floor
//     1.0 in each unit). At least two outer responses must exist: a
//     single pass can never claim the coupled fixed point. Budget 200
//     outer iterations, fail-closed; under-relaxation omega = 0.25 on
//     the temperature the buoyancy source sees (path only, the verdict
//     reads the non-relaxed response).
//
// No tolerance, gate, oracle or budget of the inner stages is modified
// by this header; the energy equation, the carrier algorithms and the
// radiosity system are consumed exactly as approved.

#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "cfdx/physics/steady_incompressible_solver.h"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

namespace detail_buoyancy_stage {
inline constexpr double kBuoyancyRelaxation = 0.25;
inline constexpr std::size_t kBuoyancyMaxOuterIterations = 200;
inline constexpr double kBuoyancyFixedPointTolerance = 1e-10;
inline constexpr const char* kThermalExpansionKey =
    "thermal_expansion_coefficient";
}  // namespace detail_buoyancy_stage

inline bool production_buoyancy_requested(const cfdx::io::CaseSetup& setup)
{
    return setup.energy_model == "buoyant";
}

// Strict gravity_vector parsing: exactly three whitespace-separated
// finite numbers ("0 0 -9.81"). Each token is parsed with the strict
// case-double parser (whole text, finite), and trailing text is
// rejected instead of being silently ignored.
inline cfdx::core::Vec3 parse_gravity_vector(const std::string& text)
{
    std::istringstream stream(text);
    double components[3] = {0.0, 0.0, 0.0};
    for (int index = 0; index < 3; ++index) {
        std::string token;
        if (!(stream >> token))
            throw std::invalid_argument(
                "production buoyancy stage: gravity_vector must declare "
                "exactly 3 components (e.g. \"0 0 -9.81\")");
        components[index] =
            parse_case_double(token, "gravity_vector component");
    }
    std::string extra;
    if (stream >> extra)
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must declare "
            "exactly 3 components; trailing text is rejected ('" + extra +
            "')");
    const cfdx::core::Vec3 gravity{components[0], components[1],
                                   components[2]};
    if (!std::isfinite(gravity.x) || !std::isfinite(gravity.y) ||
        !std::isfinite(gravity.z))
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector components must be "
            "finite");
    return gravity;
}

// Boussinesq coupling data parsed from the case, fail-closed.
struct ProductionBuoyancyCoupling {
    double beta = 0.0;                       // 1/K, strictly positive
    double reference_temperature = 0.0;       // K, finite positive
    cfdx::core::Vec3 gravity{0.0, 0.0, 0.0};  // m/s^2, non-zero
};

// Fail-closed validation of a buoyant case. A no-op for every other
// energy model: each model owns its own validation gate, and this one
// only owns the Boussinesq coupling declarations.
inline void validate_production_buoyancy_model(const cfdx::io::CaseSetup& setup)
{
    if (!production_buoyancy_requested(setup))
        return;
    if (!setup.body_force.empty())
        throw std::invalid_argument(
            "production buoyancy stage: a buoyant case must not declare a "
            "uniform body_force; the Boussinesq buoyancy source is the only "
            "body force (specify either, not both)");
    if (setup.materials.size() != 1)
        throw std::invalid_argument(
            "production buoyancy stage: a buoyant case requires exactly one "
            "declared material carrying the thermal expansion coefficient");
    const auto& material = setup.materials.front();
    const auto beta_it =
        material.extra_properties.find(detail_buoyancy_stage::kThermalExpansionKey);
    if (beta_it == material.extra_properties.end())
        throw std::invalid_argument(
            std::string("production buoyancy stage: the material must "
                        "declare '") +
            detail_buoyancy_stage::kThermalExpansionKey +
            "' in extra_properties (Boussinesq beta, 1/K) instead of "
            "silently running without buoyancy");
    const double beta = parse_case_double(
        beta_it->second, detail_buoyancy_stage::kThermalExpansionKey);
    if (!(beta > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: thermal_expansion_coefficient must "
            "be finite and strictly positive (1/K)");
    if (!std::isfinite(setup.ref_temperature) ||
        !(setup.ref_temperature > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: ref_temperature must be finite and "
            "positive (kelvin)");
    bool has_non_whitespace = false;
    for (const char character : setup.gravity_vector) {
        if (!std::isspace(static_cast<unsigned char>(character))) {
            has_non_whitespace = true;
            break;
        }
    }
    if (!has_non_whitespace)
        throw std::invalid_argument(
            "production buoyancy stage: a buoyant case must declare "
            "gravity_vector (e.g. \"0 0 -9.81\") instead of silently "
            "running without gravity");
    const cfdx::core::Vec3 gravity =
        parse_gravity_vector(setup.gravity_vector);
    const double magnitude = std::sqrt(
        gravity.x * gravity.x + gravity.y * gravity.y +
        gravity.z * gravity.z);
    if (!(magnitude > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must be non-zero; a "
            "buoyant case without gravity is a contract violation, not a "
            "silent isothermal run");
}

inline ProductionBuoyancyCoupling production_buoyancy_coupling(
    const cfdx::io::CaseSetup& setup)
{
    validate_production_buoyancy_model(setup);
    ProductionBuoyancyCoupling coupling;
    const auto& material = setup.materials.front();
    coupling.beta = parse_case_double(
        material.extra_properties.at(detail_buoyancy_stage::kThermalExpansionKey),
        detail_buoyancy_stage::kThermalExpansionKey);
    coupling.reference_temperature = setup.ref_temperature;
    coupling.gravity = parse_gravity_vector(setup.gravity_vector);
    return coupling;
}

// Boussinesq buoyancy acceleration per unit mass (m/s^2):
// a_i = -beta (T_i - T_ref) g_vec. A hot cell (T_i > T_ref) accelerates
// against gravity (upward for a downward gravity vector), the exact
// vector form of boussinesq.h's scalar convention. A cell at T_ref
// carries exactly zero source.
inline void build_cell_buoyancy_acceleration(
    std::vector<cfdx::core::Vec3>& acceleration,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    const ProductionBuoyancyCoupling& coupling)
{
    const std::size_t n_cells = temperature.size();
    acceleration.assign(n_cells, cfdx::core::Vec3{0.0, 0.0, 0.0});
    for (std::size_t cell = 0; cell < n_cells; ++cell) {
        const double delta =
            temperature(cell) - coupling.reference_temperature;
        acceleration[cell] = cfdx::core::Vec3{
            -coupling.beta * delta * coupling.gravity.x,
            -coupling.beta * delta * coupling.gravity.y,
            -coupling.beta * delta * coupling.gravity.z};
    }
}

struct ProductionBuoyancyIterationRecord {
    std::size_t outer_iteration = 0;
    std::size_t carrier_iterations = 0;
    std::size_t energy_iterations = 0;
    std::size_t s2s_outer_iterations = 0;
    // Normalized non-relaxed stationarity of the responses (infinity on
    // the first iteration: a single pass can never prove the fixed point).
    double max_temperature_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_velocity_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_pressure_change_normalized =
        std::numeric_limits<double>::infinity();
    // Normalized lag between the non-relaxed energy response and the
    // relaxed input that built the buoyancy source: at the coupled fixed
    // point the response reproduces its input exactly (the direct analog
    // of the stage-3 radiative fixed-point residual, applied to the
    // buoyancy source).
    double max_source_input_lag_normalized =
        std::numeric_limits<double>::infinity();
    bool carrier_converged = false;
    bool energy_converged = false;
    bool s2s_converged = false;
};

struct ProductionBuoyancyOuterLoopResult {
    // Blocking verdict: the exact coupled fixed point (all inner stages
    // green AND the non-relaxed responses stationary). False means no
    // claim, never a degraded PASS.
    bool converged = false;
    std::size_t outer_iterations = 0;
    // Last accepted carrier solve (the driver's flow artifacts and the
    // restart checkpoint reflect this final outer state).
    IncompressibleSolverResult flow;
    // Last energy response (and the S2S fixed point when declared).
    ProductionEnergyStageResult energy;
    bool energy_ran = false;
    bool s2s_ran = false;
    ProductionS2SOuterLoopResult s2s;
    double max_temperature_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_velocity_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_pressure_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_source_input_lag_normalized =
        std::numeric_limits<double>::infinity();
    // Empty when converged; the honest reason otherwise.
    std::string failure_reason;
    std::vector<ProductionBuoyancyIterationRecord> outer_history;
};

// The flow-energy buoyancy outer loop. U and p are the carrier solution
// fields (initial state per base_controls.initialization on the first
// outer iteration, warm start Provided afterwards; final coupled state
// on return). temperature receives the final non-relaxed energy
// response (sized to the mesh by the loop, seeded like the frozen-flow
// stage: declared initial temperature, else the case reference
// temperature).
inline ProductionBuoyancyOuterLoopResult run_production_buoyant_outer_loop(
    const cfdx::core::Mesh& mesh,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& base_controls,
    const cfdx::io::CaseSetup& setup,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    const std::string& restart_path = "")
{
    if (!production_buoyancy_requested(setup))
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: the buoyancy loop only "
            "runs for energy_model 'buoyant'; running it under another "
            "declared model is a contract violation, not a silent "
            "fallback");
    validate_production_buoyancy_model(setup);

    ProductionBuoyancyOuterLoopResult result;
    const ProductionBuoyancyCoupling coupling =
        production_buoyancy_coupling(setup);
    const double omega = detail_buoyancy_stage::kBuoyancyRelaxation;
    const std::size_t max_outer_iterations =
        detail_buoyancy_stage::kBuoyancyMaxOuterIterations;
    const double tolerance = detail_buoyancy_stage::kBuoyancyFixedPointTolerance;

    const std::size_t n_cells = mesh.n_cells();
    if (U.size() != n_cells || U.dimension() != 3 ||
        p.size() != n_cells || p.dimension() != 1)
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: solution fields do not "
            "match the mesh");
    if (temperature.size() != n_cells || temperature.dimension() != 1)
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: temperature field does not "
            "match the mesh");

    // Initial temperature: the declared initial condition, else the case
    // reference temperature (stage-1 convention, same as the frozen-flow
    // path; zero means "not declared").
    double initial_temperature = setup.ref_temperature;
    if (setup.has_initial_condition &&
        setup.initial_condition.temperature != 0.0)
        initial_temperature = setup.initial_condition.temperature;
    if (!std::isfinite(initial_temperature) ||
        !(initial_temperature > 0.0))
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: initial temperature must "
            "be finite and positive (kelvin)");
    for (std::size_t cell = 0; cell < n_cells; ++cell)
        temperature(cell) = initial_temperature;

    const FvGeometry geometry = build_fv_geometry(mesh);
    const bool s2s_requested = production_s2s_radiation_requested(setup);

    // The under-relaxed temperature the buoyancy source sees (path only;
    // the verdict reads the non-relaxed response).
    cfdx::core::Field<double, cfdx::core::Location::CELL> relaxed_temperature(
        n_cells, "buoyancy_relaxed_temperature", "K", 1);
    for (std::size_t cell = 0; cell < n_cells; ++cell)
        relaxed_temperature(cell) = initial_temperature;
    // Previous non-relaxed responses for the stationarity gates.
    cfdx::core::Field<double, cfdx::core::Location::CELL> previous_temperature(
        n_cells, "buoyancy_previous_temperature", "K", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> previous_U(
        n_cells, "buoyancy_previous_U", "m/s", 3);
    cfdx::core::Field<double, cfdx::core::Location::CELL> previous_p(
        n_cells, "buoyancy_previous_p", "Pa", 1);
    bool have_previous_response = false;

    std::vector<cfdx::core::Vec3> acceleration;

    for (std::size_t outer = 1; outer <= max_outer_iterations; ++outer) {
        ProductionBuoyancyIterationRecord record;
        record.outer_iteration = outer;

        // 1. Buoyancy source from the relaxed temperature.
        build_cell_buoyancy_acceleration(acceleration, relaxed_temperature,
                                         coupling);

        // 2. Carrier solve with the cellwise body force (the channel
        //    replaces the uniform one; the first outer iteration follows
        //    the declared initialization, afterwards the loop warm starts
        //    from the previous outer state).
        IncompressibleSolverControls controls = base_controls;
        controls.cell_body_force = acceleration;
        if (outer > 1)
            controls.initialization.mode = InitializationMode::Provided;
        const IncompressibleSolverResult flow = solve_steady_incompressible(
            mesh, U, p, velocity_bcs, pressure_bcs, controls, restart_path);
        result.flow = flow;
        record.carrier_iterations = flow.iterations;
        record.carrier_converged = flow.converged;
        if (!flow.converged) {
            result.outer_iterations = outer;
            result.failure_reason =
                "carrier flow did not converge inside the buoyancy "
                "iteration at outer iteration " + std::to_string(outer) +
                " (INCONCLUSIVE is never PASS)";
            result.outer_history.push_back(record);
            return result;
        }

        // 3. Energy response on the new authoritative flux. The energy
        //    field is re-seeded from the relaxed temperature so the
        //    response is a deterministic function of the relaxed input;
        //    the stage-3 energy-radiation outer loop (fiche 11) carries
        //    the radiative fixed point for THIS carrier when S2S is
        //    declared.
        for (std::size_t cell = 0; cell < n_cells; ++cell)
            temperature(cell) = relaxed_temperature(cell);
        if (s2s_requested) {
            const ProductionS2SOuterLoopResult s2s =
                run_production_s2s_energy_outer_loop(
                    mesh, geometry, setup, flow.authoritative_mass_flux,
                    temperature);
            result.s2s = s2s;
            result.s2s_ran = true;
            result.energy = s2s.energy;
            result.energy_ran = true;
            record.energy_iterations = s2s.energy.iterations;
            record.s2s_outer_iterations = s2s.outer_iterations;
            record.energy_converged = s2s.energy.converged;
            record.s2s_converged = s2s.converged;
            if (!s2s.energy.converged) {
                result.outer_iterations = outer;
                result.failure_reason =
                    "energy stage gates are red inside the buoyancy "
                    "iteration at outer iteration " + std::to_string(outer) +
                    " (INCONCLUSIVE is never PASS)";
                result.outer_history.push_back(record);
                return result;
            }
            if (!s2s.converged) {
                result.outer_iterations = outer;
                result.failure_reason =
                    "the energy-radiation outer loop did not reach its "
                    "fixed point inside the buoyancy iteration at outer "
                    "iteration " + std::to_string(outer) +
                    " (INCONCLUSIVE is never PASS)";
                result.outer_history.push_back(record);
                return result;
            }
        } else {
            const ProductionEnergyStageResult energy =
                run_production_energy_stage(
                    mesh, geometry, setup, flow.authoritative_mass_flux,
                    temperature);
            result.energy = energy;
            result.energy_ran = true;
            record.energy_iterations = energy.iterations;
            record.energy_converged = energy.converged;
            if (!energy.converged) {
                result.outer_iterations = outer;
                result.failure_reason =
                    "energy stage gates are red inside the buoyancy "
                    "iteration at outer iteration " + std::to_string(outer) +
                    " (INCONCLUSIVE is never PASS)";
                result.outer_history.push_back(record);
                return result;
            }
        }

        // 4. Fixed-point gates on the NON-RELAXED response:
        //    (a) the energy response must reproduce the relaxed input
        //        that built the buoyancy source (the flow was solved with
        //        the source of the input, so response == input certifies
        //        the exact coupled fixed point, the direct analog of the
        //        stage-3 radiative fixed-point residual);
        //    (b) the responses must be stationary (normalized by the gross
        //        response scale, floor 1.0 in each unit). The first
        //        response has nothing to be stationary against, so a
        //        single pass can never converge.
        {
            double max_source_input_lag = 0.0;
            double response_scale = 1.0;
            for (std::size_t cell = 0; cell < n_cells; ++cell) {
                max_source_input_lag =
                    std::max(max_source_input_lag,
                             std::abs(temperature(cell) -
                                      relaxed_temperature(cell)));
                response_scale =
                    std::max(response_scale, std::abs(temperature(cell)));
            }
            record.max_source_input_lag_normalized =
                max_source_input_lag / response_scale;
            result.max_source_input_lag_normalized =
                record.max_source_input_lag_normalized;
        }
        if (have_previous_response) {
            double max_temperature_change = 0.0;
            double temperature_scale = 1.0;
            double max_velocity_change = 0.0;
            double velocity_scale = 1.0;
            double max_pressure_change = 0.0;
            double pressure_scale = 1.0;
            for (std::size_t cell = 0; cell < n_cells; ++cell) {
                max_temperature_change = std::max(
                    max_temperature_change,
                    std::abs(temperature(cell) - previous_temperature(cell)));
                temperature_scale = std::max(
                    temperature_scale, std::abs(temperature(cell)));
                for (std::size_t d = 0; d < 3; ++d) {
                    max_velocity_change = std::max(
                        max_velocity_change,
                        std::abs(U.component_data(d)[cell] -
                                 previous_U.component_data(d)[cell]));
                    velocity_scale = std::max(
                        velocity_scale,
                        std::abs(U.component_data(d)[cell]));
                }
                max_pressure_change = std::max(
                    max_pressure_change,
                    std::abs(p(cell) - previous_p(cell)));
                pressure_scale =
                    std::max(pressure_scale, std::abs(p(cell)));
            }
            record.max_temperature_change_normalized =
                max_temperature_change / temperature_scale;
            record.max_velocity_change_normalized =
                max_velocity_change / velocity_scale;
            record.max_pressure_change_normalized =
                max_pressure_change / pressure_scale;
            result.max_temperature_change_normalized =
                record.max_temperature_change_normalized;
            result.max_velocity_change_normalized =
                record.max_velocity_change_normalized;
            result.max_pressure_change_normalized =
                record.max_pressure_change_normalized;
            if (record.max_temperature_change_normalized <= tolerance &&
                record.max_velocity_change_normalized <= tolerance &&
                record.max_pressure_change_normalized <= tolerance &&
                record.max_source_input_lag_normalized <= tolerance) {
                result.outer_iterations = outer;
                result.converged = true;  // exact coupled fixed point
                result.outer_history.push_back(record);
                return result;
            }
        }

        // 5. Under-relaxed feedback for the next carrier solve
        //    (path only): T(k+1) = T(k) + omega (R(T(k)) - T(k)).
        for (std::size_t cell = 0; cell < n_cells; ++cell)
            relaxed_temperature(cell) +=
                omega * (temperature(cell) - relaxed_temperature(cell));

        // 6. Remember the non-relaxed response for the next stationarity
        //    check.
        previous_temperature = temperature;
        previous_U = U;
        previous_p = p;
        have_previous_response = true;
        result.outer_iterations = outer;
        result.outer_history.push_back(record);
    }

    // Budget exhausted: converged=false, full history, no claim.
    result.failure_reason =
        "buoyancy outer iteration budget exhausted (" +
        std::to_string(max_outer_iterations) +
        " iterations) without the coupled fixed point";
    return result;
}

}  // namespace cfdx::physics
