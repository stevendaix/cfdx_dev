#pragma once

// Production buoyancy stage (capability #885 step 4, issue #885): the
// inverse energy->flow coupling. With energy_model "buoyant" the carrier
// flow is no longer solved once and frozen: it is iterated INSIDE an
// external flow-energy loop. Each outer iteration re-solves the steady
// incompressible carrier with a cellwise Boussinesq buoyancy source built
// from the under-relaxed temperature field, then re-solves the energy
// stage on the carrier's new authoritative mass flux (through the stage-3
// S2S energy-radiation iteration when the case declares S2S radiation).
// The verdict is the converged COUPLED FIXED POINT, never an intermediate
// Picard pass; INCONCLUSIVE is never PASS.
//
// Contract (fiche methode docs/audits/10-thermique-cht.md, APPROVED:
// Boussinesq approximation; fiche 14 for the reused energy stage; fiche
// 04 for the pressure-velocity carrier):
//   * The density perturbation exists ONLY in the gravity term
//     (rho = rho_ref [1 - beta (T - T_ref)], valid for beta*dT << 1). The
//     momentum source is the per-cell acceleration
//         a_i = -beta (T_i - T_ref) g_vec        [m/s^2]
//     with g_vec the DECLARED case gravity_vector (direction and
//     magnitude, e.g. "0 0 -9.81"): hot fluid with downward gravity
//     accelerates upward. The constant rho_ref*g_vec part is absorbed in
//     the modified pressure, matching the existing BoussinesqModel::
//     buoyancy_acceleration scalar convention (src/cfdx/physics/
//     boussinesq.h, unchanged). Full variable-density (low-Mach)
//     coupling is OUT OF SCOPE: it is not covered by any approved fiche
//     and would change the incompressible continuity constraint itself.
//   * The channel into the carrier is the solver's optional cellwise body
//     force (IncompressibleSolverControls::cell_body_force, acceleration
//     per unit mass, multiplied by the constant density internally,
//     replacing the uniform body_force when supplied). A buoyant case
//     must NOT declare a uniform case body_force: the buoyancy source is
//     the only body force, and the solver rejects a uniform and a
//     cellwise body force supplied together.
//   * Fail-closed declarations, validated before any flow solve: exactly
//     one material, beta declared through that material's
//     extra_properties ("thermal_expansion_coefficient", strict finite
//     numeric literal, strictly positive), a declared 3-component finite
//     gravity_vector (a zero vector is rejected: it would make the
//     coupling identically zero), and a finite positive ref_temperature.
//   * Outer-loop controls are fixed and documented (no case-facing knobs
//     invented): omega = 0.25 under-relaxes the temperature seen by the
//     buoyancy source (same damping as the stage-3 radiative feedback),
//     the budget is 200 outer iterations, and the acceptance tolerance
//     is 1e-10 normalized. Blocking criteria are SEPARATE and evaluated
//     on the UNRELAXED responses: temperature, velocity and pressure
//     stationarity between successive outer iterations, each normalized
//     by its own gross scale (max(1, max|field|)). The first outer
//     iteration has no previous response and can never prove the fixed
//     point.
//   * Inner gates are the unchanged production gates: the carrier must
//     converge, and the energy stage (fiche 14) and, when declared, the
//     S2S energy-radiation outer loop (fiche 11) must be green, on every
//     outer iteration. A diverged carrier or a red inner stage stops
//     the loop fail-closed with converged=false; the loop never
//     iterates on a state that is not evidence.
//   * Each outer iteration is a FULL carrier solve (warm-started from
//     the previous outer iteration's accepted state after the first
//     pass). This is the documented cost of the explicit flow-energy
//     iteration; no inner flow iteration budget is ever relaxed to
//     compensate.
//   * Convergence-at-entry: at the coupled fixed point a warm carrier
//     re-solve can end STAGNATED with every residual at the
//     machine-precision floor (diag run 39633631822: pressure,
//     continuity and corrected-flux residuals exactly 0.0, momentum
//     equation residual RHS-normalized ~1e-17) - its relative
//     criterion can no longer be re-proven against the collapsed
//     reference of the same solve. Such a solve is ACCEPTED: the
//     entry velocity and pressure fields are restored bit-exactly,
//     the previously accepted carrier (authoritative mass flux) is
//     reused, and the OUTER loop owns the coupled fixed-point
//     criteria. A genuine stagnation far above the machine floor
//     remains a fail-closed stop. Accepted-at-entry iterations are
//     recorded per iteration and on the result as
//     flow_converged_at_entry.

#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "cfdx/io/cfdx_io/case_schema.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

namespace detail_buoyancy_stage {

// Material extra_properties key carrying the Boussinesq thermal expansion
// coefficient beta (1/K). Parsed with the production strict case-double
// parser, never defaulted.
inline const std::string kThermalExpansionKey =
    "thermal_expansion_coefficient";

// Fixed outer-loop controls, documented and fail-closed (no case-facing
// knobs invented): omega damps the explicit buoyancy feedback path (same
// value as the stage-3 radiative feedback), the budget bounds the loop,
// and the acceptance tolerance matches the machine-grade coupled gates.
inline constexpr double kOuterRelaxation = 0.25;
inline constexpr std::size_t kOuterMaxIterations = 200;
inline constexpr double kOuterFixedPointTolerance = 1e-10;
// Convergence-at-entry machine floor (contract above, diag run
// 39633631822): a warm carrier re-solve at the coupled fixed point
// stalls with pressure, continuity and corrected-flux residuals
// exactly 0.0 and a momentum equation residual RHS-normalized
// ~1e-17. A genuine stagnation stays orders of magnitude above
// this floor (NACA carrier: continuity_normalized ~1e-5), so the
// distinction is unambiguous.
inline constexpr double kCarrierMachineEquilibriumResidual = 1e-14;
inline constexpr double kCarrierMachineEquilibriumMomentumRhs = 1e-12;

}  // namespace detail_buoyancy_stage

// Routing predicate of the production driver (issue #885 step 4): true
// when the case declares the buoyant energy model.
inline bool production_buoyancy_requested(const cfdx::io::CaseSetup& setup)
{
    return setup.energy_model == "buoyant";
}

// Strict gravity vector parse: the WHOLE text must be exactly three
// whitespace-separated finite numbers ("0 0 -9.81"). A missing or extra
// token, a non-numeric token or a non-finite value is rejected whole
// instead of being silently truncated or defaulted.
inline cfdx::core::Vec3 parse_gravity_vector(const std::string& text)
{
    std::istringstream stream(text);
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!(stream >> x) || !(stream >> y) || !(stream >> z))
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must be exactly 3 "
            "whitespace-separated finite numbers, declared: '" + text + "'");
    std::string extra;
    if (stream >> extra)
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector must be exactly 3 "
            "whitespace-separated finite numbers, declared: '" + text + "'");
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw std::invalid_argument(
            "production buoyancy stage: gravity_vector components must be "
            "finite, declared: '" + text + "'");
    return cfdx::core::Vec3{x, y, z};
}

// The parsed Boussinesq coupling, built once per run from the case
// declaration. Every member is fail-closed validated.
struct ProductionBuoyancyCoupling {
    // Boussinesq thermal expansion coefficient (1/K), strictly positive.
    double beta = 0.0;
    // Reference temperature T_ref (K) of the Boussinesq linearisation,
    // strictly positive.
    double reference_temperature = 0.0;
    // Declared gravity vector (m/s^2): direction AND magnitude, exactly
    // as the case declares it.
    cfdx::core::Vec3 gravity{0.0, 0.0, 0.0};
};

inline ProductionBuoyancyCoupling make_production_buoyancy_coupling(
    const cfdx::io::CaseSetup& setup)
{
    ProductionBuoyancyCoupling coupling;
    if (setup.materials.size() != 1)
        throw std::invalid_argument(
            "production buoyancy stage: the buoyant model requires exactly "
            "one material carrying the fluid properties and the "
            "'thermal_expansion_coefficient' (declared materials: " +
            std::to_string(setup.materials.size()) + ")");
    const auto& material = setup.materials.front();
    const auto beta_it =
        material.extra_properties.find(
            detail_buoyancy_stage::kThermalExpansionKey);
    if (beta_it == material.extra_properties.end())
        throw std::invalid_argument(
            "production buoyancy stage: the buoyant model requires the "
            "material extra_properties key "
            "'thermal_expansion_coefficient' (Boussinesq beta, 1/K); no "
            "default is invented");
    coupling.beta = parse_case_double(
        beta_it->second, detail_buoyancy_stage::kThermalExpansionKey);
    if (!(coupling.beta > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: thermal_expansion_coefficient must "
            "be strictly positive (declared: '" + beta_it->second + "')");
    if (!std::isfinite(setup.ref_temperature) ||
        !(setup.ref_temperature > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: ref_temperature must be finite and "
            "positive (kelvin)");
    coupling.reference_temperature = setup.ref_temperature;
    const auto& gravity_text = setup.gravity_vector;
    const bool gravity_unspecified =
        gravity_text.find_first_not_of(" \t\r\n") == std::string::npos;
    if (gravity_unspecified)
        throw std::invalid_argument(
            "production buoyancy stage: the buoyant model requires a "
            "declared gravity_vector (e.g. '0 0 -9.81'); buoyancy without "
            "gravity is not silently run as an isothermal flow");
    coupling.gravity = parse_gravity_vector(gravity_text);
    if (coupling.gravity.x == 0.0 && coupling.gravity.y == 0.0 &&
        coupling.gravity.z == 0.0)
        throw std::invalid_argument(
            "production buoyancy stage: a declared zero gravity_vector "
            "makes the buoyancy coupling identically zero; declare the "
            "physical vector or run the non-buoyant energy model");
    if (!setup.body_force.empty())
        throw std::invalid_argument(
            "production buoyancy stage: a buoyant case must not declare a "
            "uniform body_force; the Boussinesq buoyancy source is the "
            "only body force of the momentum equation, and the solver "
            "rejects a uniform and a cellwise body force supplied "
            "together");
    return coupling;
}

// Case-load validation of the buoyant model (no-op for any other declared
// energy model): all buoyancy declarations are checked fail-closed before
// any expensive flow computation.
inline void validate_production_buoyancy_model(
    const cfdx::io::CaseSetup& setup)
{
    if (!production_buoyancy_requested(setup))
        return;
    (void)make_production_buoyancy_coupling(setup);
}

// Per-cell Boussinesq buoyancy acceleration (fiche methode 10): the
// density perturbation exists ONLY in the gravity term, as
//     a_i = -beta (T_i - T_ref) g_vec        [m/s^2]
// (hot fluid, downward gravity => upward acceleration; the constant
// rho_ref*g_vec part is absorbed in the modified pressure). The output
// vector is the solver's cellwise body-force channel: acceleration per
// unit mass, multiplied by the constant density internally.
inline void build_cell_buoyancy_acceleration(
    std::vector<cfdx::core::Vec3>& acceleration,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    const ProductionBuoyancyCoupling& coupling)
{
    const std::size_t n_cells = temperature.size();
    acceleration.assign(n_cells, cfdx::core::Vec3{0.0, 0.0, 0.0});
    for (std::size_t cell = 0; cell < n_cells; ++cell) {
        const double value = temperature(cell);
        if (!std::isfinite(value) || !(value > 0.0))
            throw std::invalid_argument(
                "production buoyancy stage: temperature must be finite and "
                "positive (kelvin)");
        const double perturbation =
            coupling.beta * (value - coupling.reference_temperature);
        acceleration[cell] = cfdx::core::Vec3{
            -perturbation * coupling.gravity.x,
            -perturbation * coupling.gravity.y,
            -perturbation * coupling.gravity.z};
    }
}

struct ProductionBuoyancyOuterIteration {
    std::size_t iteration = 0;
    bool flow_converged = false;
    // True only on an outer iteration accepted at entry (see the
    // convergence-at-entry contract above): the warm carrier re-solve
    // stalled at the machine-precision floor, the entry state was
    // restored and the previously accepted carrier reused.
    bool flow_converged_at_entry = false;
    bool energy_converged = false;
    bool s2s_ran = false;
    bool s2s_converged = false;
    std::size_t flow_iterations = 0;
    std::size_t energy_iterations = 0;
    std::size_t s2s_outer_iterations = 0;
    // Stationarity of the UNRELAXED responses vs the previous outer
    // iteration, each normalized by its own gross scale
    // (max(1, max|field|)). Infinity on the first iteration: a single
    // Picard pass can never prove the coupled fixed point.
    double max_temperature_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_velocity_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_pressure_change_normalized =
        std::numeric_limits<double>::infinity();
};

struct ProductionBuoyancyOuterLoopResult {
    // Blocking verdict of the external flow-energy buoyancy iteration,
    // see the contract above: green inner gates on every iteration AND
    // temperature, velocity and pressure stationarity of the unrelaxed
    // responses within 1e-10 normalized. A red inner stage or an
    // exhausted budget can never yield converged.
    bool converged = false;
    bool fixed_point_temperature_passed = false;
    bool fixed_point_velocity_passed = false;
    bool fixed_point_pressure_passed = false;
    // True when the recorded energy response belongs to the recorded
    // final flow state (the last outer iteration completed both solves).
    bool energy_ran = false;
    std::size_t outer_iterations = 0;
    double max_temperature_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_velocity_change_normalized =
        std::numeric_limits<double>::infinity();
    double max_pressure_change_normalized =
        std::numeric_limits<double>::infinity();
    std::vector<ProductionBuoyancyOuterIteration> history;
    // Convergence-at-entry bookkeeping: true when the FINAL outer
// iteration's warm carrier re-solve was accepted at the entry state
    // (entry fields restored, previously accepted carrier reused); the
    // outer loop owns the coupled fixed-point criteria.
    bool flow_converged_at_entry = false;
    // Final outer iteration's carrier solve, energy response and, when
    // declared, energy-radiation fixed point on that carrier.
    IncompressibleSolveResult flow;
    ProductionEnergyStageResult energy;
    bool s2s_ran = false;
    ProductionS2SOuterLoopResult s2s;
    // Final energy temperature response (K).
    cfdx::core::Field<double, cfdx::core::Location::CELL> temperature;
    // The validated coupling the loop ran with.
    ProductionBuoyancyCoupling coupling;
};

// Gross-scale normalizer shared by the three stationarity gates: the raw
// change is measured against max(1, max|field|) in the field's own unit,
// so a change is never hidden by a small scale and never exaggerated by a
// large one.
inline double buoyancy_stationarity_normalized(
    double max_change, double max_magnitude)
{
    const double scale = std::max(1.0, max_magnitude);
    return max_change / scale;
}

// Convergence-at-entry support (contract above): at the coupled fixed
// point a warm carrier re-solve can stall with EVERY residual at the
// machine-precision floor - its relative criterion can no longer be
// re-proven against the collapsed reference of the same solve, but
// the equation it is asked to satisfy already IS satisfied to machine
// precision, so the solve carries no new information. Only that exact
// situation is accepted (kCarrierMachineEquilibrium*); anything above
// the floor stays a fail-closed stop.
inline bool carrier_at_machine_equilibrium(
    const IncompressibleSolveResult& flow)
{
    if (flow.history.empty())
        return false;
    const auto& last = flow.history.back();
    return std::isfinite(last.pressure_residual) &&
           std::isfinite(last.continuity_normalized) &&
           std::isfinite(last.momentum_equation_residual_rhs_relative) &&
           std::isfinite(last.corrected_flux_continuity_linf) &&
           last.pressure_residual <=
               detail_buoyancy_stage::kCarrierMachineEquilibriumResidual &&
           last.continuity_normalized <=
               detail_buoyancy_stage::kCarrierMachineEquilibriumResidual &&
           last.corrected_flux_continuity_linf <=
               detail_buoyancy_stage::kCarrierMachineEquilibriumResidual &&
           last.momentum_equation_residual_rhs_relative <=
               detail_buoyancy_stage::kCarrierMachineEquilibriumMomentumRhs;
}

inline ProductionBuoyancyOuterLoopResult run_production_buoyant_outer_loop(
    const cfdx::core::Mesh& mesh,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& base_controls,
    const cfdx::io::CaseSetup& setup,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure,
    const std::string& restart_path = {})
{
    if (!production_buoyancy_requested(setup))
        throw std::invalid_argument(
            "run_production_buoyant_outer_loop: energy_model 'buoyant' is "
            "not declared; the flow-energy buoyancy iteration only runs "
            "on a declared buoyant case");
    const ProductionBuoyancyCoupling coupling =
        make_production_buoyancy_coupling(setup);
    const bool s2s_requested = production_s2s_radiation_requested(setup);
    const FvGeometry geometry = build_fv_geometry(mesh);
    const std::size_t n_cells = mesh.n_cells();

    // The temperature seen by the buoyancy source. Seeded with the
    // declared initial temperature, else the case reference temperature
    // (the stage-1 convention of the frozen-flow path; at the reference
    // temperature the buoyancy source is identically zero, a consistent
    // cold start).
    double initial_temperature = setup.ref_temperature;
    if (setup.has_initial_condition &&
        setup.initial_condition.temperature != 0.0)
        initial_temperature = setup.initial_condition.temperature;
    if (!std::isfinite(initial_temperature) ||
        !(initial_temperature > 0.0))
        throw std::invalid_argument(
            "production buoyancy stage: initial temperature must be finite "
            "and positive (kelvin)");
    cfdx::core::Field<double, cfdx::core::Location::CELL> relaxed_temperature(
        n_cells, "temperature", "K", 1);
    for (std::size_t cell = 0; cell < n_cells; ++cell)
        relaxed_temperature(cell) = initial_temperature;

    // Previous UNRELAXED responses for the stationarity gates. Empty
    // until the first response exists.
    std::vector<double> previous_temperature;
    std::vector<double> previous_velocity;
    std::vector<double> previous_pressure;
    bool have_previous_response = false;

    ProductionBuoyancyOuterLoopResult result;
    result.coupling = coupling;
    const double omega = detail_buoyancy_stage::kOuterRelaxation;
    const std::size_t max_outer_iterations =
        detail_buoyancy_stage::kOuterMaxIterations;
    const double tolerance =
        detail_buoyancy_stage::kOuterFixedPointTolerance;

    for (std::size_t iteration = 1; iteration <= max_outer_iterations;
         ++iteration) {
        ProductionBuoyancyOuterIteration record;
        record.iteration = iteration;
        result.outer_iterations = iteration;

        // Carrier solve with the buoyancy source of the relaxed field.
        // Warm start: the first iteration keeps the driver-declared
        // initialization strategy (including a requested restart), every
        // subsequent iteration reuses the previous outer iteration's
        // accepted state unchanged.
        // Entry snapshot for convergence-at-entry (contract above): a
        // bit-exact copy of the state the carrier solve starts from, so
        // an accepted-at-entry solve can restore it unchanged.
        auto entry_velocity = velocity;
        auto entry_pressure = pressure;
        IncompressibleSolverControls controls = base_controls;
        if (iteration > 1)
            controls.initialization.mode = InitializationMode::Provided;
        build_cell_buoyancy_acceleration(
            controls.cell_body_force, relaxed_temperature, coupling);
        const auto flow = solve_steady_incompressible(
            mesh, velocity, pressure, velocity_bcs, pressure_bcs, controls,
            restart_path);
        // Convergence-at-entry (contract above): a warm re-solve that
        // ends STAGNATED with every residual at the machine-precision
        // floor is accepted at the entry state - the entry fields are
        // restored and the previously accepted carrier (authoritative
        // mass flux) reused; the outer loop owns the coupled
        // fixed-point criteria. Anything else that is not converged
        // stays a fail-closed stop.
        const bool converged_at_entry =
            !flow.converged && iteration > 1 &&
            flow.convergence_status ==
                cfdx::core::ConvergenceStatus::STAGNATED &&
            carrier_at_machine_equilibrium(flow);
        if (flow.converged) {
            result.flow = flow;
            record.flow_converged = true;
            record.flow_iterations = flow.iterations;
            result.flow_converged_at_entry = false;
        } else if (converged_at_entry) {
            velocity = entry_velocity;
            pressure = entry_pressure;
            record.flow_converged = true;
            record.flow_converged_at_entry = true;
            record.flow_iterations = flow.iterations;
            result.flow_converged_at_entry = true;
        } else {
            // A diverged or genuinely stagnated carrier inside the
            // loop is a fail-closed stop: a flow state that is not
            // evidence cannot be iterated on.
            result.flow = flow;
            result.history.push_back(record);
            return result;
        }
        // The energy response always runs on the last FULLY proven
        // carrier: the accepted solve, or the previously accepted
        // carrier when the entry state was restored.
        const auto& carrier = result.flow;

        // Energy response on the new authoritative mass flux. The guess
        // is a copy of the relaxed field, so the response is a
        // deterministic function of the field the source was built from.
        // With declared S2S radiation the stage-3 energy-radiation fixed
        // point is the response operator (fiche 11, unchanged).
        auto temperature_guess = relaxed_temperature;
        if (s2s_requested) {
            const auto s2s_result = run_production_s2s_energy_outer_loop(
                mesh, geometry, setup, carrier.authoritative_mass_flux,
                temperature_guess);
            result.s2s = s2s_result;
            result.s2s_ran = true;
            record.s2s_ran = true;
            record.s2s_converged = s2s_result.converged;
            record.s2s_outer_iterations = s2s_result.outer_iterations;
            result.energy = s2s_result.energy;
            record.energy_converged = s2s_result.energy.converged;
            record.energy_iterations = s2s_result.energy.iterations;
        } else {
            const auto energy = run_production_energy_stage(
                mesh, geometry, setup, carrier.authoritative_mass_flux,
                temperature_guess);
            result.energy = energy;
            record.energy_converged = energy.converged;
            record.energy_iterations = energy.iterations;
        }
        result.temperature = temperature_guess;
        result.energy_ran = true;
        if (!record.energy_converged ||
            (s2s_requested && !record.s2s_converged)) {
            // A red inner stage stops the loop fail-closed; no partial
            // state is ever reported as a fixed point.
            result.history.push_back(record);
            return result;
        }

        // Stationarity of the UNRELAXED responses. The first iteration
        // has no previous response and can never prove the fixed point.
        if (have_previous_response) {
            double max_temperature_change = 0.0;
            double max_temperature_magnitude = 0.0;
            double max_velocity_change = 0.0;
            double max_velocity_magnitude = 0.0;
            double max_pressure_change = 0.0;
            double max_pressure_magnitude = 0.0;
            for (std::size_t cell = 0; cell < n_cells; ++cell) {
                const double temperature_value =
                    temperature_guess(cell);
                max_temperature_change = std::max(
                    max_temperature_change,
                    std::abs(temperature_value -
                             previous_temperature[cell]));
                max_temperature_magnitude = std::max(
                    max_temperature_magnitude,
                    std::abs(temperature_value));
                for (std::size_t component = 0; component < 3;
                     ++component) {
                    const double velocity_value =
                        velocity.component_data(component)[cell];
                    max_velocity_change = std::max(
                        max_velocity_change,
                        std::abs(velocity_value -
                                 previous_velocity[3 * cell +
                                                    component]));
                    max_velocity_magnitude = std::max(
                        max_velocity_magnitude,
                        std::abs(velocity_value));
                }
                const double pressure_value = pressure(cell);
                max_pressure_change = std::max(
                    max_pressure_change,
                    std::abs(pressure_value -
                             previous_pressure[cell]));
                max_pressure_magnitude = std::max(
                    max_pressure_magnitude, std::abs(pressure_value));
            }
            record.max_temperature_change_normalized =
                buoyancy_stationarity_normalized(
                    max_temperature_change, max_temperature_magnitude);
            record.max_velocity_change_normalized =
                buoyancy_stationarity_normalized(
                    max_velocity_change, max_velocity_magnitude);
            record.max_pressure_change_normalized =
                buoyancy_stationarity_normalized(
                    max_pressure_change, max_pressure_magnitude);
            result.max_temperature_change_normalized =
                record.max_temperature_change_normalized;
            result.max_velocity_change_normalized =
                record.max_velocity_change_normalized;
            result.max_pressure_change_normalized =
                record.max_pressure_change_normalized;
            result.fixed_point_temperature_passed =
                std::isfinite(
                    record.max_temperature_change_normalized) &&
                record.max_temperature_change_normalized <= tolerance;
            result.fixed_point_velocity_passed =
                std::isfinite(
                    record.max_velocity_change_normalized) &&
                record.max_velocity_change_normalized <= tolerance;
            result.fixed_point_pressure_passed =
                std::isfinite(
                    record.max_pressure_change_normalized) &&
                record.max_pressure_change_normalized <= tolerance;
            result.history.push_back(record);
            if (result.fixed_point_temperature_passed &&
                result.fixed_point_velocity_passed &&
                result.fixed_point_pressure_passed) {
                result.converged = true;  // exact coupled fixed point
                return result;
            }
        } else {
            result.history.push_back(record);
        }

        // Under-relaxed feedback for the next carrier solve
        // (path only): T_relaxed(k+1) = T_relaxed(k)
        //                                + omega (T(k) - T_relaxed(k)).
        for (std::size_t cell = 0; cell < n_cells; ++cell)
            relaxed_temperature(cell) +=
                omega * (temperature_guess(cell) -
                         relaxed_temperature(cell));

        previous_temperature.assign(n_cells, 0.0);
        previous_velocity.assign(3 * n_cells, 0.0);
        previous_pressure.assign(n_cells, 0.0);
        for (std::size_t cell = 0; cell < n_cells; ++cell) {
            previous_temperature[cell] = temperature_guess(cell);
            for (std::size_t component = 0; component < 3; ++component)
                previous_velocity[3 * cell + component] =
                    velocity.component_data(component)[cell];
            previous_pressure[cell] = pressure(cell);
        }
        have_previous_response = true;
    }
    // Budget exhausted: converged=false, full history, no claim.
    return result;
}

}  // namespace cfdx::physics
