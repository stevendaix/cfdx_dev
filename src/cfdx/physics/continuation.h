#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

/**
 * Deterministic load-continuation controls for steady incompressible solves.
 *
 * The continuation parameter lambda scales the configured body force and
 * fixed-value velocity boundary data from 0 to 1. Pressure boundary data are
 * intentionally not scaled. Failed stages are transactional: U and p are
 * restored to the last accepted stage before the step is reduced and retried.
 * The controller never changes nonlinear or linear tolerances.
 */
struct ContinuationControls {
    bool enabled = false;
    double initial_step = 0.25;
    double minimum_step = 1e-2;
    double maximum_step = 0.5;
    double step_growth = 1.5;
    double step_reduction = 0.5;
    std::size_t max_stage_attempts = 8;
    std::size_t max_stages = 1000;
    // Experimental two-phase startup: momentum-only stages first, full
    // SST on lambda=1. An intermediate laminar solution is NEVER a final
    // turbulence qualification. Must remain an explicit opt-in.
    bool momentum_only_before_final = false;
    // Intermediate stages may use robust first order, but lambda=1 always
    // restores the user-selected final convection scheme.
    bool upwind_before_final = false;
};

inline void validate_continuation_controls(const ContinuationControls& c)
{
    if (!(c.initial_step > 0.0) ||
        !(c.minimum_step > 0.0) ||
        !(c.maximum_step >= c.initial_step) ||
        !(c.minimum_step <= c.initial_step) ||
        !(c.step_growth > 1.0) ||
        !(c.step_reduction > 0.0 && c.step_reduction < 1.0) ||
        c.max_stage_attempts == 0 ||
        c.max_stages == 0 ||
        !std::isfinite(c.initial_step) ||
        !std::isfinite(c.minimum_step) ||
        !std::isfinite(c.maximum_step) ||
        !std::isfinite(c.step_growth) ||
        !std::isfinite(c.step_reduction)) {
        throw std::invalid_argument("invalid continuation controls");
    }
}

inline bool continuation_final_physical_stage(double parameter) {
    if (!std::isfinite(parameter) || parameter < 0.0 || parameter > 1.0)
        throw std::invalid_argument("invalid continuation physical-stage parameter");
    return parameter == 1.0;
}

inline double continuation_next_target(double current, double step)
{
    if (!std::isfinite(current) || !std::isfinite(step) ||
        current < 0.0 || current > 1.0 || !(step > 0.0))
        throw std::invalid_argument("invalid continuation state");
    return std::min(1.0, current + step);
}

inline double continuation_step_after_success(
    double step, const ContinuationControls& c)
{
    validate_continuation_controls(c);
    if (!std::isfinite(step) || !(step > 0.0))
        throw std::invalid_argument("invalid continuation step");
    return std::min(c.maximum_step, step * c.step_growth);
}

inline double continuation_step_after_failure(
    double step, const ContinuationControls& c)
{
    validate_continuation_controls(c);
    if (!std::isfinite(step) || !(step > 0.0))
        throw std::invalid_argument("invalid continuation step");
    return step * c.step_reduction;
}

// Captures every persistent field affected by an SST-coupled stage.
// The solver maintains its own per-iteration retry, but a continuation
// rejection must also roll back the *whole* stage, not merely U and p.
struct ContinuationCoupledSnapshot {
    cfdx::core::Field<double, cfdx::core::Location::CELL> velocity;
    cfdx::core::Field<double, cfdx::core::Location::CELL> pressure;
    std::optional<SSTProductionState> turbulence;
    std::vector<double> mass_flux;

    ContinuationCoupledSnapshot(
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
        const SSTProductionState* sst, const std::vector<double>& flux)
        : velocity(U), pressure(p), mass_flux(flux)
    {
        if (sst) turbulence = snapshot_sst_production_state(*sst);
    }

    void restore(
        cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
        cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
        SSTProductionState* sst, std::vector<double>& flux) const
    {
        if (bool(turbulence) != (sst != nullptr))
            throw std::invalid_argument("continuation snapshot: SST runtime changed");
        U = velocity;
        p = pressure;
        if (sst) restore_sst_production_state(*sst, *turbulence);
        flux = mass_flux;
    }
};

// A previously corrected conservative flux is a valid warm start only
// when the fixed-velocity boundary contract is scaled by the SAME factor.
// This is a stage-seed transformation, never a replacement for the full
// physical flux/conservation gate at lambda=1.
inline std::vector<double> continuation_scale_mass_flux(
    const std::vector<double>& flux, double factor)
{
    if (!std::isfinite(factor) || !(factor > 0.0))
        throw std::invalid_argument("continuation: invalid face flux scale");
    std::vector<double> scaled = flux;
    for (double& f : scaled) {
        if (!std::isfinite(f))
            throw std::invalid_argument("continuation: non-finite face flux");
        f *= factor;
        if (!std::isfinite(f))
            throw std::invalid_argument("continuation: face flux overflow");
    }
    return scaled;
}

struct ContinuationStageReport {
    std::size_t stage = 0;
    std::size_t attempts = 0;
    double parameter = 0.0;
    double step = 0.0;
    bool converged = false;
    cfdx::core::ConvergenceStatus status =
        cfdx::core::ConvergenceStatus::CONTINUE;
    std::string reason;
    IncompressibleSolveResult solver_result;
    // This is the only stage at which SST and user-selected order
    // are necessarily restored. Intermediate stage convergence is NEVER
    // a certificate for the final physical problem.
    bool full_physics = false;
};

struct ContinuationSolveResult {
    bool converged = false;
    double final_parameter = 0.0;
    std::size_t total_attempts = 0;
    std::string reason;
    IncompressibleSolveResult final_solver_result;
    std::vector<ContinuationStageReport> stages;
};

inline ContinuationSolveResult solve_steady_incompressible_continuation(
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& controls,
    const ContinuationControls& continuation = {},
    const std::string& restart_path = {},
    cfdx::io::DatRestartFields restart_fields = {})
{
    validate_continuation_controls(continuation);

    ContinuationSolveResult result;
    if (!continuation.enabled) {
        result.final_solver_result = solve_steady_incompressible(
            mesh, U, p, velocity_bcs, pressure_bcs, controls,
            restart_path, restart_fields);
        result.converged = result.final_solver_result.converged;
        result.final_parameter = result.converged ? 1.0 : 0.0;
        result.total_attempts = 1;
        result.reason = result.final_solver_result.convergence_reason;
        return result;
    }

    const auto scale_velocity_bcs = [&](double parameter) {
        VelocityBoundaryConditions scaled = velocity_bcs;
        for (auto& [name, bc] : scaled) {
            (void)name;
            if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                bc.value = cfdx::core::Vec3{
                    bc.value.x * parameter,
                    bc.value.y * parameter,
                    bc.value.z * parameter};
        }
        return scaled;
    };

    double parameter = 0.0;
    double step = continuation.initial_step;
    std::size_t stage_index = 0;
    std::size_t total_attempts = 0;
    auto* sst_state = controls.sst_runtime
        ? controls.sst_runtime->state : nullptr;
    // In particular, do not read a DAT restart again after the first stage.
    std::vector<double> accepted_flux = controls.initial_face_mass_flux;
    ContinuationCoupledSnapshot last_accepted(
        U, p, sst_state, accepted_flux);
    IncompressibleSolveResult last_solver_result;

    while (parameter < 1.0) {
        if (++stage_index > continuation.max_stages)
            throw std::runtime_error(
                "continuation exceeded the configured maximum number of stages");

        std::size_t attempts = 0;
        bool accepted = false;

        while (!accepted) {
            const double target = continuation_next_target(parameter, step);
            if (++attempts > continuation.max_stage_attempts) {
                last_accepted.restore(U, p, sst_state, accepted_flux);
                result.final_parameter = parameter;
                result.total_attempts = total_attempts;
                result.final_solver_result = last_solver_result;
                throw std::runtime_error(
                    "continuation exhausted the configured stage-attempt budget");
            }
            ++total_attempts;

            const ContinuationCoupledSnapshot previous(
                U, p, sst_state, accepted_flux);
            IncompressibleSolverControls stage_controls = controls;
            stage_controls.body_force = cfdx::core::Vec3{
                controls.body_force.x * target,
                controls.body_force.y * target,
                controls.body_force.z * target};
            const bool final_physical = continuation_final_physical_stage(target);
            if (!final_physical && continuation.momentum_only_before_final) {
                stage_controls.sst_runtime = nullptr;
                stage_controls.turbulent_viscosity = 0.0;
                stage_controls.cell_turbulent_viscosity.clear();
            }
            if (!final_physical && continuation.upwind_before_final)
                stage_controls.convection_scheme = ConvectionScheme::UPWIND;
            const auto stage_velocity_bcs = scale_velocity_bcs(target);
            // Keep the conservative face seed consistent with each stage's
            // fixed-velocity BCs. When a stage succeeds, its authoritative
            // corrected flux becomes the seed for the next stage.
            if (!accepted_flux.empty()) {
                stage_controls.initial_face_mass_flux = continuation_scale_mass_flux(
                    accepted_flux, parameter > 0.0 ? target / parameter : target);
            }

            // Stage 0 initializes the ramp from the declared strategy. Every
            // later attempt inherits the last accepted state (Provided), so a
            // restart is re-read only on the first attempt of the first stage
            // and cannot silently overwrite a converged stage.
            const bool first_attempt_of_first_stage =
                (stage_index == 1 && attempts == 1);
            if (!first_attempt_of_first_stage) {
                stage_controls.initialization.mode = InitializationMode::Provided;
                stage_controls.initialization.restart_path.clear();
            }

            IncompressibleSolveResult solver_result;
            try {
                const auto stage_restart_path =
                    first_attempt_of_first_stage ? restart_path : std::string{};
                const auto stage_restart_fields =
                    first_attempt_of_first_stage
                        ? restart_fields : cfdx::io::DatRestartFields{};
                solver_result = solve_steady_incompressible(
                    mesh, U, p, stage_velocity_bcs, pressure_bcs,
                    stage_controls, stage_restart_path, stage_restart_fields);
            } catch (const std::runtime_error& error) {
                previous.restore(U, p, sst_state, accepted_flux);
                if (step <= continuation.minimum_step) {
                    throw std::runtime_error(
                        std::string("continuation stage failed at minimum step: ") +
                        error.what());
                }
                step = std::max(
                    continuation.minimum_step,
                    continuation_step_after_failure(step, continuation));
                continue;
            } catch (...) {
                previous.restore(U, p, sst_state, accepted_flux);
                throw;
            }

            last_solver_result = solver_result;
            const bool stage_ok =
                solver_result.converged &&
                solver_result.convergence_status ==
                    cfdx::core::ConvergenceStatus::CONVERGED;

            if (stage_ok) {
                parameter = target;
                accepted_flux.clear();
                if (solver_result.authoritative_mass_flux.size() == mesh.n_faces()) {
                    accepted_flux.resize(mesh.n_faces());
                    for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                        accepted_flux[f] = solver_result.authoritative_mass_flux(f);
                }
                last_accepted = ContinuationCoupledSnapshot(
                    U, p, sst_state, accepted_flux);
                accepted = true;
                result.stages.push_back(ContinuationStageReport{
                    stage_index, attempts, parameter, step, true,
                    solver_result.convergence_status,
                    solver_result.convergence_reason, solver_result});
                result.stages.back().full_physics = final_physical;
                step = continuation_step_after_success(step, continuation);
            } else {
                previous.restore(U, p, sst_state, accepted_flux);
                result.stages.push_back(ContinuationStageReport{
                    stage_index, attempts, target, step, false,
                    solver_result.convergence_status,
                    solver_result.convergence_reason, solver_result});
                result.stages.back().full_physics = final_physical;
                if (step <= continuation.minimum_step) {
                    result.final_solver_result = solver_result;
                    result.final_parameter = parameter;
                    result.total_attempts = total_attempts;
                    result.reason =
                        "continuation stage failed at minimum step";
                    return result;
                }
                step = std::max(
                    continuation.minimum_step,
                    continuation_step_after_failure(step, continuation));
            }
        }
    }

    result.converged = true;
    result.final_parameter = parameter;
    result.total_attempts = total_attempts;
    result.final_solver_result = last_solver_result;
    result.reason = "all continuation stages converged";
    return result;
}

} // namespace cfdx::physics