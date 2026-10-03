#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
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

    const auto initial_U = U;
    const auto initial_p = p;
    const auto scale_velocity_bcs = [&](double parameter) {
        VelocityBoundaryConditions scaled = velocity_bcs;
        for (auto& [name, bc] : scaled) {
            (void)name;
            if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                bc.value = bc.value * parameter;
        }
        return scaled;
    };

    double parameter = 0.0;
    double step = continuation.initial_step;
    std::size_t stage_index = 0;
    std::size_t total_attempts = 0;
    IncompressibleSolveResult last_solver_result;

    while (parameter < 1.0) {
        if (++stage_index > continuation.max_stages)
            throw std::runtime_error(
                "continuation exceeded the configured maximum number of stages");

        const double target = continuation_next_target(parameter, step);
        std::size_t attempts = 0;
        bool accepted = false;

        while (!accepted) {
            if (++attempts > continuation.max_stage_attempts) {
                U = initial_U;
                p = initial_p;
                throw std::runtime_error(
                    "continuation exhausted the configured stage-attempt budget");
            }
            ++total_attempts;

            const auto previous_U = U;
            const auto previous_p = p;
            IncompressibleSolverControls stage_controls = controls;
            stage_controls.body_force = controls.body_force * target;
            const auto stage_velocity_bcs = scale_velocity_bcs(target);

            IncompressibleSolveResult solver_result;
            try {
                solver_result = solve_steady_incompressible(
                    mesh, U, p, stage_velocity_bcs, pressure_bcs,
                    stage_controls, {}, {});
            } catch (const std::runtime_error& error) {
                U = previous_U;
                p = previous_p;
                if (step <= continuation.minimum_step) {
                    throw std::runtime_error(
                        std::string("continuation stage failed at minimum step: ") +
                        error.what());
                }
                step = std::max(
                    continuation.minimum_step,
                    continuation_step_after_failure(step, continuation));
                continue;
            }

            last_solver_result = solver_result;
            const bool stage_ok =
                solver_result.converged &&
                solver_result.convergence_status ==
                    cfdx::core::ConvergenceStatus::CONVERGED;

            if (stage_ok) {
                parameter = target;
                accepted = true;
                result.stages.push_back(ContinuationStageReport{
                    stage_index, attempts, parameter, step, true,
                    solver_result.convergence_status,
                    solver_result.convergence_reason, solver_result});
                step = continuation_step_after_success(step, continuation);
            } else {
                U = previous_U;
                p = previous_p;
                result.stages.push_back(ContinuationStageReport{
                    stage_index, attempts, target, step, false,
                    solver_result.convergence_status,
                    solver_result.convergence_reason, solver_result});
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
