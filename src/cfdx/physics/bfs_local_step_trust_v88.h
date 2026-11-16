#pragma once
// V88 opt-in pre-SST trust-region check of the *increment* of U, rather
// than only the magnitude of U. A failed trial is recomputed by the existing
// SIMPLE pressure/velocity solve after full nonlinear rollback; this helper
// neither clips velocities nor alters conservative Rhie-Chow face fluxes.
// This is a STARTUP/progression safeguard, never a physical convergence gate.
#include <cmath>
#include <stdexcept>
#include <string>
namespace cfdx::physics {
struct BfsLocalStepTrialV88 {
    bool may_continue_sst = false;
    double max_delta_velocity_ratio = 0.0;
    double pressure_increment_span_over_qin = 0.0; // Diagnostic only; gauge-invariant.
    std::string reason;
};
inline BfsLocalStepTrialV88 check_bfs_local_step_v88(
    double maximum_cell_velocity_increment,
    double pressure_increment_span,
    double inlet_speed,
    double density,
    double allowed_delta_u_over_uin = 0.5)
{
    if (!(std::isfinite(maximum_cell_velocity_increment) &&
          maximum_cell_velocity_increment >= 0.0 &&
          std::isfinite(pressure_increment_span) && pressure_increment_span >= 0.0 &&
          std::isfinite(inlet_speed) && inlet_speed > 0.0 &&
          std::isfinite(density) && density > 0.0 &&
          std::isfinite(allowed_delta_u_over_uin) &&
          allowed_delta_u_over_uin > 0.0 && allowed_delta_u_over_uin <= 2.0))
        throw std::invalid_argument("V88 invalid or nonfinite local-step trial");
    BfsLocalStepTrialV88 trial;
    trial.max_delta_velocity_ratio = maximum_cell_velocity_increment / inlet_speed;
    trial.pressure_increment_span_over_qin =
        (pressure_increment_span / (0.5 * density * inlet_speed)) / inlet_speed;
    if (!std::isfinite(trial.max_delta_velocity_ratio) ||
        !std::isfinite(trial.pressure_increment_span_over_qin))
        throw std::invalid_argument("V88 local-step diagnostic overflow");
    trial.may_continue_sst = trial.max_delta_velocity_ratio <= allowed_delta_u_over_uin;
    trial.reason = trial.may_continue_sst ? "delta_velocity_preflight_only_not_certified"
                                          : "local_velocity_increment_envelope";
    return trial;
}
} // namespace cfdx::physics