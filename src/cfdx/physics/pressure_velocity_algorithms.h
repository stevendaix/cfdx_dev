#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/flux.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

enum class PressureVelocityAlgorithm {
    SIMPLE,
    SIMPLEC,
    PISO,
    PIMPLE,
    // Fractional-step projection: momentum predictor followed by one or
    // more pressure projections. This is intentionally distinct from PISO:
    // each projection operates on the current conservative face flux.
    FRACTIONAL_STEP,
    // Fully coupled pressure-based solve. Momentum and continuity are
    // assembled in one block system with an implicit Rhie-Chow pressure block.
    COUPLED
};

enum class CoupledSchurModel { BlockLocal, PCD, LSC, BFBT, SIMPLE, SIMPLEC };

inline const char* to_string(CoupledSchurModel model) noexcept {
    switch (model) {
        case CoupledSchurModel::BlockLocal: return "block_local";
        case CoupledSchurModel::PCD: return "pcd";
        case CoupledSchurModel::LSC: return "lsc";
        case CoupledSchurModel::BFBT: return "bfbt";
        case CoupledSchurModel::SIMPLE: return "simple";
        case CoupledSchurModel::SIMPLEC: return "simplec";
    }
    return "unknown";
}

inline bool parse_coupled_schur_model(const std::string& name, CoupledSchurModel& out) {
    if (name == "block_local") { out = CoupledSchurModel::BlockLocal; return true; }
    if (name == "pcd") { out = CoupledSchurModel::PCD; return true; }
    if (name == "lsc") { out = CoupledSchurModel::LSC; return true; }
    if (name == "bfbt") { out = CoupledSchurModel::BFBT; return true; }
    if (name == "simple") { out = CoupledSchurModel::SIMPLE; return true; }
    if (name == "simplec") { out = CoupledSchurModel::SIMPLEC; return true; }
    return false;
}

struct CouplingControls {
    double alpha_u = 0.7;
    double alpha_p = 0.3;
    int n_outer_correctors = 1;
    int n_pressure_correctors = 2;
    // EXPERIMENTAL, opt-in: for PISO/PIMPLE, rebuild H/A from the velocity
    // corrected by the previous pressure solve, using the FROZEN relaxed
    // momentum matrix. This incorporates neighbour velocity corrections
    // neglected by a frozen-predictor pressure-only loop (Issa, 1986).
    // Physical and temporal-order qualification is required before promotion.
    bool piso_neighbor_correction = false;
    int n_fractional_steps = 2;
    // Coupled linear solve controls. The coupled matrix is indefinite but the
    // Rhie-Chow pressure block makes the collocated formulation nonsingular
    // after pressure gauge fixing.
    std::size_t coupled_max_iterations = 2000;
    // Default tightened from 1e-10. At 1e-10 the coupled FGMRES true
    // relative residual floors at ~9.94e-11 with an absolute momentum
    // residual floor of ~2.4e-10 (stable within 0.02% over 8000 outer
    // iterations on the NACA0012 n64 case, diag/naca-coupled-fix4b),
    // which leaves the iteration-1-relative momentum metric ~11x above
    // the 1e-8 gate: the previous default could never pass the
    // convergence gates. At 1e-12 the same case converges with all
    // gates green (n64 and n128 both in 67 outer iterations,
    // diag/naca-coupled-fix4c). See issue #663.
    double coupled_linear_tolerance = 1e-12;
    // Explicit pressure-Schur approximation for the coupled production path.
    // Unsupported models must be rejected; there is no implicit fallback.
    CoupledSchurModel schur_model = CoupledSchurModel::PCD;
    // Research-only, opt-in. Uses -Mp^-1 Fp Kp^-1 with a positive, lumped
    // pressure mass. The default remains the qualified legacy -Fp^-1 action.
    // No physical performance claim or automatic fallback is implied.
    bool pcd_lumped_mass_commutator_experimental = false;
    // Diagnostic A/B only: align the nested pressure operators' gauge column
    // with the monolithic eliminated gauge column. No production default change.
    bool pcd_consistent_gauge_experimental = false;
    // Research diagnostics only: exact row scaling of Fp and guarded nested
    // actions. These change only the Schur preconditioner, not FV equations.
    bool pcd_fp_row_equilibration_experimental = false;
    bool pcd_fp_quality_guard_experimental = false;
    bool pcd_fp_second_pass_experimental = false;
    // Experimental Fp-only strong-component anchor. The coupled FV matrix,
    // physical pressure gauge and nonlinear gates are unchanged.
    bool pcd_fp_recirculation_anchor_experimental = false;
    // Opt-in certified nonsymmetric Fp deflation. Must carry a LEFT/RIGHT
    // basis file whose coefficient fingerprint matches the CURRENT Fp.
    // The matrix-free solve changes ONLY the PCD action, never the coupled FV
    // system. No automatic basis reuse when numerical coefficients change.
    std::string pcd_fp_bisided_basis_file_experimental;
    // V15: select a CURRENT, independently certified basis from a directory
    // keyed by the exact Fp CSR coefficient fingerprint.  No stale fallback.
    std::string pcd_fp_bisided_basis_directory_v15_experimental;
    // V15: invert Fp + W Z^T on the FULL auxiliary right-hand side.  Unlike
    // the V14 projected action, this preserves the deficient Schur directions
    // so that the preconditioner need not map left-null vectors to zero.
    bool pcd_fp_full_rank_action_v15_experimental = false;
    // V16: if a converged deflated action narrowly fails the ORIGINAL
    // projected residual certificate, refine the same operator at 10x and
    // 100x tighter inner tolerances. No relaxed acceptance or Schur fallback.
    bool pcd_fp_refine_projected_v16_experimental = false;
    // Opt-in, SMALL coupled matrices only (<=4096 unknowns): retry a failed
    // Schur/PCD outer solve with a bounded-fill monolithic ILUT.  The
    // original saddle matrix, physical residual and tolerances are unchanged;
    // acceptance requires the true original-system residual independently.
    // NOT a parallel/GPU preconditioner, and never a silent fallback.
    bool coupled_monolithic_ilut_fallback_v16_experimental = false;

    // Preconditioner-only A/B: pressure-space Dirichlet data on actual
    // inflow faces, including mixed-in/out farfield patches. The physical
    // momentum/continuity BCs and the Schur matrix A are NOT modified.
    // Off by default until NACA/N9/BFS full physical A/B qualification.
    bool pcd_pressure_inflow_dirichlet_experimental = false;
    // Shared diagnostic across BlockLocal, PCD, LSC, BFBt, SIMPLE and
    // SIMPLEC Schur. Off by default; no silent substitution on rejection.
    bool schur_action_growth_guard_experimental = false;
    // Diagnostic A/B only: align the nested pressure operators' gauge column
    // with the monolithic eliminated gauge column. No production default change.
    bool pcd_consistent_gauge_experimental = false;
    // Research diagnostics only: exact row scaling of Fp and guarded nested
    // actions. These change only the Schur preconditioner, not FV equations.
    bool pcd_fp_row_equilibration_experimental = false;
    bool pcd_fp_quality_guard_experimental = false;
    bool pcd_fp_second_pass_experimental = false;
    // Experimental Fp-only strong-component anchor. The coupled FV matrix,
    // physical pressure gauge and nonlinear gates are unchanged.
    bool pcd_fp_recirculation_anchor_experimental = false;
    // Opt-in certified nonsymmetric Fp deflation. Must carry a LEFT/RIGHT
    // basis file whose coefficient fingerprint matches the CURRENT Fp.
    // The matrix-free solve changes ONLY the PCD action, never the coupled FV
    // system. No automatic basis reuse when numerical coefficients change.
    std::string pcd_fp_bisided_basis_file_experimental;
    // V15: select a CURRENT, independently certified basis from a directory
    // keyed by the exact Fp CSR coefficient fingerprint.  No stale fallback.
    std::string pcd_fp_bisided_basis_directory_v15_experimental;
    // V15: invert Fp + W Z^T on the FULL auxiliary right-hand side.  Unlike
    // the V14 projected action, this preserves the deficient Schur directions
    // so that the preconditioner need not map left-null vectors to zero.
    bool pcd_fp_full_rank_action_v15_experimental = false;

    // Preconditioner-only A/B: pressure-space Dirichlet data on actual
    // inflow faces, including mixed-in/out farfield patches. The physical
    // momentum/continuity BCs and the Schur matrix A are NOT modified.
    // Off by default until NACA/N9/BFS full physical A/B qualification.
    bool pcd_pressure_inflow_dirichlet_experimental = false;
    // Shared diagnostic across BlockLocal, PCD, LSC, BFBt, SIMPLE and
    // SIMPLEC Schur. Off by default; no silent substitution on rejection.
    bool schur_action_growth_guard_experimental = false;
    // Serial, small-subdomain, experimental RAS/weighted Schwarz for COUPLED
    // with BlockLocal Schur. Explicit opt-in; setup failure is fatal rather
    // than silently replaced by AMG. No MPI / coarse-space scaling claim.
    bool coupled_schwarz_experimental = false;
    std::size_t coupled_schwarz_cells_per_subdomain = 1;
    std::size_t coupled_schwarz_overlap_hops = 1;
    std::size_t coupled_schwarz_max_local_unknowns = 64;
    bool coupled_schwarz_weighted = false;
};

inline void validate_coupling_controls(const CouplingControls& c)
{
    if (!(c.alpha_u > 0.0 && c.alpha_u <= 1.0) ||
        !(c.alpha_p > 0.0 && c.alpha_p <= 1.0) ||
        c.n_outer_correctors < 1 || c.n_pressure_correctors < 1 ||
        c.n_fractional_steps < 1 || c.coupled_max_iterations == 0 ||
        !(c.coupled_linear_tolerance > 0.0))
        throw std::invalid_argument("pressure-velocity controls: invalid relaxation/corrector count");
    if (c.coupled_schwarz_experimental &&
        (c.coupled_schwarz_cells_per_subdomain == 0 ||
         c.coupled_schwarz_overlap_hops > 8 ||
         c.coupled_schwarz_max_local_unknowns < 4 ||
         c.coupled_schwarz_max_local_unknowns > 256 ||
         c.coupled_schwarz_max_local_unknowns % 4 != 0 ||
         c.coupled_schwarz_cells_per_subdomain > c.coupled_schwarz_max_local_unknowns / 4))
        throw std::invalid_argument("COUPLED Schwarz: invalid subdomain/overlap controls");
}

struct AdaptiveRelaxationControls {
    bool enabled = false;
    double min_alpha_u = 0.2;
    double max_alpha_u = 0.9;
    double min_alpha_p = 0.1;
    double max_alpha_p = 0.5;
    double improvement_threshold = 0.05;
    double degradation_threshold = 0.10;
    // Small additive increases make the controller conservative near a stable
    // fixed point; degradation remains twice as fast so overshoot is recovered
    // more quickly than it is introduced.
    double increase_step = 0.005;
    double decrease_step = 0.02;
    // A severe residual jump immediately following an accepted relaxation
    // increase is treated as an overshoot event rather than waiting for the
    // normal multi-window hysteresis. Residual jumps without a causal alpha
    // increase are not classified as adaptive overshoots.
    double severe_degradation_ratio = 2.0;
    // Number of adaptation windows during which increases are frozen after a
    // severe-degradation rollback.
    std::size_t recovery_cooldown_windows = 2;
    // Adapt against a short history window rather than consecutive noisy
    // iterations. This prevents the controller from reacting to the
    // pressure/velocity coupling oscillation of SIMPLE itself.
    std::size_t adaptation_window = 4;
    // Number of consecutive window decisions required before changing a
    // relaxation factor. This is deliberately explicit: the controller is
    // hysteretic and must not react to a single SIMPLE oscillation.
    std::size_t required_trend_windows = 2;
};

inline void validate_adaptive_relaxation_controls(const AdaptiveRelaxationControls& c)
{
    if (c.min_alpha_u <= 0.0 || c.max_alpha_u < c.min_alpha_u || c.max_alpha_u > 1.0 ||
        c.min_alpha_p <= 0.0 || c.max_alpha_p < c.min_alpha_p || c.max_alpha_p > 1.0 ||
        c.improvement_threshold < 0.0 || c.degradation_threshold < c.improvement_threshold ||
        !(c.increase_step > 0.0) || !(c.decrease_step > 0.0) ||
        c.increase_step > (c.max_alpha_u - c.min_alpha_u) ||
        c.decrease_step > (c.max_alpha_u - c.min_alpha_u) ||
        c.increase_step > (c.max_alpha_p - c.min_alpha_p) ||
        c.decrease_step > (c.max_alpha_p - c.min_alpha_p) ||
        !(c.severe_degradation_ratio > 1.0 + c.degradation_threshold) ||
        c.recovery_cooldown_windows == 0 ||
        c.adaptation_window == 0 || c.required_trend_windows == 0)
        throw std::invalid_argument("invalid adaptive relaxation controls");
}

// V21: a full-flux V19 pressure correction drives continuity close to
// round-off before U is consistent. Its tiny continuity residual is NOT
// a useful merit for adapting alpha_p. Use the reassembled momentum defect
// and the independently reconstructed face-flux mismatch instead. These are
// diagnostics of coupling, not substitutes for the final physical gates.
inline double v21_pressure_relaxation_merit(
    bool full_flux_step, double continuity, double momentum_rhs_relative,
    double flux_velocity_mismatch_normalized)
{
    if (!std::isfinite(continuity) || continuity < 0.0)
        throw std::invalid_argument("V21 adaptive pressure merit: invalid continuity");
    // Preserve the historical criterion (and its required inputs) when V19
    // split flux is not active. An unrelated optional diagnostic must never
    // invalidate a previously supported SIMPLE/PISO relaxation path.
    if (!full_flux_step)
        return continuity;
    if (!std::isfinite(momentum_rhs_relative) || momentum_rhs_relative < 0.0 ||
        !std::isfinite(flux_velocity_mismatch_normalized) ||
        flux_velocity_mismatch_normalized < 0.0)
        throw std::invalid_argument("V21 adaptive pressure merit: invalid coupling residual");
    return std::max(momentum_rhs_relative, flux_velocity_mismatch_normalized);
}

inline double adapt_relaxation_factor(
    double alpha,
    double previous_residual,
    double residual,
    double min_alpha,
    double max_alpha,
    const AdaptiveRelaxationControls& controls)
{
    validate_adaptive_relaxation_controls(controls);
    if (!std::isfinite(alpha) || !std::isfinite(previous_residual) ||
        !std::isfinite(residual) || previous_residual < 0.0 || residual < 0.0)
        throw std::invalid_argument("invalid adaptive relaxation state");
    const double safe_alpha = std::clamp(alpha, min_alpha, max_alpha);
    if (!controls.enabled || previous_residual == 0.0)
        return safe_alpha;
    const double ratio = residual / previous_residual;
    if (ratio <= 1.0 - controls.improvement_threshold)
        return std::min(max_alpha, safe_alpha + controls.increase_step);
    if (ratio >= 1.0 + controls.degradation_threshold)
        return std::max(min_alpha, safe_alpha - controls.decrease_step);
    // Hysteresis band: retain the current relaxation and filter small residual noise.
    return safe_alpha;
}

inline double adapt_relaxation_factor_windowed(
    double alpha,
    double reference_metric,
    double current_metric,
    double min_alpha,
    double max_alpha,
    const AdaptiveRelaxationControls& controls,
    std::size_t& improvement_streak,
    std::size_t& degradation_streak)
{
    validate_adaptive_relaxation_controls(controls);
    if (!std::isfinite(alpha) || !std::isfinite(reference_metric) ||
        !std::isfinite(current_metric) || reference_metric <= 0.0 ||
        current_metric < 0.0)
        throw std::invalid_argument("invalid windowed adaptive relaxation state");

    const double safe_alpha = std::clamp(alpha, min_alpha, max_alpha);
    if (!controls.enabled)
        return safe_alpha;

    const double ratio = current_metric / reference_metric;
    if (ratio <= 1.0 - controls.improvement_threshold) {
        ++improvement_streak;
        degradation_streak = 0;
    } else if (ratio >= 1.0 + controls.degradation_threshold) {
        ++degradation_streak;
        improvement_streak = 0;
    } else {
        improvement_streak = 0;
        degradation_streak = 0;
    }

    if (improvement_streak >= controls.required_trend_windows) {
        improvement_streak = 0;
        return std::min(max_alpha, safe_alpha + controls.increase_step);
    }
    if (degradation_streak >= controls.required_trend_windows) {
        degradation_streak = 0;
        return std::max(min_alpha, safe_alpha - controls.decrease_step);
    }
    return safe_alpha;
}

inline bool adaptive_relaxation_severe_degradation(
    double previous_metric,
    double current_metric,
    double previous_alpha,
    double current_alpha,
    const AdaptiveRelaxationControls& controls)
{
    validate_adaptive_relaxation_controls(controls);
    if (!std::isfinite(previous_metric) || !std::isfinite(current_metric) ||
        !std::isfinite(previous_alpha) || !std::isfinite(current_alpha) ||
        previous_metric <= 0.0 || current_metric < 0.0 ||
        previous_alpha <= 0.0 || current_alpha <= 0.0)
        throw std::invalid_argument("invalid adaptive relaxation guard state");

    // A residual jump is only an adaptive-relaxation overshoot when the
    // responsible relaxation factor actually increased between two accepted
    // nonlinear states. This prevents ordinary SIMPLE oscillations from
    // repeatedly consuming the bounded retry budget.
    if (!(current_alpha > previous_alpha))
        return false;
    return current_metric / previous_metric >= controls.severe_degradation_ratio;
}

inline double relaxed_value(double old_value, double computed_value, double alpha)
{
    if (!std::isfinite(alpha) || !std::isfinite(old_value) || !std::isfinite(computed_value) ||
        alpha <= 0.0 || alpha > 1.0)
        throw std::invalid_argument("relaxed_value: alpha must be in (0,1]");
    return old_value + alpha * (computed_value - old_value);
}

// SIMPLEC uses the signed entries stored in the assembled momentum matrix.
// CFDX stores neighbour entries as A_PN = -a_PN with a_PN >= 0, so the
// consistent diagonal A_P - sum(a_PN) is A_P + sum(A_PN). Keeping this
// sign conversion explicit prevents the segregated and algebraic Schur paths
// from drifting away from the assembled matrix convention.
inline double simplec_consistent_diagonal(
    double diagonal,
    double signed_offdiag_sum)
{
    if (!std::isfinite(diagonal) || !std::isfinite(signed_offdiag_sum) ||
        diagonal <= 0.0)
        throw std::invalid_argument(
            "simplec_consistent_diagonal: invalid momentum coefficients");
    // Matrix off-diagonals are signed (normally non-positive): convert
    // them back to the positive neighbour-coefficient convention of SIMPLEC.
    const double denominator = diagonal + signed_offdiag_sum;
    if (!(denominator > 0.0) || !std::isfinite(denominator))
        throw std::invalid_argument(
            "simplec_consistent_diagonal: non-positive/non-finite denominator");
    return denominator;
}

// Absolute-value form is equivalent only when all stored off-diagonals
// are non-positive. Kept as an independent diagnostic, not as the contract.
inline double simplec_abs_sum_diagnostic_denominator(
    double diagonal,
    double abs_offdiag_sum)
{
    if (!std::isfinite(diagonal) || !std::isfinite(abs_offdiag_sum) ||
        diagonal <= 0.0 || abs_offdiag_sum < 0.0)
        throw std::invalid_argument(
            "simplec_abs_sum_diagnostic_denominator: invalid coefficients");
    return diagonal - abs_offdiag_sum;
}

// Discrete coefficient convention used throughout the pressure-velocity path:
//
//   A_P   = integrated FV momentum diagonal
//   rAU   = 1/A_P
//   dAU   = V*rAU = V/A_P
//
// The volume is intentionally NOT folded into rAU. rAU belongs to the
// momentum algebra; dAU is the physical velocity response to a pressure
// gradient. Pressure correction, Rhie-Chow, and coupled Schur coefficients
// must use dAU (or an equivalent V*rAU expression) so they remain consistent
// with the cell velocity reconstruction and continuity flux. Keep this
// distinction explicit because 1/A_P and V/A_P are both valid quantities,
// but they belong to different discrete operators.

inline double piso_correction_gain(double diagonal, double neighbor_sum)
{
    if (!std::isfinite(diagonal) || !std::isfinite(neighbor_sum) ||
        diagonal <= 0.0)
        throw std::invalid_argument(
            "piso_correction_gain: diagonal must be positive");
    const double denominator = diagonal - neighbor_sum;
    // A non-positive effective diagonal is an invalid algebraic operator, not
    // a small denominator to be silently regularized. Clamping it here can
    // manufacture an arbitrarily large PISO gain and conceal a sign error.
    if (!(denominator > 0.0) || !std::isfinite(denominator))
        throw std::invalid_argument(
            "piso_correction_gain: effective diagonal must be finite and positive");
    return 1.0 / denominator;
}


// The fractional-step pressure projection is deliberately explicit about its
// contract. It is a projection of the conservative face flux, not an alias for
// PISO. Keeping this helper here also makes the algorithm choice visible to
// callers without embedding policy in the finite-volume transport layer.
inline bool is_fractional_step_algorithm(PressureVelocityAlgorithm a)
{
    return a == PressureVelocityAlgorithm::FRACTIONAL_STEP;
}

inline bool is_coupled_algorithm(PressureVelocityAlgorithm a)
{
    return a == PressureVelocityAlgorithm::COUPLED;
}

}  // namespace cfdx::physics
