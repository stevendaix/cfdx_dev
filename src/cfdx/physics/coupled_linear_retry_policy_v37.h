#pragma once

// The outer U/p under-relaxation shrinks only AFTER the COUPLED saddle-point
// linear solve. On the SAME rejected outer state, that shrink cannot improve
// an already failed linear solve unless it changes the linear preconditioner
// (SIMPLEC Schur) or an explicit per-retry pseudo-time operator (V31/V32).
// This guard saves identical expensive Krylov replays; it does not change
// linear tolerances, residuals, or the nonlinear convergence certificate.
namespace cfdx::physics {

struct CoupledLinearRetryV37 {
    bool shrink_can_change_linear_problem = false;
    bool reason_is_coupled_linear_failure = false;
    bool changed_predictor_requires_picard_fallback = false;

    bool retry_can_help() const noexcept {
        return reason_is_coupled_linear_failure &&
            (shrink_can_change_linear_problem ||
             changed_predictor_requires_picard_fallback);
    }
    bool fail_fast() const noexcept {
        return reason_is_coupled_linear_failure && !retry_can_help();
    }
};

} // namespace cfdx::physics
