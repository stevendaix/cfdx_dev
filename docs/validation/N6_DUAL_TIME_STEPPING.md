# N6 dual-time stepping

Dual-time stepping separates physical-time discretisation from pseudo-time convergence of the nonlinear system.

The implementation is intentionally delivered as a reusable numerical contract first; production Navier-Stokes wiring is a follow-up integration step.


## Mathematical contract

Backward Euler:
\[
R(U^{n+1}) = U^{n+1}-U^n-\Delta t F(U^{n+1})=0.
\]

BDF2:
\[
R(U^{n+1}) = \frac32 U^{n+1}-2U^n+\frac12 U^{n-1}-\Delta t F(U^{n+1})=0.
\]

The pseudo-time iterations solve this physical residual. They do not advance physical time.

## Implementation

- src/cfdx/physics/dual_time_stepping.h defines the physical residual, convergence controls, diagnostics and failure type.
- src/cfdx/physics/dual_time_driver.h provides the deterministic pseudo-time loop.
- tests/unit/test_dual_time_stepping.cpp verifies implicit BDF2 and backward-Euler roots and deterministic non-convergence failure.
- Pseudo-time adaptation is bounded and based only on residual evolution.
- No convergence is accepted merely because the iteration budget was reached.

## Scope boundary

This PR is a reusable numerical contract, not a claim that the complete incompressible/compressible Navier–Stokes transient application is already wired to dual time. The current production transient integration point is not a single mature NS loop, so that coupling should be a separate PR with field-level residuals, nonlinear solver integration, checkpoint history and transient MMS.

## Literature direction

Recent work supports BDF1/BDF2 dual-time integration, adaptive pseudo-time strategies, and physics-specific positivity/entropy safeguards. See DUAL_TIME_STEPPING_BIBLIOGRAPHY.md. The CFDX core deliberately stays deterministic; ML-based local pseudo-time prediction remains research-only.
