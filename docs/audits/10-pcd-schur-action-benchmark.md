# PCD Schur-action A/B benchmark on the coupled Ghia system

## Objective

Measure the cost and linear-solver behavior of alternative PCD Schur-complement actions on the actual assembled coupled pressure-velocity systems from the Ghia lid-driven cavity qualification. This is a diagnostic/benchmark method, not a change to the production preconditioner or a physical qualification.

## Method

1. Capture the matrix, right-hand side, Rhie–Chow pressure-Schur diagonal, and linearization velocity from the production coupled assembly. The optional callback is null by default and must not alter solver behavior.
2. Replay selected captured stages through the baseline and alternative Schur actions using the same matrix, right-hand side, initial guess, outer Krylov controls, and convergence tolerance.
3. Recompute the true residual independently as ||b - A x|| / ||b||. Report outer iterations, nested-solve iterations and time, and setup time separately.
4. Require the baseline action to reproduce the production linear convergence target on every captured stage. Require converged alternatives to agree with the baseline solution within the stated relative comparison gate.
5. Keep physical acceptance separate: the existing Ghia qualification test remains the sole physical-case gate, with its existing oracle, tolerance, and timeout unchanged.

## Interpretation limits

- A faster action on captured stages is evidence about linear-solver cost, not by itself evidence of nonlinear convergence or physical correctness.
- A non-baseline action that fails to converge is recorded as a failed measurement, not silently substituted or treated as a pass.
- No result from this benchmark authorizes changing the production PCD action. Such a change requires an explicit follow-up, full coupled qualification, and no regression in the physical oracle.
- The branch-only benchmark workflow is evidence collection and is not intended to be merged into master as-is. The test registration and read-only assembly hook require separate review before any production integration.

## Traceability

- PR #966: PCD Schur-action A/B benchmark on real Ghia coupled systems.
- Issue #914: audit of Fp versus Kp Schur actions.
- Issue #390: convergence robustness.

## Reference

- Kay, D., Loghin, D., and Wathen, A. J. (2002). “A preconditioner for the steady-state Navier–Stokes equations.” *SIAM Journal on Scientific Computing*, 24(1), 237–256. DOI: [10.1137/S1064827500384937](https://doi.org/10.1137/S1064827500384937). This provides background on pressure-convection-diffusion preconditioning; it does not validate the particular CFDX action or this benchmark's measured outcomes.
