# NACA refinement divergence: exact-HEAD diagnostic protocol

Status: diagnostic plan; no production solver behavior is changed by this document.

Tracked by [issue #986](https://github.com/stevendaix/cfdx_dev/issues/986).

## Scope

Investigate both long-validation failures:

- `test_naca0012_qualification`
- `test_naca0012_wall_shear_audit`

The second may be a cascade from the first. Treat it as independent only if its own solve/gradient path produces a reproducible failure after the NACA flow solve is qualified.

Historical evidence from an older HEAD reported N64 converged, followed by N128 early divergence (`momentum_rhs_relative=17217.3`, `continuity=1.36922`, `alpha_u=0.6`). This is a lead, not proof for the current HEAD. Pin and record the current master SHA and exact run before drawing conclusions.

## A. Capture the production trajectory without changing it

Use the existing production SIMPLE path, existing mesh fixtures, initialisation, boundary conditions, discretisation, linear-solver settings and acceptance gates. Run N64, N128 and N256 independently where possible so an early failure on one mesh does not suppress evidence from the others.

Record at each nonlinear iteration:

- iteration index and mesh cell count;
- momentum RHS-relative and reference-relative residuals;
- pressure linear residual and linear solver convergence/iteration count;
- continuity L1, L-infinity and normalized residuals;
- velocity and pressure infinity-norm changes;
- effective `alpha_u` and `alpha_p`;
- finite-value checks and the first failing cell/face if available.

The existing NACA test enables `iteration_trace` at frequency 200, which is too coarse to locate an onset around iteration 11. For the diagnostic-only run, capture every iteration or use the already-retained nonlinear history; do not change the production algorithm. Keep raw logs and a machine-readable per-iteration table as artifacts.

### Checks before attributing the failure to relaxation

1. Verify mesh topology, cell volumes, face areas/orientation and local geometric coefficients on N64/N128.
2. Check the first cells where momentum residuals grow for finite/sign-consistent momentum coefficients and convection/diffusion contributions.
3. Check the pressure-correction matrix, pressure gauge and fixed-pressure boundary handling against the actual face-flux correction.
4. Independently evaluate the continuity/flux identity before and after pressure correction.
5. Compare momentum residual growth against pressure residual, continuity and linear-solver status. A pressure solve that reports convergence is not sufficient if the assembled correction is inconsistent with the flux operator.
6. Compare the wall-shear audit's failure point with the flow solve verdict. Do not call a downstream audit failure an independent root cause without evidence.

## B. Classify nonlinear convergence states; do not alter the solution

Add or validate a convergence classification that distinguishes:

- `CONVERGED`: all existing authoritative gates pass;
- `STAGNATED`: insufficient improvement over the configured stagnation window without sustained explosive growth;
- `DIVERGED`: sustained, material residual growth over a bounded window, or non-finite state/metrics;
- `MAX_ITERATIONS`: finite trajectory reaches the configured budget without convergence or a stronger classification.

The exact thresholds/window must be justified against existing convergence controls and synthetic traces, not tuned to make NACA pass. Do not classify a single noisy increase as divergence. Preserve the first triggering iteration, metric, value, comparison window and reason in the result.

Required unit tests: convergent decay, flat stagnation, monotone growth, exponential growth, oscillatory-but-bounded residuals, transient spike followed by recovery, zero reference metric, NaN/Inf, and exact precedence between `DIVERGED` and `MAX_ITERATIONS`.

The classification is diagnostic only: it must not change `alpha_u`, `alpha_p`, the linear solver, residual normalization, convergence tolerances, force oracle, iteration budget or fallback behavior.

## Acceptance

- N64 retains its existing result and authoritative gates.
- N128 divergence is detected at its measured onset and not mislabeled as ordinary stagnation.
- N256 is run independently if the harness permits, and its verdict is reported separately.
- The wall-shear audit is attributed as independent or cascading based on its own evidence.
- Release, DebugSanitizers, focused NACA tests and exact-HEAD `validation-total` results are reported. A quick green is not sufficient for merge.
- No tolerance, oracle or budget is weakened; no test is skipped.

## Separate follow-up: adaptive relaxation (not part of A/B)

Only after A/B is complete, compare fixed `alpha_u` values on the same N64/N128 meshes and solver path. Present per-iteration residuals, continuity/flux identity, linear-solver status and force results for each run. Adaptive relaxation is justified only if the evidence isolates relaxation as causal and the flux/matrix checks pass.

Any production change must be proposed as a separate diff and reviewed before implementation. No adaptive relaxation or other production solver change is authorized by this diagnostic work.
