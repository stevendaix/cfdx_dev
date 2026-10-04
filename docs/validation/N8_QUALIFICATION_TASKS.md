# N8 Qualification Task List

> Status: **qualification in progress**. This document is a work list, not a qualification claim.

N8 implementation is substantially merged. The remaining work is to establish reproducible numerical evidence on production-relevant CFD operators and to close the documented qualification boundaries.

## 1. FGMRES gate

- [x] Merge and validate PR #598.
- [x] Confirm `test_fgmres_solver` is present in CTest.
- [x] Confirm the N8 campaign discovers and requires the dedicated FGMRES test.
- [x] Confirm the variable-preconditioner test passes with an independently recomputed true residual.

**Exit evidence:** explicit FGMRES dispatch + dedicated variable-preconditioner qualification test + green CI.

## 2. Production physical campaign

Cases:
- [ ] Couette
- [ ] Poiseuille
- [ ] Ghia cavity, Re=100
- [ ] Controlled skew/non-orthogonal case

For each selected solver/preconditioner combination:
- [ ] convergence status
- [ ] nonlinear iterations
- [ ] Krylov iterations
- [ ] independently recomputed true residual
- [ ] continuity/conservation metrics
- [ ] physical QoIs
- [ ] setup cost
- [ ] solve cost
- [ ] memory/hierarchy diagnostics where applicable

**Exit evidence:** machine-readable campaign report with reproducible case/solver configuration and numerical results.

## 3. Pressure-velocity algorithm matrix

Explicitly execute and retain evidence for:
- [ ] SIMPLE
- [ ] SIMPLEC
- [ ] PISO
- [ ] PIMPLE
- [ ] Fractional Step
- [ ] COUPLED

Do not treat API availability alone as qualification.

**Exit evidence:** every claimed production algorithm has at least one representative converged case with objective numerical checks.

## 4. Krylov solver matrix

- [x] CG
- [x] BiCGStab
- [x] GMRES
- [x] FGMRES

For each solver:
- [ ] true residual check
- [ ] convergence/failure behavior
- [ ] appropriate matrix class
- [ ] documented applicability restrictions

## 5. Schur / LSC / BFBT qualification

Current algebraic tests are infrastructure evidence, not final CFD qualification.

- [ ] Build representative pressure-velocity CFD matrices.
- [ ] Compute exact Schur reference where tractable.
- [ ] Compare SIMPLE/SIMPLEC, LSC and BFBT against the exact reference.
- [ ] Measure approximation error and conditioning.
- [ ] Include null-space handling where applicable.
- [ ] Establish an acceptance envelope from representative matrices.
- [ ] Identify the FP64 numerical floor from measured results.
- [ ] Test numeric-value refresh without rebuilding the graph.
- [ ] Document failure modes instead of weakening thresholds.

**Exit evidence:** representative CFD matrix dataset + quantitative comparison + justified acceptance envelope.

## 6. AMG difficult-matrix qualification

- [ ] Isotropic baseline.
- [ ] Strong anisotropy sweep.
- [ ] Strong coefficient-scaling sweep.
- [ ] Increasing matrix size.
- [ ] Native AMG.
- [ ] Smoothed Aggregation AMG.
- [ ] Record hierarchy diagnostics.
- [ ] Record setup/solve cost.
- [ ] Record convergence degradation and failure modes.
- [ ] Separate diagnostic limits from hard qualification gates.

**Exit evidence:** reproducible sweep and documented robustness/failure envelope.

## 7. Production preconditioner matrix

Qualify the production paths actually exposed by CFDX:

- [ ] Native AMG
- [ ] Smoothed Aggregation AMG
- [ ] Native FieldSplit
- [ ] Coupled Block Schur
- [ ] MGR

For each applicable path:
- [ ] representative operator
- [ ] true residual
- [ ] convergence behavior
- [ ] setup/solve cost
- [ ] documented limitations

## 8. Production solver policy

- [ ] Explicit pressure solver policy.
- [ ] Explicit momentum solver policy.
- [ ] Explicit scalar solver policy.
- [ ] Explicit coupled-system policy.
- [ ] No silent solver substitution.
- [ ] No silent CPU/GPU fallback.
- [ ] No automatic relaxation of numerical tolerances.
- [ ] Unsupported combinations fail explicitly and diagnostically.

**Exit evidence:** documented production policy and tests proving unsupported configurations do not silently change method.

## 9. MPI boundary

- [ ] Run representative Schur/AMG production cases in serial.
- [ ] Run the same cases under MPI where supported.
- [ ] Compare convergence and numerical results.
- [ ] Check deterministic/reproducible reductions where required.
- [ ] Measure scaling only if the implementation is mature enough to support the claim.
- [ ] Otherwise document the exact N8 MPI boundary and follow-up work.

**Exit evidence:** explicit qualification boundary; no implied MPI qualification from unit tests alone.

## 10. Final qualification evidence

- [ ] Full N8 physical campaign passes.
- [ ] All mandatory gates are present in CTest.
- [ ] No validation case is disabled merely to obtain a green build.
- [ ] No tolerance inflation.
- [ ] No silent fallback.
- [ ] Numerical maturity passes on the exact HEAD.
- [ ] Numerical maturity audit passes on the exact HEAD.
- [ ] CFDX CI passes on the exact HEAD.
- [ ] Machine-readable evidence artifacts are retained.
- [ ] Issue #461 is updated with implemented vs qualified status.
- [ ] N8 is promoted to `QUALIFIED` only after all mandatory evidence is reviewed.

## Recommended execution order

1. **#598 FGMRES** — closed and verified.
2. **Production physical campaign** — establish the baseline evidence matrix.
3. **LSC/BFBT on representative CFD matrices** — N8 Schur qualification.
4. **AMG difficult-matrix campaign** — N8 AMG robustness envelope.
5. **Production solver/preconditioner policy** — close method-selection gaps.
6. **MPI boundary** — qualify or explicitly bound.
7. **Final exact-HEAD campaign + audit** — only then close N8.

## Qualification rule

Passing the orchestration script is **not** sufficient for N8 qualification. The campaign must produce numerical evidence that is representative of the production CFD operators and satisfies the acceptance criteria without tolerance inflation, disabled validation, or silent method substitution.
