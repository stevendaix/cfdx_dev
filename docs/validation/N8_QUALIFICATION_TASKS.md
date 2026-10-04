# N8 Qualification Task List

> Status: **qualification in progress**. This document is a work list, not a qualification claim.

N8 implementation is substantially merged. The remaining work is to establish reproducible numerical evidence on production-relevant CFD operators and to close the documented qualification boundaries.

> **Re-baselined 2026-10-04.** Section boxes below were reconciled against the merged tree. Items merged after this list was written were checked only where the current source tree evidences them; sections that remain unchecked say what is missing. Three corrections matter for reading the rest of this document:
>
> - **LSC and BFBt are implemented, not missing** (`src/cfdx/core/linalg/lsc_bfbt_schur.h`, #562), verified against the exact-Schur oracle (#565) and measured for conditioning and FP64 error floor (#568). Only **PCD** is unimplemented.
> - **The generic Schur layer is not production-selectable.** No file under `src/` or `apps/` includes it and no case-file setting chooses a Schur approximation; production still uses the separate `CoupledBlockSchur*` pair. Section 7 cannot close until that changes.
> - **SIMPLE/SIMPLEC lacked a pressure null-space policy and a graph-signature guard on `update_values`**; #634 gives it both, matching the exact and LSC/BFBt approximations.

## 1. FGMRES gate

- [x] Merge and validate PR #598.
- [x] Confirm `test_fgmres_solver` is present in CTest.
- [x] Confirm the N8 campaign discovers and requires the dedicated FGMRES test.
- [x] Confirm the variable-preconditioner test passes with an independently recomputed true residual.

**Exit evidence:** explicit FGMRES dispatch + dedicated variable-preconditioner qualification test + green CI.

## 2. Production physical campaign

Cases:
- [x] Couette
- [x] Poiseuille
- [x] Ghia cavity, Re=100
- [x] Controlled skew/non-orthogonal case

All four execute in the campaign and appear in `results[]` with their ctest status and output. **They do not yet emit structured numerical records.** `test_n8_pressure_velocity_matrix` emits `MODEL_CONFIG`/`MODEL_PLAN`/`MODEL_RESULT`; Poiseuille, Ghia and skew emit no records the report parses, so they contribute pass/fail and elapsed time only. Couette is covered through that matrix rather than by `test_couette_quick`, which the campaign does not require.

Each is run at its `--quick` variant; the full non-quick variants are gated behind `CFDX_ENABLE_LONG_VALIDATION`, so this campaign does not exercise the long-validation acceptance paths.

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
- [x] SIMPLE
- [x] SIMPLEC
- [x] PISO
- [x] PIMPLE
- [x] Fractional Step
- [x] COUPLED

All six run through `test_n8_pressure_velocity_matrix` (no `--quick`), which announces each model's requested configuration, records the plan the dispatcher resolved, and reports `solver_converged`, `iterations` and `gates_failed` per model. `model_resolution` reports `tallies_agree` between the announced and observed counts. COUPLED is exercised under three preconditioners (`BlockSchur`, `MGR` and `PCD`), giving eight configured models.

`PCD` passes this matrix, with the same gates as the other coupled preconditioners. Its action and reference-cell gauge were corrected after the case was first exercised; see [PCD_SCHUR_QUALIFICATION.md](PCD_SCHUR_QUALIFICATION.md#production-path-evidence). Its operator ordering remains unqualified on a non-zero-convection case.

API availability alone was not treated as qualification: the evidence is the per-model resolved plan and verdict.

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
- [x] Compute exact Schur reference where tractable. (Implicit exact-Schur oracle with a preconditioner-free CG action, #492; independent dense-assembly cross-check, #565.)
- [x] Compare SIMPLE/SIMPLEC, LSC and BFBT against the exact reference. (`test_schur_approximation_comparison`.)
- [x] Measure approximation error and conditioning. (Infinity-norm `cond(Auu)`, oracle discrepancy, measured FP64 backward-error floor, per-method `algebra_error`/`exact_schur_error`, #568. Approximation quality is deliberately left diagnostic: no envelope exists.)
- [x] Include null-space handling where applicable. Delivered for the exact and LSC/BFBt approximations (`test_lsc_bfbt_schur_null_space`, #575) and for SIMPLE/SIMPLEC (#634). The Schur layer still offers only the constant/mean-zero projector — there is no pinned-pressure policy there, so that part remains open.
- [ ] Establish an acceptance envelope from representative matrices. Blocked on the first item; the controlled family measures LSC/BFBt action error at 4.8–8.6 against exact, which is why no threshold may be invented yet.
- [x] Identify the FP64 numerical floor from measured results. (Machine epsilon and attainable backward-error floor reported per case, #568; #482 established the `tol >> eps*cond(A)` rule.)
- [x] Test numeric-value refresh without rebuilding the graph. (`n8_schur_benchmark_lifecycle`: value-only change on an unchanged graph, `hierarchy_builds` held constant, one `numeric_update`, true residual independently recomputed.) Graph-*change* rejection is guarded for the exact, LSC/BFBt and SIMPLE/SIMPLEC approximations.
- [x] Document failure modes instead of weakening thresholds.

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

- [x] Native AMG
- [x] Smoothed Aggregation AMG
- [x] Native FieldSplit
- [x] Coupled Block Schur
- [x] MGR

Each has a production-path test the campaign requires and runs: `test_amg_preconditioner_qualification`, `test_advanced_preconditioners`, `test_schur_preconditioner`, `test_coupled_block_schur_amg` and `test_mgr_preconditioner`.

**This closes the "the path is exercised" question only.** The per-path items below are not met:

For each applicable path:
- [ ] representative operator
- [ ] true residual
- [ ] convergence behavior
- [ ] setup/solve cost
- [ ] documented limitations

The campaign records true residual, Krylov iterations, setup/solve time and NNZ for the coupled Schur *algebraic* family (`test_n8_schur_production_benchmark`), which is not the production preconditioner set, and records no per-path cost or limitation summary for the five paths above. Separately, the generic `SchurApproximation` layer (exact oracle, SIMPLE/SIMPLEC, LSC/BFBt) remains unreachable from production: no case-file or solver-options setting selects it, so the approximations qualified in section 5 are not the ones a production case can run.

## 8. Production solver policy

- [ ] Explicit pressure solver policy.
- [ ] Explicit momentum solver policy.
- [ ] Explicit scalar solver policy.
- [ ] Explicit coupled-system policy.
- [x] No silent solver substitution. `audit_linear_plan` compares each announced linear request against the plan the dispatcher resolved, for the Krylov method *and* the preconditioner of both the coupled and the pressure sub-problem. `substitutions`, `models_without_plan` and `structure_mismatches` are wired into the campaign status, making `explicit_linear_request_must_be_honored` the one enforced gate in the report (#620).
- [ ] No silent CPU/GPU fallback.
- [ ] No automatic relaxation of numerical tolerances.
- [ ] Unsupported combinations fail explicitly and diagnostically.

The remaining three are unverified for N8 specifically. `PreconditionerModel::LSC` exists in the model enum and catalogue but `make_scalar_preconditioner` rejects block preconditioners and the scalar dispatcher throws `not available in scalar dispatch`; that is an explicit failure, but the enum entry advertising an unavailable path is itself an unresolved item rather than evidence for this section.

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
- [x] All mandatory gates are present in CTest. The campaign is registered as `test_n8_physical_qualification` (`CMakeLists.txt`), all 22 required tests resolve against `ctest -N`, and the audit validator runs in the `cfdx-maturity` workflow.
- [ ] No validation case is disabled merely to obtain a green build.
- [ ] No tolerance inflation. Both are claims about development history rather than properties of the current tree, so they are not checked here. The report asserts `changes_numerical_tolerances: false` and `disables_validation: false` in its own policy block, but a self-assertion is not independent evidence; this item stays open until reviewed against the diffs.
- [ ] No silent fallback. No silent *solver substitution* is enforced (section 8); silent CPU/GPU fallback is unverified.
- [ ] Numerical maturity passes on the exact HEAD.
- [ ] Numerical maturity audit passes on the exact HEAD.
- [ ] CFDX CI passes on the exact HEAD.
- [x] Machine-readable evidence artifacts are retained.
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
