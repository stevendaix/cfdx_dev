# SST production-solver implementation tracker

**Status date:** 2026-10-19  
**Overall status:** IN PROGRESS — RANS production selection remains fail-closed.  
**Scope:** Track implementation and qualification evidence for end-to-end SST in CFDX. Update this file in each turbulence PR; do not mark a gate complete without exact-HEAD evidence. Numerical-core implementation is governed by [the RANS/SST method fiche](../audits/09-turbulence-rans.md).

## Current baseline

- [x] Production driver rejects unsupported non-laminar turbulence selections instead of silently solving laminar equations (merged PR #941).
- [x] SST transport component APIs and dynamic blending exist; they are not yet coupled to the production outer loop.
- [x] Cellwise momentum-viscosity assembly overload exists and has uniform-coefficient equivalence coverage.
- [x] Production solver consumes a cellwise turbulent-viscosity field consistently in predictor, final nonlinear residual, and independent momentum conservation reconstruction (PR #967 merged; end-to-end SST remains disabled).
- [ ] SST state `k`, `omega`, wall distance and turbulence residuals are owned by the production solve and participate in rollback/convergence.
- [x] Add a non-production SST case preflight for explicit positive initial `k`/`omega` and validated inlet turbulence intensity/length scale; derive the standard inlet estimates without enabling the production selector.
- [x] Represent explicit per-patch `k`/`omega` boundary-condition types in the persistent case schema and HDF5 round-trip, with compatibility validation.
- [ ] Define restart semantics and wire validated scalar boundary contracts into the production transport solve.
- [ ] Production-path transport MMS on three or more meshes.
- [ ] Frozen DNS channel qualification (Moser–Kim–Mansour), with mesh/y+ sensitivity and independent balance checks.
- [ ] Flat-plate and NACA0012 turbulent validation after channel qualification.

## Work packages

### T0 — Cellwise effective viscosity plumbing

**Intent:** Establish one authoritative cellwise dynamic-viscosity vector,
`mu_eff[c] = rho * (nu + nu_t[c])`, and use it in every momentum equation path.

- [x] Add an optional cellwise kinematic turbulent-viscosity input to solver controls.
- [x] Reject wrong-sized, negative, or non-finite cell values before mutating the solution (unit regression added).
- [x] Reject ambiguous simultaneous nonzero scalar and cellwise turbulent-viscosity inputs (unit regression added).
- [x] Use the same effective-viscosity vector for all three momentum predictors.
- [x] Use that vector when reconstructing the final momentum equations used for nonlinear acceptance.
- [x] Use the same vector in the independent post-solve momentum flux/balance audit.
- [x] Regression: uniform cellwise viscosity is algebraically equivalent to the legacy scalar-viscosity path (existing assembly regression).
- [x] Regression: invalid cellwise viscosity is rejected.
- [x] Exact-head CI/build/test evidence reviewed for PR #967 before merge.

**Merged prerequisite:** [#967 — cellwise effective viscosity](https://github.com/stevendaix/cfdx_dev/pull/967) (merged).  
**Known risk:** this plumbing alone does not enable SST and must not bypass the production RANS guard.

### T1 — SST fields and boundary contract

- [x] Add preflight validation for explicit initial `k`/`omega` and inlet intensity/length-scale inputs, plus a standard inlet estimate helper (unit tests added; not wired into production selection).
- [x] Add explicit persisted SST scalar-boundary metadata for fixed values, zero gradients, zero wall `k`, and SST wall `omega`; reject missing, duplicate, and patch-incompatible declarations (component tests; not yet consumed by production transport assembly).
- [ ] Define persistent production state for `k`, `omega`, and derived `nu_t`.
- [ ] Define initialization and restart contracts; reject missing/inconsistent fields.
- [ ] Define physically consistent wall/inlet/outlet boundary values for both transports.
- [ ] Compute wall distance from the production mesh and validate positivity/finiteness.
- [ ] Compute boundary-aware velocity gradients and SST strain-rate invariants.
- [ ] Compute `F1`/`F2` from the current accepted state and documented SST equations.
- [ ] Add state rollback for every turbulence field, not just `U` and `p`.

### T2 — Coupled outer iteration and convergence

- [ ] Assemble/solve `k` and `omega` using the current authoritative mass flux.
- [ ] Apply declared under-relaxation and positivity/boundedness handling without clipping away failures.
- [ ] Recompute `nu_t` from accepted `k`, `omega`, strain and wall-distance fields.
- [ ] Reassemble momentum using the updated cellwise `nu_t`.
- [ ] Include normalized turbulence residuals and physical boundedness checks in convergence acceptance.
- [ ] Ensure rejected iterations roll back `U,p,k,omega,nu_t,phi` atomically.
- [ ] Add production E2E test proving SST selection invokes the coupled path and emits model/iteration evidence.
- [ ] Keep unsupported models fail-closed.

### T3 — Transport verification and physical qualification

- [ ] Introduce a reviewed manufactured-source seam through the actual SST transport assembly.
- [ ] Verify each transport equation using independently derived exact source and boundary data.
- [ ] Run at least three mesh levels and report volume-weighted L1, L2, Linf and observed order.
- [ ] Freeze the DNS dataset/version/hash, normalization, sampling plane and comparison metrics before running.
- [ ] Validate turbulent channel: mean velocity, Reynolds stresses where available, friction velocity, y+ distribution, mass balance and pressure-gradient/bulk-flow balance.
- [ ] Demonstrate mesh/y+ sensitivity and solver-iteration convergence.
- [ ] Validate flat plate, then turbulent NACA0012 only after preceding gates pass.

## Related work

- PR #941 — production guard and audit baseline (merged).
- PR #970 — SST strain-rate utility, merged as `374c8fb4` (not production coupling).
- PR #971 — explicit boundary-value Gauss gradient API, merged as `257070a6` (Dirichlet face-value contract only; not full wall treatment).
- PR #972 — closed as superseded by #973; its tracker correction is included in this branch.
- Issue #943 — end-to-end SST coupling plan (created by the prior work; verify issue state/number if the tracker is moved).

### T1 implementation increment — production-owned state contract

- [x] Introduce an explicit `SSTProductionState` bundle for cellwise `k`, `omega`, `nu_t`, wall distance and strain rate.
- [x] Validate field dimensions and finite/physical value ranges before the bundle is used.
- [x] Add explicit value snapshot/restore helpers and unit regressions covering every owned field.
- [ ] Wire this state into `solve_steady_incompressible`; this PR does not claim production-loop ownership yet.
- [ ] Define initialization and restart semantics from case/HDF5 data; no implicit defaults or guessed values.
- [ ] Extend the outer nonlinear transaction to include this state and the authoritative face mass flux.

This is an enabling contract only. It does not call SST transport from the production loop, enable the RANS selector, or qualify the physical model.

### T1 implementation increment — atomic carrier/SST transaction primitive

- [x] Add a transaction primitive that snapshots/restores the carrier velocity, pressure, authoritative face mass flux, and all fields in `SSTProductionState` as one unit.
- [x] Test exact rollback across all state groups and verify commit retains candidate fields.
- [ ] Connect this primitive to the production nonlinear retry/acceptance loop; the current increment does not call it from `solve_steady_incompressible`.
- [ ] Define HDF5/restart initialization and construct the state from the actual production mesh.
- [ ] Keep the production turbulence selector fail-closed until transport, viscosity refresh, and convergence gates are integrated and verified.

## Evidence log

| Date | Change / evidence | Decision |
|---|---|---|
| 2026-10-19 | Source audit: production solver used one scalar `mu_eff` in predictor, final residual, and independent conservation reconstruction; cellwise assembly overload already exists. | T0 first. Do not enable RANS yet. |
| 2026-10-19 | PR #967 merged as `03ac6f6c`; cellwise effective-viscosity plumbing is integrated. | T0 complete; no RANS enablement implied. |
| 2026-10-19 | PRs #970/#971 passed exact-HEAD CI and numerical maturity workflows, then merged as `374c8fb4` and `257070a6`. | Strain-rate and explicit Dirichlet-face-gradient prerequisites integrated; physical and end-to-end SST gates remain open. |
| 2026-10-19 | Added `sst_case_contract.h` to validate explicit initial k/omega and inlet intensity/length scale, with unit coverage. | Component/preflight only; production guard unchanged. |
| 2026-10-19 | Added typed per-patch SST scalar BC metadata and HDF5 round-trip coverage for inlet fixed values and wall `k=0`/SST-omega wall treatment. Compatibility validator rejects incomplete/duplicate or patch-incompatible contracts. | Persistent metadata only; restart semantics and production application remain open. |

## Merge/promotion rules

- No gate is complete merely because a component unit test passes.
- No merge while required CI is pending or unclear.
- Do not widen tolerances, disable tests, silently fall back to laminar, or substitute a laminar case for DNS evidence.
- Keep the production turbulence guard until all required coupling and qualification gates pass.
