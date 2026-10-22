# SST production-solver implementation tracker

**Status date:** 2026-10-21  
**Overall status:** IN PROGRESS — RANS production selection remains fail-closed.  
**Scope:** Track implementation and qualification evidence for end-to-end SST in CFDX. Update this file in each turbulence PR; do not mark a gate complete without exact-HEAD evidence. Numerical-core implementation is governed by [the RANS/SST method fiche](../audits/09-turbulence-rans.md).

## Current baseline

- [x] Production driver rejects unsupported non-laminar turbulence selections instead of silently solving laminar equations (merged PR #941).
- [x] SST transport APIs, dynamic blending, and explicit opt-in production outer-loop coupling are present in PR #1032; production E2E and physical qualification remain open.
- [x] Cellwise momentum-viscosity assembly overload exists and has uniform-coefficient equivalence coverage.
- [x] Production solver consumes a cellwise turbulent-viscosity field consistently in predictor, final nonlinear residual, and independent momentum conservation reconstruction (PR #967 merged; end-to-end SST remains disabled).
- [x] The solve accepts an explicit opt-in `SSTProductionRuntime`; SST state is validated and snapshotted per outer attempt, transport is called with the current pressure-corrected mass flux, and `nu_t` refreshes `mu_eff` for the next predictor (PR #1032 pending merge; exact-head CI green).
- [x] Add a non-production SST case preflight for explicit positive initial `k`/`omega` and validated inlet turbulence intensity/length scale; derive the standard inlet estimates without enabling the production selector.
- [x] Represent explicit per-patch `k`/`omega` boundary-condition types in the persistent case schema and HDF5 round-trip, with compatibility validation.
- [ ] Connect persisted case scalar-boundary metadata and DAT state restoration to the application driver's runtime construction; the current runtime requires the caller to supply validated state and boundary contracts explicitly.
- [ ] Production-path transport MMS on three or more meshes.
- [ ] Frozen DNS channel qualification (Moser–Kim–Mansour), with mesh/y+ sensitivity and independent balance checks.
- [ ] Flat-plate and NACA0012 turbulent validation after channel qualification.

## Latest implementation increment

- [x] Add `initialize_sst_production_state_from_mesh(...)`: production mesh geometry supplies exact BVH wall distance and boundary-aware strain rate to SST state initialization.
- [x] Add a regression for computed positive wall distance, nonzero strain from fixed inlet data, initialized `k`/`omega`, derived `nu_t`, and invalid velocity-field shape.
- [x] Production outer-loop ownership/advance is wired through opt-in SSTProductionRuntime in PR #1032. Production E2E and application case/restart wiring remain open.

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
- [x] Define a validated in-memory production-state bundle for `k`, `omega`, and derived `nu_t` (PR #993); production-loop ownership is still open.
- [x] Define component-level initialization and DAT restart contracts; reject missing/inconsistent fields (PRs #1006, #1007, #1008, #1009). Production lifecycle wiring and trajectory-equivalence qualification remain open.
- [x] Define persisted scalar-boundary metadata and runtime adapters, including resolved-wall omega face values (PRs #998, #1000); case-to-live-production transport wiring remains open.
- [x] Compute wall distance from the production mesh and validate positivity/finiteness (PR #1024; exact-HEAD unit/CI evidence passed).
- [x] Provide a boundary-aware SST strain-rate reconstruction using explicit fixed-value velocity patch data; retain the existing zero-gradient overload unchanged.
- [ ] Wire production velocity BCs into the boundary-aware strain-rate path and qualify gradient accuracy on representative meshes.
- [ ] Compute `F1`/`F2` from the current accepted state and documented SST equations.
- [x] Provide a reusable atomic transaction primitive covering `U`, `p`, authoritative `phi`, and all SST state fields (PR #995); connecting it to the production nonlinear retry loop remains open.

### T2 — Coupled outer iteration and convergence

- [x] The opt-in production path calls the transactional `k`/`omega` transport update with the current authoritative mass flux (PR #1032); production E2E and physical qualification remain open.
- [ ] Apply declared under-relaxation and positivity/boundedness handling without clipping away failures.
- [x] Refresh wall distance, boundary-aware strain rate, blending and `nu_t` from the current carrier field before the first predictor and during each SST transport step (PR #1032).
- [x] Rebuild `mu_eff[c] = rho * (nu + nu_t[c])` after an accepted coupled iteration for the next momentum predictor (PR #1032).
- [x] Add an independent SST convergence/boundedness hard gate without mixing dimensional transport residuals into the normalized carrier metric (PR #1032).
- [x] Per-attempt SST snapshot/restore is connected to carrier retry/rejection and unexpected-exception paths; nested rollback regression is covered by PR #1030.
- [ ] Add production E2E coverage proving the opt-in runtime executes transport, refreshes viscosity, and restores all coupled fields on a real solve rejection. Current PR adds runtime contract coverage only; do not mark this gate complete yet.
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
- [x] Wire the validated state into `solve_steady_incompressible` through the explicit opt-in runtime (PR #1032); case/restart construction and production E2E tests remain open.
- [ ] Wire the existing initialization/DAT restore helpers into application case/HDF5 lifecycle and prove continuous-versus-restarted production trajectory equivalence; no implicit defaults or guessed values.
- [x] Extend outer retry/rejection handling to restore SST state alongside carrier rollback and authoritative face mass flux (PR #1032; nested regression in #1030).

This is an enabling contract only. It does not call SST transport from the production loop, enable the RANS selector, or qualify the physical model.

### T1 implementation increment — atomic carrier/SST transaction primitive

- [x] Add a transaction primitive that snapshots/restores the carrier velocity, pressure, authoritative face mass flux, and all fields in `SSTProductionState` as one unit.
- [x] Test exact rollback across all state groups and verify commit retains candidate fields.
- [x] Connect SST snapshot/restore to production retry/acceptance in PR #1032; nested commit-then-outer-reject regression is in #1030.
- [ ] Define HDF5/restart initialization and construct the state from the actual production mesh.
- [x] Keep the general production turbulence selector fail-closed. The opt-in runtime does not enable the general RANS selector or imply physical qualification.

## Evidence log

| Date | Change / evidence | Decision |
|---|---|---|
| 2026-10-19 | Source audit: production solver used one scalar `mu_eff` in predictor, final residual, and independent conservation reconstruction; cellwise assembly overload already exists. | T0 first. Do not enable RANS yet. |
| 2026-10-19 | PR #967 merged as `03ac6f6c`; cellwise effective-viscosity plumbing is integrated. | T0 complete; no RANS enablement implied. |
| 2026-10-19 | PRs #970/#971 passed exact-HEAD CI and numerical maturity workflows, then merged as `374c8fb4` and `257070a6`. | Strain-rate and explicit Dirichlet-face-gradient prerequisites integrated; physical and end-to-end SST gates remain open. |
| 2026-10-19 | Added `sst_case_contract.h` to validate explicit initial k/omega and inlet intensity/length scale, with unit coverage. | Component/preflight only; production guard unchanged. |
| 2026-10-19 | Added typed per-patch SST scalar BC metadata and HDF5 round-trip coverage for inlet fixed values and wall `k=0`/SST-omega wall treatment. Compatibility validator rejects incomplete/duplicate or patch-incompatible contracts. | Persistent metadata only; restart semantics and production application remain open. |
| 2026-10-21 | PR #1024 merged as `63e195ac`; mesh-derived SST initializer uses production mesh wall distance and boundary-aware strain rate. Exact-head DebugSanitizers, Release, audit gate, numerical-maturity, and N10 fixture workflows passed. | Mesh-derived initialization prerequisite complete; live production-loop ownership, SST transport coupling, turbulence convergence, and physical qualification remain open. |

## Merge/promotion rules

- No gate is complete merely because a component unit test passes.
- No merge while required CI is pending or unclear.
- Do not widen tolerances, disable tests, silently fall back to laminar, or substitute a laminar case for DNS evidence.
- Keep the production turbulence guard until all required coupling and qualification gates pass.


### T1 implementation increment — validated DAT state hydration bridge

- [x] Read explicit DAT `k` and `omega` into candidate fields and reject checkpoints without `omega`.
- [x] Stage carrier `U`/`p` and SST fields on temporaries; malformed/incomplete restart must leave the caller's state unchanged.
- [x] Recompute strain-rate magnitude from the restored velocity and recompute `nu_t` from restored `k`/`omega`, current-mesh wall distance and SST blending.
- [ ] Integrate the bridge into the production solve lifecycle and prove continuous-vs-restarted trajectory equivalence.
- [ ] Recompute wall distance through the production wall-distance provider on restart and validate provider/operator equivalence.
- [ ] Keep the RANS/SST selector fail-closed until full coupled-loop and physical qualification gates pass.

This helper is a restart-state bridge only; it is not yet called by the production driver and does not qualify end-to-end SST.

### T1 implementation increment — explicit SST restart scalars

- [x] Extend the native DAT checkpoint format to version 4 with a stable, named `omega` scalar field alongside the existing `k` field.
- [x] Validate omega field shape and finiteness before writing; reject a requested omega field when the checkpoint does not contain it.
- [x] Add a round-trip regression for iteration/time, `k`, and `omega` on the same identified mesh.
- [ ] Wire this format into production SST restart initialization and verify restart-equivalent solution trajectories.
- [ ] Recompute derived `nu_t`, wall distance and strain rate from validated current fields/mesh after restart; do not treat cached derived values as authoritative.
- [ ] Add interrupted/corrupt-checkpoint atomicity coverage before declaring the production restart path complete.

This increment defines persistent scalar names only. It does not serialize the full `SSTProductionState`, establish restart-equivalent production behavior, or enable the RANS selector.

### T1 implementation increment — runtime boundary adapter (merged in PR #998)

- [x] Translate persisted SST fixed-value and zero-gradient conditions into the existing scalar transport boundary API.
- [x] Translate `WALL_K_ZERO` to a zero fixed-value condition for `k`.
- [x] Preserve `WALL_OMEGA_SST` as an explicit unresolved wall-operator obligation; never silently approximate it as fixed-value or zero-gradient.
- [x] Add unit regressions for inlet/outlet mapping and fail-closed wall-omega behavior.
- [ ] Apply this adapter in the production case path only after a qualified SST wall-omega operator exists.
- [ ] Couple boundary handling to production transport and restart state; the production RANS selector remains fail-closed.


### T2 implementation increment — resolved-wall omega face operator (PR #1000 merged)

- [x] Implement the documented resolved-wall value \(\omega_w=6\nu/(\beta_1 y_1^2)\) with finite-positive validation.
- [x] Build explicit per-face omega values from wall patch names and owner-cell-to-face normal distances; reject missing patches, degenerate faces, and patches with no faces.
- [x] Thread the optional face-value contract through SST transport and the transactional update without changing default behavior for existing callers.
- [x] Add a unit regression for the formula and invalid geometry inputs.
- [x] Add a mesh-level integration regression proving only requested wall faces receive finite values and non-wall faces remain unset.
- [x] Build and run targeted numerical tests and the required CI checks on exact HEAD `ff8da5a931e52389e32d777d3a6a9ce2a8f3dcda`.
- [x] Keep the production RANS selector fail-closed until state initialization/restart, full momentum coupling, MMS and turbulent-channel DNS qualification are complete.


### T1 implementation increment — explicit production-state initialization (PR #1006 merged)

- [x] Build uniform initial (k) and (omega) fields only from the explicit case scalar fields validated by the SST case contract.
- [x] Require mesh-sized, finite positive wall distance and finite non-negative strain rate supplied by the production caller; do not invent either field.
- [x] Compute initial (F_2) from the initialized fields and refresh (
u_t) using the SST closure before validating the complete state.
- [x] Add unit coverage for valid initialization and invalid wall-distance/strain-rate input.
- [ ] Wire the initializer to production case loading after restart semantics and wall-distance/strain-rate producers are defined.
- [ ] Add per-cell restart read/write and exact restart-equivalence coverage.
- [ ] Keep the production RANS selector fail-closed; this helper is not a coupled production solve or physical qualification.

| 2026-10-20 | PR #1000 merged as `62166fec`: resolved wall-(omega) face operator and mesh-level face assignment regression; all applicable required workflows on exact head passed. Long validation, thermal/radiation regression and physical SST campaigns were skipped by path filters. | Wall operator component complete; no RANS enablement. |
| 2026-10-20 | Began explicit production-state initializer: requires case-declared initial (k,omega), valid SST boundary metadata, and caller-supplied mesh-sized wall-distance/strain-rate fields; initializes closure (
u_t) and rejects invalid input. | Initialization contract only; production caller/restart/coupled iteration remain open. |


### T1 implementation increment — boundary-aware SST strain rate (PR #1015 merged)

- [x] Add an explicit patch-name to fixed velocity-value overload for SST strain-rate reconstruction.
- [x] Preserve zero-gradient extrapolation for patches without a fixed-value entry and exclude EMPTY extrusion faces from the gradient sum.
- [x] Reject non-finite velocity and boundary values, invalid face ownership, missing boundary patch assignment and invalid cell volumes.
- [x] Add regression coverage proving a fixed-value inlet changes the computed strain invariant and non-finite boundary values fail closed.
- [x] Verify an affine velocity field with fixed-value boundary data on an orthogonal box mesh (PR #1017).
- [x] Preserve fixed velocity boundary values when reconstructing strain rate from a DAT-restored velocity (PR #1018).
- [ ] Connect the production driver's declared velocity BCs to this overload for the live solve.
- [ ] Verify gradient accuracy with boundary-aware manufactured fields on representative non-orthogonal meshes before using it for SST production coupling.

This remains a prerequisite utility only. Restart reconstruction now accepts explicit fixed-value velocity patch data, but the live production nonlinear loop still does not call SST transport or use SST-updated viscosity.


### T1 implementation increment — restart boundary-data correctness (PR #1018 merged)

- [x] Expose an explicit fixed-velocity-boundary map on the DAT SST restoration helper.
- [x] Reconstruct restored strain rate with the boundary-aware gradient overload instead of silently assuming zero-gradient boundaries.
- [x] Add a DAT restart regression with non-zero inlet velocity and assert finite, positive restored strain rate.
- [ ] Wire restart restoration into the production driver's actual case/restart lifecycle.
- [ ] Prove continuous-versus-restarted production trajectory equivalence, including mass flux and all SST-owned fields.

Evidence: PR #1017 merged the affine-field boundary reconstruction test; PR #1018 merged the restart-boundary fix. The production RANS selector remains fail-closed because this work repairs restart reconstruction only; it does not integrate SST into the production outer iteration.


### T1 implementation increment — transport-field shape validation (PR #1020 merged)

- [x] Validate scalar dimension as well as cell/face count for `k`, `omega`, strain rate, wall distance, `nu_t`, and authoritative face mass flux before a transactional transport update.
- [x] Add regressions rejecting vector-valued `k` and vector-valued face mass flux.
- [x] Exact-HEAD CI passed on PR head `4aaa801a8bbf1133b5c608c2bf13f72aa647a932`: Release, DebugSanitizers, GUI/headless, production solver E2E, numerical maturity suite, N10 public fixtures, audit gate, and all applicable N-jobs succeeded. Path-filtered long validation and total-validation workflows were skipped and are not claimed as evidence.
- [x] Merge PR #1020 as `a878408494b5b28d22251adb53adfbf710b8aa65`.

**Decision:** the helper now rejects malformed field shapes at its API boundary. This is input-contract hardening only; it does not integrate SST transport into `solve_steady_incompressible` or enable the production RANS selector.

## Next increment — production E2E and case/restart lifecycle

PR #1032 adds the opt-in production outer-loop coupling: it refreshes mesh/carrier-derived state, advances SST transport using the current authoritative mass flux, restores SST on retry/rejection, gates convergence on transport status and boundedness, and refreshes `mu_eff` for the next momentum predictor. Remaining work:

1. Add production E2E coverage proving transport executes, viscosity refreshes, and coupled fields restore on a real solver retry/rejection.
2. Construct `SSTProductionRuntime` from persisted case boundary metadata and initialized or DAT-restored state in the application driver; resolve geometric `WALL_OMEGA_SST` face values explicitly and fail closed when unsupported.
3. Prove continuous-versus-restarted production trajectory equivalence, including authoritative mass flux and all SST-owned fields.
4. Complete transport MMS on at least three mesh levels, frozen DNS channel qualification and mesh/y+ sensitivity, then flat-plate and NACA0012 validation.
5. Keep the general RANS/SST selector fail-closed until physical qualification gates pass.

