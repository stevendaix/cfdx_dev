# SST production-solver implementation tracker

**Status date:** 2026-10-19  
**Overall status:** IN PROGRESS — RANS production selection remains fail-closed.  
**Scope:** Track implementation and qualification evidence for end-to-end SST in CFDX. Update this file in each turbulence PR; do not mark a gate complete without exact-HEAD evidence.

## Current baseline

- [x] Production driver rejects unsupported non-laminar turbulence selections instead of silently solving laminar equations (merged PR #941).
- [x] SST transport component APIs and dynamic blending exist; they are not yet coupled to the production outer loop.
- [x] Cellwise momentum-viscosity assembly overload exists and has uniform-coefficient equivalence coverage.
- [x] Production solver consumes a cellwise turbulent-viscosity field consistently in predictor, final nonlinear residual, and independent momentum conservation reconstruction (implementation in this PR; CI still pending).
- [ ] SST state `k`, `omega`, wall distance and turbulence residuals are owned by the production solve and participate in rollback/convergence.
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
- [ ] Exact-head CI/build/test evidence reviewed.

**Current implementation PR:** [#967 — cellwise effective viscosity](https://github.com/stevendaix/cfdx_dev/pull/967) (draft; exact-HEAD CI pending).  
**Known risk:** this plumbing alone does not enable SST and must not bypass the production RANS guard.

### T1 — SST fields and boundary contract

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
- PR #958 — SST strain-rate utility (open; not production coupling).
- PR #964 — explicit boundary-value Gauss gradient API (open; prerequisite, Dirichlet face-value contract only).
- Issue #943 — end-to-end SST coupling plan (created by the prior work; verify issue state/number if the tracker is moved).

## Evidence log

| Date | Change / evidence | Decision |
|---|---|---|
| 2026-10-19 | Source audit: production solver used one scalar `mu_eff` in predictor, final residual, and independent conservation reconstruction; cellwise assembly overload already exists. | T0 first. Do not enable RANS yet. |
| 2026-10-19 | PR #967 opened. Audit gate first run correctly flagged a missing method-fiche link in the PR description; the description now links `docs/audits/09-turbulence-rans.md`. Re-run the gate on the next exact HEAD. | Keep draft until the gate and all required CI pass. |

## Merge/promotion rules

- No gate is complete merely because a component unit test passes.
- No merge while required CI is pending or unclear.
- Do not widen tolerances, disable tests, silently fall back to laminar, or substitute a laminar case for DNS evidence.
- Keep the production turbulence guard until all required coupling and qualification gates pass.
