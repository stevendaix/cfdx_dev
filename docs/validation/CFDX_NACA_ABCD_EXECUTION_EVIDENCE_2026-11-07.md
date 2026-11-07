# CFDX - NACA XY, coupled operators, aerodynamics, SST and thermal follow-up (2026-11-07)

**Source baseline:** CFDX V5 **after** the 7-file NACA XY/SST separation patch of 2026-11-07. This patch is NOT automatically applicable to arbitrary GitHub PR HEADs. All changes in this follow-up are additive or opt-in; the historic NACA grid and numeric tolerance settings remain unchanged.

## Important separation of cases

- `test_naca0012_geometry.cpp`: planar XY extruded, one Z cell, `EMPTY` front/back, *laminar* reference at Re=1000, `alpha=0`, finite trailing-edge thickness and a 20c outer O-grid. N64/N128/N256 are a numerical qualification family, **not** the NASA SST airfoil.
- `TUR-NACA0012` in `docs/validation/TUR-NACA0012_QUALIFICATION_CASE.md`: existing dedicated turbulent RANS case, PLANNED/BLOCKED without physically converged full SST campaigns.
- NASA TMR `2DN00`/`naca0012numerics_val`: modified NACA geometry and farfield near 500c, nominal Re=6e6 / Mach=0.15 with specified angle(s) and SST formulation; `scripts/import_tmr_naca0012.py` exists but import alone is not SST validation. The TMR reference site moved in February 2026 to https://tmbwg.github.io/turbmodels/ . See https://tmbwg.github.io/turbmodels/naca0012numerics_val.html and https://tmbwg.github.io/turbmodels/naca0012_val.html . The experimental reference is not a measurement of every wall shear curve; distinguish experimental and cross-code comparisons.
- Thermal/CHT, natural convection and P1/S2S/DOM are **separate** physics and benchmark families. A planar NACA `EMPTY` cap must not be silently replaced by WALL in a thermal driver.

## Lot A - discrete truth, conservative flux and checkpoint ownership

**Implemented:**

1. `IncompressibleSolverControls::committed_iteration_callback` is invoked only **after** `transaction.commit()` and after publication of `authoritative_mass_flux`. The original `iteration_output_callback` is deliberately preserved for instrumentation/intentional rejection in the existing test suite; it remains a *tentative* callback and MUST NOT persist restart state.
2. NACA audit rebuilds per-cell `sum(phi_owner - phi_neighbour)` independently from the authoritative face-flux vector, evaluates boundary net and owner/neighbour cancellation. It does **not** reconstruct Rhie-Chow flux from U and does **not** prove `J*v = dR/dx * v`.
3. NACA planar state gate checks `Uz` and flux on EMPTY caps on each *accepted* iteration, not just final output.
4. Checkpoints are saved after a committed iteration, and label accepted-but-NOT-physically-converged states appropriately.

**Remaining P0:** directional finite-difference Jacobian of the **production** FV residual including pressure gradient, Rhie-Chow, deferred nonorthogonal terms, Dirichlet/Neumann/EMPTY boundaries and the correct pressure gauge; perform separate operator/continuity tests for SIMPLE/SIMPLEC/COUPLED. A block-only algebraic certificate is not a proof of this identity.

## Lot B - NACA convergence (no tolerance relaxations)

Actual first accepted iterations (same `alpha_u=0.6`, `alpha_p=0.25`; pressure Krylov is accurate but nonlinear momentum remains far from converged):

| Historical grid | Cells | Momentum relative | Pressure linear residual | Normalized continuity | Max local mass defect | Worst XY momentum cell and defect |
|---|---:|---:|---:|---:|---:|---|
| N64 | 3,096 | 0.00726261 | 1.3681e-11 | 0.000116295 | 0.00508457 | x, cell 191, 0.209823 |
| N128 | 12,336 | 0.0112036 | 4.77048e-11 | 0.000063378 | 0.00332560 | x, cell 386, 0.636310 |
| N256 | 49,248 | 0.0200948 | 3.87961e-11 | 0.000033968 | 0.00184808 | x, cell 770, 1.518707 |

All report `Uz=0` and `EMPTY face flux=0`. All are **NONCONVERGED** because only 1 nonlinear iteration was requested. N64 also ran 10 accepted iterations; an accepted-state checkpoint was created at iteration 10 and correctly marked nonconverged. Do NOT call these numbers a reproduced N256 long-run cause or a proof of N256 convergence.

A script `scripts/analyze_naca_nonlinear_history.py` accepts accepted-iteration CSV, refuses incomplete/nonfinite rows, isolates worst XY momentum component and cell and detects a late momentum plateau without promoting a linear pressure PASS. Next test: 50/100/long N256 iterations with spatially localized face-by-face momentum balance around identified cells; compare pressure traction, tangential viscous stress, Rhie-Chow and nonorthogonal treatment with the same discretization.

## Lot C - grid family and force qualification

The old V5 O-grid uses `r(eta)=eta^1.8`. The *observed* first-layer refinement factors are 3.2344 (N64/N128) and 3.2656 (N128/N256), while tangential/radial indices double. An ordinary one-parameter GCI is therefore **blocked** by the pre-existing `naca_mesh_similarity_gate.py`.

**New opt-in alternative:** `scripts/generate_2d_validation_meshes.py naca0012 --naca-radial-stretch 6.3 --output-dir <dir>`. The positive exponential mapping is `r(eta)=expm1(s*eta)/expm1(s)` (same `s` at all levels). It generates `_es6p3` suffixes and a separate mesh-quality JSON, never overwriting the historical meshes. The family quality gate now verifies that the stretch parameters agree across grids.

| Alternative grid | First layer min | Max 2D nonorth | Inverted cells | First-iteration momentum residual |
|---|---:|---:|---:|---:|
| N64 | 0.00842246 | 54.38 deg | 0 | 0.00455516 |
| N128 | 0.00423676 | 55.36 deg | 0 | 0.00620434 |
| N256 | 0.00218479 | 55.78 deg | 0 | 0.0100114 |

First-layer refinement ratios: 1.988 and 1.939; the existing single-ratio geometry gate accepts this **candidate** family. However its worst nonorthogonality is larger (~55.8 deg vs ~39.4 deg on the baseline), so its robustness/force accuracy is NOT established. The improvement of one-iteration momentum residual is not evidence of nonlinear or spatial convergence.

The updated `scripts/naca_three_grid_gci.py` requires independent quality records, three *physically converged* results, separate finite Cd/Cd-pressure/Cd-viscous/Cl/Cm, and the drag partition. It computes only Richardson **candidates** after these gates; experimental validation is a different criterion. The force normalization by `q*c*span` and moment normalization by `q*c^2*span` and planar checks were already added in the preceding 7-file NACA patch; do not reapply these edits blindly.

**Still necessary:** independent face traction/reaction check and Cp/Cf agreement; same geometry+BC per grid; wall-resolved mesh and 3-grid GCI only on *converged* solutions. Geometry gate is necessary, not sufficient, for the asymptotic GCI regime.

## Lot D - SST, CHT, natural convection, radiation

- SST: `require_sst_empty_wls_dimension()` now refuses opt-in wall WLS gradients with `EMPTY` caps if `sst_gradient_dimension != 2`. WLS-disabled legacy path unchanged; this is a configuration check, not validation of SST equations. The separate `TUR-NACA0012` requires NASA/TMR compatible geometry, Re/Mach/alpha, wall spacing/y+, pressure+shear loads, turbulent transport budgets and model version; **NOT PASS**.
- ChT: local CHT validation targeted **10 passed**; does not establish all CHT interfaces.
- S2S: targeted **8 passed**; P1 and DOM not independently qualified here.
- Vahl Davis: `Ra=1e3` at `4x4`, **PREQUALIFICATION_PASS** with heat-balance ~1.52e-14 and Nu_hot~Nu_cold~1.297805; this is not spatial grid convergence or textbook agreement.
- The user-space Ghia and NACA physical runs still require full CTest DebugSanitizers and long qualification; no claims about the current GitHub HEAD.

## Executed tests (local baseline, not full CI)

- `test_steady_incompressible_solver`: **33/33** including rejected-and-retried attempt, single callback after commit.
- `test_sst_closure_input_contract`: **6/6**.
- `test_naca_planar_contract`: passed.
- `test_naca0012_qualification --quick N64`: passed.
- `pytest` initial GCI/mesh/controlled campaign/force suite: **23 passed**; nonlinear history analyzer **3 passed**; new exponent stretch/mesh/GCI/history set: **18 passed** (overlapping subsets, do not add totals naively).
- `test_cht_validation`: **10/10**; `test_s2s_radiation_vv`: **8/8**; natural convection Ra=1e3 4x4: one **PREQUALIFICATION_PASS**.
- Baseline mesh quality gate returns code 2 (EXPECTED fail-closed) and new exponential grid geometry passes.
- Real N64/N128/N256 1-iteration runs return nonzero (EXPECTED non-convergence, physical gate preserved). One N64 10-iteration checkpoint run also returns nonzero with an accepted state checkpoint at iteration 10.

## Commands

```bash
python3 scripts/generate_2d_validation_meshes.py naca0012 --naca-radial-stretch 6.3 --output-dir build/naca_es63
python3 scripts/naca_mesh_similarity_gate.py build/naca_es63/naca0012_mesh_quality_es6p3.json
python3 scripts/naca_mesh_similarity_gate.py build/validation_meshes/naca0012_mesh_quality.json # expected BLOCKED

# Uses CFDX test_naca0012_qualification compiled from the *patched* sources.
CFDX_NACA_AUDIT_ACCEPTED=1 CFDX_NACA_PHYSICAL_HISTORY_CSV=build/n64_trace.csv \
  build/test_naca0012_qualification --single build/naca_es63/naca0012_laminar_n64_es6p3.h5 0.6 0.25 1
# expected nonzero with 1 iteration; examine the complete trace rather than marking success
python3 scripts/analyze_naca_nonlinear_history.py build/n64_trace.csv --output build/n64_diagnostic.json

# GCI is conditional on a complete converged 3-grid summary; no fabricated values.
python3 scripts/naca_three_grid_gci.py build/naca_summary.json \
  --mesh-quality build/naca_es63/naca0012_mesh_quality_es6p3.json \
  --output build/naca_gci.json
```

## Promotion gate

A: no `Jacobian` claim until directional check against independent production FV residual; no conservation PASS from face topology alone. B: every grid meets nonlinear, face-local conservation, pressure, forces, and planar gates; N256 cannot be promoted after one or 100 iterations solely on linear convergence. C: three solutions on a similar verified family in the asymptotic spatial range; Cp/Cf separated; GCI + independent reference comparison. D: real SST/CHT/radiation benchmarks with saved data and reproducible conditions; a test registry labelled READY does not equal PASS. Zero relaxation of tolerances.
