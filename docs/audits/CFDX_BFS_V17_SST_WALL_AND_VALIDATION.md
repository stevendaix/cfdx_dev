# CFDX V17 — BF Step: SST, wall operators and physical truth

**Status:** opt-in experimental implementation; **not** OpenFOAM-equivalent;
**not** converged/qualified on Driver–Seegmiller; no tolerances changed.

## Source lineage and test configurations

- Exact V16 reference: `CFDX_BFS_V16_SOURCES_COMPLETES.zip`.
- Cumulative source reference: user archive `cfdx_dev-master (2)(1).zip`.
- Diagnostic mesh: reconstructed six-block BFS, 712 cells, 2907 faces.
  This is not the native OpenFOAM `blockMesh` mesh nor the NASA validation mesh.
- Both compared HDF5 cases use exactly the same mesh, inlet, wall patch names,
  Reynolds number, SST tolerances and molecular boundary-diffusion setting.
  Only `v17_openfoam_wall_operator_experimental` differs.
- Build: `-DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON`; the executable is
  a diagnostic build rather than a representative optimized-performance binary.

## Physical and numerical changes

1. An explicit opt-in `v17_openfoam_wall_operator_experimental=enabled`.
   Requires steady, non-pseudo-momentum SST, all WALL patches with stationary
   no-slip U, `kqRWallFunction` (zero-gradient k), `omegaWallFunction` and
   `sst_molecular_wall_diffusion=true`.
2. The SAME candidate U, k, omega, nu_t is used to derive wall-face Spalding
   `nut` for momentum diffusion, wall-face omega for the SST scalar diffusion,
   and wall-owner production G for k and omega sources. Owner production and
   omega references use area-weighting at cells sharing multiple wall faces.
3. G is a **kinematic production rate** (m^2/s^3). The transport source is
   rho times this rate; V8 and V15 wall diagnostics now use consistent units.
4. The wall closure is refreshed for every nonlinear SST candidate and
   after the SST transaction. The independent returned-state residual and
   row-source budget use the refreshed wall values, not static pre-SST faces.
5. The original molecular-only wall mode remains available and unchanged by
   default. V17 does **not** claim the exact OpenFOAM wall-cell constraint,
   corner blending, gradient or `nutUSpaldingWallFunction` iteration semantics.
6. A separate accepted-state evidence parser counts actual V13 transaction
   commits and demands a V14 same-operator contraction ledger before claiming
   any three-state contraction. Neither means physical validation.
7. A field audit computes local Pk, Dk, nut/nu and the conversion
   `p_static = p_star - 2*rho*k/3`; raw p_star must NOT be used as NASA Cp.

## Reproducible observations

All source data and complete diagnostics are in the V17 proof archive.

| Diagnostic | Observed |
|---|---|
| Control: molecular-wall SST, no V17 wall functions | First SST iteration converged in 411; one outer state committed, global U merit 1.78341, NOT converged |
| V17 active wall treatment, same mesh/other settings | FGMRES outer linear system converged; SST rejected after 800 internal iterations; no accepted outer state |
| V17 active + SST nonlinear line-search | First SST attempt rejected after 7 internal iterations; no accepted outer state |
| V17 active + SST scalar relaxation 0.3 | Reduced early oscillations but no demonstrated convergence; run was stopped, NOT qualified |
| V16 accepted-field hotspot | cell 426, k=543.52854587 m2/s2, Pk/Dk~3.6809, limiter inactive |
| Modified vs static pressure at cell 426 | Delta Cp ~0.370951, so raw p* comparisons would be wrong |
| 3 committed U,p,k,omega states | NOT demonstrated |
| NASA reference fields, 3-grid spatial study, wall-function equivalence | NOT demonstrated |

The local hotspot source estimates are NOT the actual FV cell budget without
advection, diffusion and volume integration. An exceptionally large k by itself
is a warning to investigate, not a proof of numerical incorrectness.

## Governing numerical acceptance contract

- Keep original SST absolute and independent returned-state residual criteria.
- Reject/rollback all SST failures and never count them as outer acceptance.
- The FGMRES/ILUT small-system diagnostic does not establish scalable
  COUPLED/PCD/AMG/MGR convergence on 20,540 cells or beyond.
- For each Schur variant, report the actual linear plan, KSP true residual,
  four-field nonlinear merit, positivity, conservation and wall operator ID.
- Before any Driver–Seegmiller comparison require fully converged states,
  native reference boundary conditions, and 3-grid checks; only then evaluate
  x_r/h=6.26±0.10, Cf, Cp with correct pressure origin, velocity profiles and
  Reynolds stresses against independent NASA data.

## Commands

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build --target cfdx_production_solver test_sst_wall_operator_v17
ctest --test-dir build --output-on-failure -R test_sst_wall_operator_v17
python3 -m pytest -q tests/python/test_audit_bfs_sst_hotspots_v17.py \
  tests/python/test_audit_bfs_commits_v17.py
python3 scripts/audit_bfs_sst_hotspots_v17.py --vtu result_1.vtu \
  --output hotspots.json
python3 scripts/audit_bfs_commits_v17.py --log bfs_run.log \
  --output committed_history.json
```

NASA data and OpenFOAM case are external standards, not bundled numerical results.
The option is experimental and must never be enabled in routine physically
certified workflows without further independent tests.
