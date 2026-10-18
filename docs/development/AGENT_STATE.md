# Agent State

This file records what the code agent has done and what remains. It is the
handoff document for the next session.

## Active work item

**Issue #812 — GUI Phase 3: Application facade ownership** (open).

P3-A completed by merged PR #817. P3-B through P3-E are DONE.

### P3-B (DAT/restart production contract) — DONE
### P3-C (post-processing capability registry) — DONE
### P3-D (monitor history contract) — DONE
### P3-E (workflow status) — DONE

## Current focus: GUI Phase 3 qualification (issue #812)

### Step 1 — Validate the #812 test base — CLOSED

The HDF5 format-version investigation is resolved. The root cause was **not**
the object-header version byte (it is `0x54` / v1 in both failing and working
files). The real defect was in `write_case_cfdx_h5` (C++): it opened the output
file with `H5F_ACC_RDONLY` instead of `H5F_ACC_TRUNC`, so every Python-written
`case.cfdx.h5` was silently truncated to an empty file on the C++ side. That
empty file then failed `H5Fopen` and surfaced as "CFDX case load failed".

The fix (commit `c37639ac`) changes the open mode to `H5F_ACC_TRUNC`. The
`H5Tis_variable_str` fix in `case_read_attr_str` (commit `729e19a4`) remains
correct and necessary: h5py 3.x writes variable-length C_S1 by default, and
reading those as fixed-size produced garbage `format_version`/`schema_version`
strings.

### Step 2 — Qualify the full GUI workflow — DONE

`tests/python/test_application_production_workflow_e2e.py` drives the real
production solver through the complete cycle: Open case -> Configure ->
Validate -> Run -> Monitor -> Checkpoint -> Reload -> Restart -> Results. Both
tests pass.

### Step 3 — Post-processing contract — DONE

`tests/python/test_application_postprocessing_contract.py` qualifies the
post-processing availability contract against real VTU frames. The registry
correctly marks `field:U`, `derived:u_magnitude`, and `mesh_quality:*` as
`availability="conditional"`.

### Step 4 — Monitoring and status — DONE

`src/cfdx/python/cfdx/execution.py` was rewritten to read the solver's own
`convergence.json` verdict instead of inferring convergence from the process
exit code. `ProcessResult` carries `output_dir`; `SolverRunner.start()` gains
an `output_dir` parameter; `_complete` now distinguishes CONVERGED/STOPPED
based on `converged: True/False` with a fallback to the old exit-code behavior.

### Step 5 — Consolidate Application / Workbench — NOT STARTED

Once the contracts are validated, route all business logic through the
`Application` API and its shared state rather than direct HDF5 file access
or Workbench-specific paths.

### SU2 adapter boundary-patch mismatch — CLOSED

`test_n9_coupled_schur_case_selection` failed with "case has no boundary
condition for mesh patch 'boundary'". The SU2 adapter wrote a single generic
`boundary` patch while the solver matched BCs against named patches (airfoil,
farfield, etc.). Fixed in commit `c37639ac`:

- `_meshio_to_cfdx` now deduplicates faces by canonical vertex set so the CSR
  contract holds (one face entry per canonical face, owner + neighbour).
- 2D edge faces are promoted to non-degenerate 3D triangles by adding a
  synthetic midpoint vertex per edge with a unique tiny z offset.
- The adapter-supplied `boundary_patch_assignment` is resolved onto CSR face
  indices and emitted as named patches (one per marker name); duplicate patch
  names are merged.
- `case_hdf5_io.cpp`: `write_case_cfdx_h5` used `H5F_ACC_RDONLY` instead of
  `H5F_ACC_TRUNC` (the root cause of the "CFDX case load failed" symptom that
  masked the topology bugs above).

All 6 Schur models in `test_n9_coupled_schur_case_selection` now pass.

## Blocked dependency

Issue #943 (SST coupling) remains a dependency for removing the turbulence
guard in #949. PR #957 ("make SST transport updates transactional") was
merged to `origin/master` on 2026-10-18 and is now in the working tree.

## Pre-existing failures (unchanged, environment-related)

| Test | Root cause | Status |
| --- | --- | --- |
| `test_packaging.py::test_package_builds` | Packaging/build artifact check | Pre-existing, not a code defect |
| `test_production_solver_numerics_e2e.py::test_production_solver_explicit_cfdx_case_e2e` | Stale assertion: test expects `("convection", "numerics.convection.upwind")` but the explicit SU2 config resolves to `("temporal", "numerics.temporal.euler_implicit")` | Pre-existing contract drift |
| `test_n9_poiseuille_body_force_production` | Stale assertion | Pre-existing |
| `test_runner.py::test_runner_streams_stdout_and_stderr` | Intermittent failure in full runs due to test-ordering/cleanup; passes when run alone | Pre-existing |

## Environment

- Python: `/usr/bin/python3` (3.12)
- Ruff: 0.13.2 via `python3 -m ruff`
- Build artifacts: `build/cfdx_production_solver`, `build/validation_meshes/channel_n16.h5`