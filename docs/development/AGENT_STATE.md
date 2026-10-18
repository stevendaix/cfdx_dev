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

### Step 1 — Validate the #812 test base — IN PROGRESS

Four tests fail on a clean checkout with the production solver environment
variables set. The failures are reproducible and have been triaged:

| Test | Root cause | Status |
| --- | --- | --- |
| `test_n9_coupled_schur_case_selection` | Python-written `case.cfdx.h5` uses HDF5 2.0 object-header version; C++ solver (HDF5 1.8.7) rejects it | **Real defect** — h5py `copy()` rewrites the file to version 1 and the solver accepts it |
| `test_n9_poiseuille_body_force_production` | Same HDF5 format version mismatch in the test's own `_write_case` helper | **Real defect** — h5py `copy()` rewrites the file to version 1 and the solver accepts it |
| `test_production_solver_explicit_cfdx_case_e2e` | Same HDF5 format version mismatch **plus** a contract drift: the test asserts `("convection", "numerics.convection.upwind")` but the explicit SU2 config now resolves to `("temporal", "numerics.temporal.euler_implicit")` | **Two defects** — HDF5 format + stale assertion |
| `test_production_solver_non_converged_exit_is_persisted` | Test asserts `returncode == 1` for a non-converged run, but the solver returns `0` with `converged: true` | **Solver behavior change** — the channel mesh converges in one iteration |

The production solver binary at `build/cfdx_production_solver` is linked
against `libhdf5_serial.so.103` (HDF5 1.8.7, Debian 11). Files written by
h5py 3.16 (HDF5 2.0) use object-header version 8, which the 1.8.7 reader
rejects as `invalid HDF5 format version`. The validation meshes in
`build/validation_meshes/` are written by
`scripts/generate_2d_validation_meshes.py` with version 1 and are accepted.

### Step 2 — Qualify the full GUI workflow — IN PROGRESS

A new test, `tests/python/test_application_production_workflow_e2e.py`,
drives the real production solver through the complete cycle:

    Open case -> Configure -> Validate -> Run -> Monitor -> Checkpoint
    -> Reload -> Restart -> Results

It verifies the artifact contract (execution.json, convergence.json,
restart.dat), the restart path through `--restart`, the Application facade
controller construction, and the monitor history parser's ability to ingest
the solver's own `INCOMPRESSIBLE_ITER` stdout. Both tests pass.

### Step 3 — Post-processing contract — NOT STARTED

The registry already marks `field:U` and `derived:u_magnitude` as
`conditional` because the production solver's VTU frames carry pressure
only. The remaining work is to make the GUI surface this limitation
explicitly rather than presenting unavailable fields as available.

### Step 4 — Monitoring and status — NOT STARTED

The GUI must distinguish at minimum: running, converged, stopped without
convergence, solver error, and conditional/unavailable. The convergence
history parser already captures the production solver's vocabulary; the
remaining work is wiring the status display to the evidence rather than
to the presence of a process or file.

### Step 5 — Consolidate Application / Workbench — NOT STARTED

Once the contracts are validated, route all business logic through the
`Application` API and its shared state rather than direct HDF5 file
access or Workbench-specific paths.

## Blocked dependency

Issue #943 (SST coupling) remains a dependency for removing the turbulence
guard in #949. PR #957 ("make SST transport updates transactional") was
merged to `origin/master` on 2026-10-18 and is now in the working tree.

## Pre-existing failures (unchanged, environment-related)

- `test_n9_coupled_schur_case_selection` — see Step 1
- `test_n9_poiseuille_body_force_production` — see Step 1
- `test_production_solver_e2e::test_production_solver_non_converged_exit_is_persisted` — see Step 1
- `test_production_solver_numerics_e2e::test_production_solver_explicit_cfdx_case_e2e` — see Step 1

## Environment

- Python: `/usr/bin/python3` (3.12)
- Ruff: 0.13.2 via `python3 -m ruff`
- Build artifacts: `build/cfdx_production_solver`, `build/validation_meshes/channel_n16.h5`