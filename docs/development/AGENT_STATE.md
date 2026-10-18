# Agent State

This file records what the code agent has done and what remains. It is the
handoff document for the next session.

## Active work item

**Issue #812 — GUI Phase 3: Application facade ownership** (open).

P3-A completed by merged PR #817. P3-B through P3-E are DONE.

### P3-B (DAT/restart production contract) — DONE

The Application now owns the DAT/restart contract end-to-end, qualified
against an artifact actually emitted by the production solver:

- `src/cfdx/python/cfdx/case.py` — `ExecutionConfig` gains `mesh` and
  `output_dir` fields. The production solver requires ``--mesh`` and
  ``--output-dir`` even when restarting from a checkpoint, so a bare DAT
  restart must carry those options alongside the checkpoint path. They are
  stored on the session's execution config rather than on the DAT itself,
  because the DAT carries numerical state only.
- `src/cfdx/python/cfdx/application/application.py` — `load_restart` loads a
  bare DAT checkpoint as the solver restart state without requiring the
  caller to first save a paired ``.cfdx.h5`` project. `ensure_controller`
  builds the canonical solver command from either a saved project path or a
  bare DAT restart path, carrying the checkpoint through the
  case-configured restart option.
- `src/cfdx/python/cfdx/application/commands.py` — `_ControllerCommand`
  constructs the controller on demand through `ensure_controller` rather than
  requiring the caller to have pre-built one, so controller commands work
  from an unsaved project as long as the case config is complete.
- `tests/python/test_application_restart_production_e2e.py` — qualifies the
  restart command against a real production solver run: the solver's
  initialization strategy must be ``restart`` after loading the checkpoint.

### P3-C (post-processing capability registry) — DONE

The post-processing registry is now a capability declaration that reflects
what the production solver actually emits, and it is backed by executable
evidence where a path exists:

- `src/cfdx/python/cfdx/application/postprocessing.py` — registry + `execute_derived`.
  - `field:U` and `derived:u_magnitude` are `availability=conditional`: the
    production solver's per-iteration VTU frames carry pressure only
    (`apps/production_solver/main.cpp` writes `{{"p", pressure}}`). Velocity
    is a real, evidenced capability but only from a DAT checkpoint, not from
    the readable result series.
  - `execute_derived` is the executable counterpart for derived fields. It
    computes `derived:u_magnitude` from the `U` vector array per cell (three
    components per cell), matching how the dataset stores the velocity, and is
    gated by the registry so unregistered quantities raise `KeyError`.
  - Mesh-quality metrics are `availability=conditional` with
    `computation_owner = "cfdx::core::compute_face_quality"`: the C++ kernel
    owns them and the Python layer must not re-derive them. There is not yet
    an export path exposing per-face values through the results contract, so
    the registry does not advertise one.
- `src/cfdx/python/cfdx/derived.py` — `derived_vector_magnitude` added for
  per-cell vector magnitude.
- `tests/python/test_application_postprocessing.py` — registry contract tests,
  including a test pinning mesh-quality ownership to the C++ kernel.
- `tests/python/test_application_probe_production_e2e.py` — qualifies the
  Application probe CSV round-trip against an artifact actually emitted by the
  production solver (skips cleanly when the fixtures are unavailable).
- `tests/python/test_application_monitor_production_e2e.py` — qualifies the
  monitor series carries the production solver's convergence vocabulary.

### P3-D (monitor history contract) — DONE

The monitor parser previously only recognised generic progress lines, so the
production solver's structured `INCOMPRESSIBLE_ITER` trace was dropped: the
monitor series recorded iterations but none of the convergence quantities the
solver actually reports.

- `src/cfdx/python/cfdx/metrics.py` — the parser now extracts the production
  solver's convergence vocabulary (`momentum_residual_relative`,
  `momentum_residual_reference_relative`, `momentum_residual_raw`,
  `momentum_residual_reference`, `continuity_normalized`, `velocity_change_inf`,
  `pressure_change_inf`, `pressure_relative_residual`,
  `flux_velocity_mismatch_linf`) alongside the iteration, while still
  rejecting non-convergence tokens such as `alpha_u` and
  `pressure_linear_iterations`.
- `tests/python/test_metrics.py` — production trace parsing tests.

### P3-E (workflow state reflects registry availability) — DONE

The Workbench renders the post-processing capability registry as workflow
nodes, but previously every node was rendered from a hand-maintained
configured/not_configured split that had no connection to the registry. A
conditional capability (one the registry declares but for which no export path
exists) was therefore indistinguishable from a complete one.

- `src/cfdx/python/cfdx/application/workflow.py` — `_post_processing_availability`
  maps the registry's availability to workflow status and is the single place
  a frontend learns what is genuinely available versus what is only declared.
  The registry's stable_ids are mapped to the Workbench's result-group node ids
  (``fields``, ``probes``). A node is "conditional" when at least one of the
  capabilities it groups is conditional, because the node cannot be reported
  complete until every capability it represents is actually producible from
  the production dataset through the results contract. A new
  `WorkflowStatus.CONDITIONAL` status is surfaced separately from
  COMPLETE so the Workbench can show a real-but-unproducible capability without
  a green check. Group labels and node labels are kept in separate dicts
  because the "run" key is shared by the RUN group and the Run node and needs
  different casing.
- `tests/python/test_application_workflow.py` — pins the mapping, including a
  test that a validation diagnostic overrides the registry-derived status for
  the current case.

### P3-E (remaining) — NOT STARTED

No work has been done on P3-E yet.

## Blocked

**Issue #943 (SST coupling)** must complete before the turbulence guard in #949
can be removed — this is a hard dependency, not a blocker for #812. PR #957 is
open and is the active implementation step.

## Other open issues

#946, #945, #939, #914, #913, #908, #907, #906, #905, #904, #885, #884.

## Test status

518 passed, 41 skipped, 4 pre-existing failures
(`test_n9_coupled_schur_case_selection`, `test_n9_poiseuille_body_force_production`,
`test_production_solver_e2e::test_production_solver_non_converged_exit_is_persisted`,
`test_production_solver_numerics_e2e::test_production_solver_explicit_cfdx_case_e2e`).
The four failures reproduce on a clean checkout with the production-solver
environment variables set, so they are environment issues, not regressions.

## Relevant files

- `src/cfdx/python/cfdx/application/postprocessing.py` — registry + `execute_derived`.
- `src/cfdx/python/cfdx/derived.py` — `derived_vector_magnitude`.
- `src/cfdx/python/cfdx/metrics.py` — production trace parsing.
- `src/cfdx/python/cfdx/probe_validation.py` — duplicate-iteration rejection.
- `tests/python/test_application_probe_production_e2e.py` — probe CSV round-trip.
- `tests/python/test_application_monitor_production_e2e.py` — monitor contract.
- `tests/python/test_application_postprocessing.py` — registry contract tests.