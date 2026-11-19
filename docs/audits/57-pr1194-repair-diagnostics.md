# PR #1194 — repair and temporary diagnostic campaign

## Purpose

This follow-up PR is based on PR #1194 and has two deliberately separate goals:

1. Repair integration defects that prevent the existing regression and SST BFS case from running.
2. Run a bounded GitHub Actions campaign and retain the complete logs, generated cases, solver outputs, residual histories, JSON evidence, build metadata and CTest output as downloadable artifacts.

The workstation result mentioned in the parent PR discussion (converged k-omega SST solution, reported 6–7% physical discrepancy and long runtime) is not available in this repository snapshot. This PR does not claim to reproduce it. The CI campaign uses the repository's generated TH03 debug case and is explicitly diagnostic, not experimental validation.

## Repairs

- Update ScalarEquation aggregate initializers in the TVD deferred-defect regression test to account for the implicit_relaxation boolean added to the struct.
- Serialize sst_adaptive_relaxation as the boolean string true, because the production solver's raw-setting contract accepts only true or false. The previous generator emitted enabled, causing the SST case to fail before entering the solver.

No physical convergence tolerance, residual threshold, iteration acceptance criterion or physical reference value is relaxed.

## CI campaign

The temporary workflow .github/workflows/bfs-pr1194-repair-diagnostics.yml:

- builds the production solver and the repaired regression target in Release with HDF5;
- runs the focused TVD regression and BFS smoke test;
- runs a same-mesh comparison of SIMPLE, SIMPLEC, PCD-auto, BlockLocal-AMG, LSC-AMG and BFBT-AMG;
- runs the PR-specific --bfs-converged, --bfs-fieldsplit and --bfs-pimple modes on a common generated case;
- applies a 600-second wall-clock cap per solver execution, retaining timeout exit codes and partial output;
- uploads the complete diagnostic directory as a 30-day artifact, even when solver runs fail or time out.

The workflow is intentionally diagnostic: solver non-convergence and timeouts are recorded as negative evidence and do not get silently converted into success. Build/test failures and solver outcomes are distinguishable in the uploaded evidence.

## Interpretation rules

- EXIT_ZERO is process status, not a physical convergence certificate.
- The generated debug mesh is not a Driver–Seegmiller validation mesh.
- A small mesh can diagnose configuration, dispatch, stability and cost, but cannot establish grid-independent separation or reattachment.
- No comparison to the unavailable workstation result is valid until the original case, mesh, inlet profile, initial fields, SST settings and solver controls are recovered or checked into an immutable reproducibility package.
- Keep all failed variants and timeouts in the report. Do not select only the best run.

## Artifact contents

The artifact contains the exact commit and runner identity, CMake configuration, build log, focused CTest log and exit code, each campaign log and output directory, generated HDF5 cases, solver run directories, direct-mode logs/exit codes, CMake cache, compile commands and CTest temporary output.

## Removal

This workflow is temporary. After the diagnostic campaign has been reviewed and the required evidence has been transferred into a permanent audit report, remove the temporary workflow or replace it with a smaller permanent regression.
