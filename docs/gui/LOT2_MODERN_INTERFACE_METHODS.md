# CFDX Workbench: interface and modern-method diagnostics (lot 2)

**Base:** `CFDX_Workbench_moderne_2026-11-08.zip` (previous Workbench lot).
**Scope:** case editing, cockpit, execution evidence, and method discovery.
**Do not infer solver qualification from this delivery.** This is a Python/Qt
Workbench improvement, not a new native Newton-Krylov or Anderson backend.

## What is implemented

| Component | Change | Important guarantee |
| --- | --- | --- |
| `boundary_editor.py` | Typed, lossless boundary-condition editing; inlet/outlet/wall/symmetry presets | Selecting a different field no longer silently discards the previous condition |
| `boundary_setup.py` | Validates fields, condition types, finite vector/scalar values | Staged multi-field patches are blocked from production execution until their native mapping is proved |
| `setup_panel.py` | Searchable patch list, contextual controls, presets, explicit destructive-change confirmation | Presets do not silently commit case changes; selecting an empty field resets previous values |
| `convergence_diagnostics.py` | Robust Theil-Sen trend on a bounded sample window | No trend can set a scientific convergence verdict; missing samples are never invented |
| `run_center_panel.py` | Grouped metrics, log-scale plot when pyqtgraph is installed, tabs for logs, batched UI refresh | True solver convergence evidence stays separate from graph trends and process state |
| `execution.py` | Captures boolean `converged` from the solver's `convergence.json` when available | GUI shows `NOT VERIFIED` if the convergence artifact is missing or unreadable |
| `modern_methods.py` and `methods_panel.py` | Read-only catalog of 11 numerical techniques | Adapter mappings and research candidates never imply native integration or success |
| `workbench.py` | Docked Numerical Methods Explorer | Uses the existing application and session layer |
| `prototypes/gui-nextgen/modern_lot2.html` | Self-contained French interaction preview | An optional local JSON load never launches a solver or fabricates physical data |
| `.github/workflows/cfdx-gui.yml` | Adds tests to the Qt/Xvfb regression stage | GitHub CI is **not** run by this local delivery |

## Numerical interpretation

The monitoring panel separates named quantities into:

1. nonlinear/equation residuals;
2. conservation/mass and face-flux mismatch;
3. linear solver residuals;
4. solution updates;
5. other solver-reported metrics.

The slope is a **trend estimate** in log10-decades per iteration, computed
with the median of pairwise slopes (Theil-Sen) on a bounded recent window.
It is only a visual clue for `DECREASING`, `STAGNATING`, `GROWING`,
`OSCILLATING`, and `INSUFFICIENT`. It is **not** a scientific convergence
criterion and does not relax any solver threshold or change the CFD algorithm.

The method catalog documents SIMPLE/SIMPLEC/COUPLED, PCD, AMG, FGMRES,
WLS, TVD, Newton-Krylov, Anderson, and nonlinear GMRES. The latter three
are explicitly marked **research candidates without certified native integration**.

## Known limitations and follow-up

- Qt, PyVista, and pyqtgraph were unavailable in the local verification
  environment. The graphical test is therefore skipped locally; CI on
  a machine with a working Qt/OpenGL stack is still required.
- The existing CFDX single-field boundary schema cannot translate all
  combinations of fields and conditions. This lot stores extra fields
  **losslessly but blocks the run** until a native translation is qualified.
- The session layer still has a *legacy* fallback that calls a clean solver
  exit `CONVERGED` when `convergence.json` is absent. The new cockpit displays
  `NOT VERIFIED` for that situation. Replacing that legacy state in the
  overall application requires a broader compatibility migration.
- Native mesh rendering at 1/10/50 million cells, streaming, real 3D
  comparisons, performance budgets, and solver-level nonlinear acceleration
  are not qualified by this change.
- The HTML prototype is a separate UX preview. Any visual shapes or UI labels
  in it are illustrative, not actual CFDX physics.

## Local verification evidence (2026-11-08)

- New / setup / execution / Workbench targeted suite: **64 passed, 19 skipped**.
- Additional application / visualisation / mesh / results suite:
  **209 passed, 7 skipped**.
- Python `compileall` on modified modules: passed.
- JavaScript syntax check with Node.js `new Function(...)`: passed.
- Patch application on an untouched copy of the **previous Workbench lot**:
  verified separately before delivery.

Tests were run with `PYTHONPATH=src/cfdx/python` from the repository root.
A pass here does not certify C++ production convergence or physical results.

## Reproduce and inspect

Install optional GUI extras using instructions in `docs/gui.md` and launch
`cfdx-gui`. In the Workbench, open **Case Setup** to edit a boundary, **Run**
to inspect residual telemetry, or **Numerical Methods Explorer** to view
read-only method guidance.

For the self-contained prototype, open
`prototypes/gui-nextgen/modern_lot2.html` in a web browser.

## Verified technical references

- Qt: Model/View Programming - https://doc.qt.io/qt-6/model-view-programming.html
- PyVista: Using Common Filters - https://docs.pyvista.org/examples/01-filter/using_filters.html
- PETSc: SNES nonlinear solvers / nonlinear GMRES / Anderson - https://petsc.org/release/manual/snes/
- PETSc: KSP and preconditioners - https://petsc.org/release/manual/ksp/
- hypre: BoomerAMG design - https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html
- ParaView: parallel / remote visualization - https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html

## Next engineering gates

1. Execute real GUI tests under a supported Qt/OpenGL stack, including
   keyboard navigation, focus, screen reader metadata, snapshot changes,
   and several successive case switches.
2. Agree on a native boundary condition schema and lossless mapping for
   (velocity, pressure, temperature, k, omega, epsilon) with conservative
   and dimensional checks for every patch type. Do not unblock the staged
   multi-field cases before this integration passes.
3. Introduce a distinct *completed but unverified* process/session state
   after updating legacy clients and regression expectations.
4. Wire true nonlinear operator diagnostics, conservation and independent
   reference-error tests to the qualification center, not to UI estimates.
5. Benchmark bounded HDF5, chunking, streamed surfaces and VTK interactivity
   on realistic meshes with measured RAM, FPS and render latency.
6. Only then trial Anderson/NGMRES/Newton-Krylov against a native full nonlinear
   residual on N64/N128/N256 and other cross-domain qualification cases.
