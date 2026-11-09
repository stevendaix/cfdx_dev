# CFDX Workbench - Lot 3: OLED black and task-oriented UX

Base: `CFDX_Workbench_interface_methodes_lot2.zip` (user-provided previous lot).
This change modifies the Qt Workbench UI and its read-only method catalog;
it **does not** change the native C++ solver, physical results, convergence
thresholds, or qualification certificates.

## First-run experience and theme rules

- Black is the first-run default, independent of the OS (`#000000` main canvas).
- Graphite and light remain opt-in appearance preferences.
- The chosen theme is stored as `appearance/theme` in Qt `QSettings`;
  missing or unknown preferences fall back to black, without changing case data.
- All three themes share semantic colors for panels, inputs, text, borders,
  selection, status, menus, tables, tabs, scrollbars and focus.
- QPalette disabled roles and focus indicators are explicit. Contrast ratios
  have unit tests across all three themes for text and actionable boundaries.
- VTK scene background, probe and convergence chart chrome follow the active
  theme; scalar ranges, derived fields, numerical inputs and colormaps do not.
- The old standalone `CFDXMainWindow` also uses black by default.

## UI improvements in the actual Qt code

- Reorganized top shell: File / Edit / Visualization / View menus, short
  primary toolbar, Run emphasis, appearance dropdown, and keyboard navigation.
- Five workspaces: Setup, Mesh, Run, Analyze, Validate. Shortcuts `Ctrl+1` to
  `Ctrl+5` switch contexts. Existing docks and object names remain stable.
- An actionable empty-state dashboard guides the user to setup, mesh, run and
  qualification. Buttons connect to the **real** workbench workspace actions.
- `View > Reset panel layout` uses the pre-restoration dock snapshot; the
  action does not reset the case, running solver, or chosen appearance.
- The `Ctrl+K` palette exposes appearance choices and the layout reset, in
  addition to the established capability registry with missing prerequisites.
- Methods Explorer: searched methods, documented mapping-vs-research filter,
  evidence-aware two-method comparison. No read-only exploration sets a solver
  option or asserts physical qualification.
- Interactive HTML design review: `prototypes/gui-nextgen/modern_lot3.html`;
  first-use black, optional graphite/light, keyboard command palette, and
  optional **real** JSON registry/monitor imports. Not a CFD solver client.

## Verification and acceptance

Run the following on the exact patched checkout:

```bash
PYTHONPATH=src/cfdx/python pytest -q \
  tests/python/test_theme.py \
  tests/python/test_modern_methods.py \
  tests/python/test_workbench_modern_gui.py \
  tests/python/test_workbench_lot3_prototype.py
python -m compileall -q src/cfdx/python/cfdx
```

`cfdx-gui.yml` now includes the theme and prototype checks and requests black
and light screenshots, but those CI jobs have **not** run in this local lot.
GUI tests (Qt/OpenGL) are skipped when optional dependencies are not installed.
A JavaScript syntax check cannot replace a rendered browser/interaction test.

## Scientific and architectural boundaries

- No new CFD numerical algorithm is implemented by this UI lot.
- Anderson, NGMRES and Newton-Krylov remain documented research candidates,
  not qualified native production algorithms.
- A process finishing with code zero is not independent nonlinear convergence
  evidence. The existing Run Center still surfaces `NOT VERIFIED` when the
  solver's convergence artifact is missing.
- The legacy cross-application session state fallback identified in lot 2
  remains unresolved; it must be migrated with compatibility tests.
- Real OpenGL screenshot comparison, Qt interaction checks, GPU performance,
  and benchmarked industrial mesh visualization still require a GUI-equipped
  runner and cannot be inferred from headless Python tests.

## Bibliography (authoritative implementation references)

- Qt QMainWindow dock state/toolbar persistence:
  https://doc.qt.io/qt-6/qmainwindow.html
- Qt QSettings durable configuration:
  https://doc.qt.io/qt-6/qsettings.html
- Qt QPalette and active/disabled color groups:
  https://doc.qt.io/qt-6/qpalette.html
- PyVista set_background (render-only, not scalar colormap):
  https://docs.pyvista.org/api/plotting/_autosummary/pyvista.plotter.set_background
- PETSc SNES for future *solver-level* nonlinear acceleration qualification:
  https://petsc.org/release/manual/snes/
- ParaView visualization principles for future large-mesh work:
  https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html
