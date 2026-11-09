# CFDX Next-generation Workbench — interactive prototype

This folder contains a **visual and interaction prototype**, not a replacement GUI and not a CFDX runtime client.

## Run locally

From the repository root:

```bash
python -m http.server 8000 --directory prototypes/gui-nextgen
```

Open http://localhost:8000.

No npm install, build step, browser extension, or network access is required.

## What to explore

- **Viewport-first workspace:** the simulation is the main surface; controls appear as contextual overlays.
- **Workspace contexts:** Setup, Mesh, Run and Analyze update the active workspace context and viewport title.
- **Movable contextual panels:** drag the workflow and selection inspector by their headers; show/hide them from the workspace controls.
- **Command palette:** press `Ctrl+K` (or click the search field) to filter and execute demo actions: panel visibility, comparison, theme switch, and clearly labelled non-connected case/validation actions.
- **Run cockpit:** residual-history hierarchy, iteration/time and solver status appear together.
- **Comparison layout:** switch to a two-viewport comparison concept.
- **Theme:** switch between graphite/cyan dark and cool neutral light themes; the selected theme is persisted locally in the browser.
- **Run button:** deliberately does not launch a solver.

## Data and integration boundary

All mesh geometry, residual curves, temperatures, cell counts and solver status in this prototype are **illustrative UI data**. They are not read from CFDX and must not be interpreted as solver results. The comparison view is a layout demonstration, not yet a synchronized scientific comparison.

The prototype deliberately does not:
- load or write `.cfdx.h5`, `.dat.h5` or VTU files;
- launch or control the solver;
- implement physics, mesh editing, or numerical post-processing;
- duplicate `CFDXSession`, `Application`, execution state, or capability logic.

The production integration must reuse the existing application facade, session, events, capability registry, renderer adapters, persistence and numerical contracts.

## Review status

This is a visual/interaction proposal for explicit product review. It is not a production GUI, not a qualified workflow and not a migration of the current PySide6 Workbench.

## Lot 2: modern workspaces and real-data import

Open `modern_lot2.html` for the latest interaction concept. It adds
numerical-method discovery, patch editing mockup, Mesh/Run/Analyze/Validate
navigation, dark/light mode and *optional* local import of a monitor or
qualification JSON file. Data and the native solver are not modified.
All content remains a browser preview; see
`docs/gui/LOT2_MODERN_INTERFACE_METHODS.md` for the integrated Python code.

## Lot 3: OLED black by default and accessible command palette

Open `modern_lot3.html`. Its first-use background is `#000000` (OLED black),
with alternative graphite and light themes. Theme choice is saved locally by
browser preference, and the search button / `Ctrl+K` opens a keyboard-friendly
palette. Monitor and registry JSON can still be inspected, without inventing
CFD results. The native Qt Workbench's default theme is also black; unlike this
HTML preview, it runs real CFDX case and solver actions.
