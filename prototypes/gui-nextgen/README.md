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
- **Command palette:** press `Ctrl+K` (or click the search field) to filter demo commands.
- **Run cockpit:** residual-history hierarchy, iteration/time and solver status appear together.
- **Comparison layout:** switch to a two-viewport comparison concept.
- **Theme:** switch between graphite/cyan dark and cool neutral light themes.
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
