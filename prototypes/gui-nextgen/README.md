# CFDX Next-generation Workbench — interactive prototype

This folder contains an **interaction and visual-design prototype**, not a replacement GUI and not a CFDX runtime client.

## Run locally

From the repository root:

```bash
python -m http.server 8000 --directory prototypes/gui-nextgen
```

Open http://localhost:8000.

No npm install, build step, browser extension, or network access is required.

## What to explore

- **Viewport-first workspace:** the simulation remains the main surface; controls appear as contextual overlays.
- **Workspace switching:** Setup, Mesh, Run, Analyze change the working context without pretending to create or destroy CFD state.
- **Movable panels:** drag panel headers; resize panels from their lower-right corner; hide/show panels from the top command strip.
- **Workflow canvas:** inspect the case as connected stages instead of a permanently dominant tree.
- **Contextual inspector:** selection changes the visible properties and available actions.
- **Command palette:** press `Ctrl+K` (or click the command field) to search demo commands.
- **Run cockpit:** residual history, iteration/time and solver state are shown together.
- **Multi-view comparison:** switch to a two-viewport comparison concept from the top-right control.

## Data and integration boundary

All mesh geometry, residual curves, temperatures, cell counts and solver status in this prototype are **illustrative UI data**. They are not read from CFDX and must not be interpreted as solver results.

The prototype deliberately does not:
- load or write `.cfdx.h5`, `.dat.h5` or VTU files;
- launch or control the solver;
- implement physics, mesh editing, or numerical post-processing;
- duplicate `CFDXSession`, `Application`, execution state, or capability logic.

The production integration must reuse the existing application facade, session, events, capability registry, renderer adapters, persistence and numerical contracts. This prototype is a way to review the interaction model before choosing a production frontend technology.

## Review status

This is a visual/interaction proposal for explicit product review. It is not a production GUI, not a qualified workflow and not a migration of the current PySide6 Workbench.
