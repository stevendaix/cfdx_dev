# CFDX GUI — Phase 2 Architecture and Delivery Plan

Status: implementation baseline for the post-Workbench phase  
Parent: `cfdx/gui`  
Related: #168, #426, #518  
Supersedes the sequencing in `docs/gui/ARCHITECTURE.md`; it does not replace its architectural boundaries.

## 1. Purpose

The original Workbench plan (GUI-0, PR A-I) has been substantially implemented. The next phase must not add isolated widgets or duplicate CFD semantics. It must consolidate the existing GUI into stable application contracts and complete the professional CFD workflow.

The governing rule remains:

> The GUI is a client of the application model. It is never a second CFD model.

The same application contracts must remain usable by GUI, TUI, CLI and automation.

## 2. Current baseline

The repository already contains the main Workbench layers:

- shared `CFDXSession` and application state;
- explicit application commands/events;
- Workbench shell and dock management;
- setup/property/validation panels;
- run center and monitoring;
- results browser and renderer adapters;
- 3-D mesh/result interaction;
- probes and probe validation/export;
- case + DAT persistence and restart workflows;
- diagnostics and multi-tab setup diagnostics;
- optional GUI packaging and Xvfb-based GUI CI.

PR #665 is the current GUI packaging/CI baseline. The remaining work is therefore **contract hardening and workflow qualification**, not another GUI rewrite.

## 3. Target ownership model

```text
                         GUI / TUI / CLI
                               |
                               v
                    Application facade
                               |
             +-----------------+-----------------+
             |                 |                 |
          Project            Setup           Execution
             |                 |                 |
             +-----------------+-----------------+
                               |
                         CFDXSession
                               |
          +------------------+------------------+
          |                  |                  |
       Persistence       Solver/runner       Results
          |                  |                  |
     HDF5 / DAT             C++              VTU
                               |
                       Post-processing
```

### Non-negotiable rules

1. Physics, numerical schemes, restart semantics and validation contracts live below the GUI.
2. GUI controls invoke application commands; they do not mutate solver objects directly.
3. GUI availability is derived from backend/application capabilities.
4. Unsupported capabilities are explicit; no silent fallback or fake GUI state.
5. HDF5 case configuration and numerical DAT restart state remain distinct.
6. Derived visualisation quantities are never presented as authoritative solver state.
7. TUI/CLI and GUI must observe identical case semantics.

## 4. Phase 2 domains

### P2.1 Project and case lifecycle

Establish one coherent lifecycle:

```text
Project
  ├── Case configuration (.cfdx.h5)
  ├── Mesh / topology references
  ├── Physics / numerics configuration
  ├── Results / VTU
  ├── Restart state (.dat.h5)
  └── provenance / compatibility metadata
```

Required contracts:

- open/create/save project;
- case identity and revision;
- dirty/clean state;
- compatibility validation;
- paired DAT handling;
- checkpoint/restart state;
- provenance;
- explicit errors for incompatible files.

Do not introduce a GUI-only project manifest until the repository-wide project/reproducibility contract is defined.

### P2.2 Setup framework

Move toward schema-driven setup rather than widget-driven semantics:

```text
Case schema
    -> property metadata
    -> generic editor
    -> command
    -> validation
    -> ApplicationState
    -> GUI/TUI/CLI
```

The setup domains are:

- general;
- geometry/mesh;
- physics;
- materials;
- boundary conditions;
- initialization;
- numerics;
- solver;
- convergence;
- time integration;
- advanced options.

The GUI must not expose a solver model merely because a kernel or enum exists. Capability exposure must follow the actual supported production path.

### P2.3 Execution and monitoring

Keep the Run Center thin and contract-driven:

- run/pause/resume/stop;
- iteration/time state;
- residuals;
- CFL/time-step information;
- monitor series;
- logs;
- diagnostics;
- checkpoint state;
- restart status.

The production STOP -> checkpoint -> reload -> restart path is an acceptance item, not merely a GUI feature.

### P2.4 Results and display-object model

Formalize visualization as data, not widget callbacks:

```text
ResultDataset
  └── DisplayObject
       ├── surface
       ├── contour
       ├── slice
       ├── vector/glyph
       ├── streamline
       ├── volume
       └── probe
```

A display object should own only presentation parameters and references to authoritative result data. It must be serializable/reproducible without embedding Qt objects.

### P2.5 Post-processing registry

Create one backend/application registry for:

- native fields;
- derived fields;
- mesh-quality metrics;
- point probes;
- surface/volume integrations;
- engineering quantities;
- reports.

The registry should expose:

- stable identifier;
- human-readable name;
- dimensional/unit information;
- input requirements;
- availability/capability;
- computation owner;
- validation status.

The GUI should enumerate the registry instead of maintaining independent lists of derived quantities.

### P2.6 Diagnostics and validation

Use one validation/diagnostic model across setup, execution and results:

```text
Diagnostic
  ├── severity
  ├── code
  ├── domain
  ├── message
  ├── source
  ├── affected object
  └── remediation
```

Distinguish:

- configuration error;
- unsupported capability;
- runtime error;
- numerical warning;
- validation failure;
- informational diagnostic.

Do not turn numerical validation into a GUI-only pass/fail indicator.

### P2.7 Renderer boundary

The renderer remains an adapter:

```text
Application result/display model
            |
            v
       Renderer adapter
            |
       PyVista / VTK
```

No Qt widget should contain CFD-specific mesh, field or solver algorithms.

### P2.8 Front-end parity

The application layer must be independently usable from:

- GUI;
- TUI;
- CLI;
- automation/agents.

A feature is not architecturally complete if it can only be driven through a GUI widget.

## 5. Delivery sequence

The next implementation work should follow vertical slices, not another large migration.

| Slice | Scope | Exit criterion |
| --- | --- | --- |
| P2-A | Architecture/contract audit | authoritative ownership and gaps documented |
| P2-B | Project/case lifecycle | open/save/compatibility/restart contract tested |
| P2-C | Setup schema/property model | representative setup domains share one contract |
| P2-D | Execution/restart E2E | real solver STOP/checkpoint/reload/restart proven |
| P2-E | Results/display-object model | display objects independent of Qt widgets |
| P2-F | Post-processing registry | probes + derived/quality quantities use one registry |
| P2-G | Diagnostics unification | setup/run/results share diagnostic contract |
| P2-H | Renderer hardening | headless contract tests + controlled OpenGL smoke |
| P2-I | GUI production E2E | real case from setup through results/restart |
| P2-J | TUI/CLI parity | representative workflows execute without GUI |

Each slice should be a small PR with focused tests. Avoid parallel PRs modifying the same Workbench core unless they are explicitly stacked.

## 6. Validation matrix

GUI validation must be layered.

### Layer 0 — headless application

No Qt/PySide6/PyVista dependency:

- state;
- commands;
- events;
- setup validation;
- persistence contracts;
- restart compatibility;
- post-processing contracts.

### Layer 1 — Qt Workbench

Under Xvfb:

- shell startup;
- navigation;
- setup;
- property editing;
- diagnostics;
- run controls;
- results selection.

### Layer 2 — renderer

Controlled Qt/VTK environment:

- renderer creation;
- mesh display;
- field display;
- selection;
- display-object updates.

### Layer 3 — production E2E

Real CFDX workflow:

```text
create/open case
  -> setup
  -> validate
  -> run
  -> monitor
  -> checkpoint
  -> stop
  -> reload
  -> restart
  -> write results
  -> inspect results
  -> export
```

A Layer 1 green CI job must never be described as production workflow qualification.

## 7. Immediate priorities

Based on the current repository state and #426:

1. finish the numerical DAT/HDF5 stop/restart contract;
2. audit existing post-processing, mesh-quality and derived-field capabilities before adding duplicates;
3. formalize the project/case lifecycle;
4. formalize schema-driven setup/property contracts;
5. create the display-object/post-processing registry;
6. unify diagnostics;
7. add production-solver GUI E2E;
8. prove GUI/TUI/CLI semantic parity.

## 8. What should not be done

- Do not rewrite the Workbench.
- Do not move CFD algorithms into Qt.
- Do not create GUI-specific physics enums or validation rules.
- Do not create a second restart format.
- Do not duplicate probe/derived-field implementations in GUI.
- Do not weaken numerical tests to make GUI CI green.
- Do not treat screenshots or widget smoke tests as solver validation.
- Do not make OpenGL availability a hidden prerequisite for headless application tests.
- Do not create another long-lived GUI integration branch containing unreviewed feature work.

## 9. Tracking policy

#426 remains the professional-workflow integration tracker.

Individual implementation PRs should:

1. reference #426 and the relevant domain issue;
2. identify the authoritative owner;
3. state which layer is changed;
4. include focused automated tests;
5. distinguish implementation, verification and validation;
6. report GUI CI separately from production-solver validation.

This document is the phase-2 sequencing baseline. The original `docs/gui/ARCHITECTURE.md` remains the historical Workbench architecture and migration boundary.

## 10. Definition of Done

The GUI phase is complete only when:

- application semantics are shared by GUI/TUI/CLI;
- project/case/restart lifecycle is coherent;
- setup is contract-driven;
- execution and monitoring are contract-driven;
- results and display objects are independent of Qt;
- post-processing has one authoritative registry;
- diagnostics are shared;
- GUI CI is green;
- renderer validation is controlled;
- production solver E2E is demonstrated;
- restart continuity is quantitatively verified;
- documentation matches implementation.

No checkbox is closed merely because a GUI widget exists.
