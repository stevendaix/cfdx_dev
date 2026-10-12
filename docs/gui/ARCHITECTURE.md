# CFDX GUI architecture

Status: proposed implementation baseline  
Parent branch: `cfdx/gui`  
Scope: GUI-0 / PR A — design system and application architecture

## 1. Purpose

CFDX already has a functional optional PySide6 GUI, a shared `CFDXSession`, an
`ExecutionController`, mesh/result adapters, and contracts for case I/O and
validation. This document defines the target architecture for evolving that
GUI into a Workbench without duplicating CFD state or coupling Qt widgets to
numerical implementation details.

The first implementation rule is preservation: existing setup, execution,
monitoring, and 3D behavior remains usable while responsibilities move behind
small application-facing contracts.

## 2. Non-negotiable boundaries

```text
Qt views and actions
        |
        v
Application facade / commands / events / capabilities
        |
        +--> Case model and validation
        +--> ExecutionController and MonitorSeries
        +--> Results and renderer adapters
        |
        v
CFDXSession (single source of truth)
        |
        +--> HDF5/DAT/VTU persistence
        +--> solver runner and C++ numerical kernels
```

The GUI may know about domain concepts such as `Boundary`, `Field`, `Material`,
`Monitor`, `Result`, `MeshEntity`, and `ValidationMessage`. It must not know how
the solver assembles matrices, computes gradients, writes HDF5, or executes
MPI.

### State ownership

| Concern | Authoritative owner | GUI responsibility |
| --- | --- | --- |
| Case data and revisions | `CFDXSession.case` and session revisions | Display and request changes |
| Run state, iteration, time | `CFDXSession` / `ExecutionController` | Render state and enable actions |
| Solver metrics and monitor history | `ExecutionController.monitor_series` | Subscribe and plot |
| Mesh entities and stable IDs | mesh/application adapters | Display selection |
| Results and fields | result discovery/renderer adapters | Request display |
| Validation | shared validation/application layer | Group, filter, and navigate messages |
| Layout and user preferences | GUI persistence layer | Restore safe UI state only |

No widget should become a second model for any row in this table.

## 3. Target package layout

The target layout is intentionally incremental. It does not require moving the
current modules in one large change.

```text
python/cfdx/
├── application/
│   ├── __init__.py
│   ├── application.py       # facade used by GUI, TUI, and future CLI flows
│   ├── state.py             # immutable/read-only application snapshots
│   ├── commands.py          # explicit user operations and undo metadata
│   ├── events.py            # typed state/change notifications
│   ├── validation.py        # application-level validation aggregation
│   ├── capabilities.py      # backend-advertised feature registry
│   ├── project.py           # project/case identity and persistence boundary
│   ├── workflow.py          # setup/run/results steps and status
│   ├── execution.py         # application adapter around ExecutionController
│   └── results.py            # result/display-object application contracts
└── gui/
    ├── main_window.py       # Workbench shell and composition root
    ├── shell/                # toolbar, menus, status bar, persistence
    ├── navigation/           # workflow/case tree and selection model
    ├── setup/                # setup pages and property editors
    ├── run/                  # run center, logs, monitors
    ├── results/              # results tree and display controls
    ├── properties/           # contextual property panels
    ├── viewport/             # renderer adapter host, never CFD logic
    ├── dialogs/               # file/validation/preferences dialogs
    └── widgets/               # reusable presentational widgets
```

During migration, the existing `cfdx.gui` remains the compatibility entry point
for `cfdx-gui`. New code should depend on application contracts rather than
importing the monolithic window implementation.

## 4. Application state and events

`ApplicationState` is a read-only snapshot assembled from the authoritative
session, controller, validation report, selection, capabilities, and project
metadata. It is a view model for rendering, not a replacement for
`CFDXSession`.

The minimum state groups are:

- `project`: path, name, dirty flag, and persistence status;
- `workflow`: ordered steps and `complete/warning/error/not_configured` status;
- `simulation`: session state, iteration, time, revision, restart requirement;
- `selection`: selected stable entity ID and domain kind;
- `execution`: latest metrics, monitor series summary, and error;
- `results`: available frames, fields, and active display object;
- `capabilities`: features exposed by the backend;
- `diagnostics`: blocking errors and non-blocking warnings.

Events are notifications that a snapshot may have changed. They do not carry
mutable widget objects and do not instruct a view to manipulate the solver.
The first event vocabulary is:

- `ApplicationStateChanged`
- `ProjectChanged`
- `SelectionChanged`
- `ValidationChanged`
- `ExecutionChanged`
- `ResultsChanged`
- `CapabilitiesChanged`

Qt signals may adapt these events at the GUI boundary, but the application
contracts remain usable by TUI and CLI clients.

## 5. Commands and undo/redo

User intent is represented by explicit commands. A command validates its input,
invokes the application/session API, and emits the resulting state change.
Widgets only construct commands and render their outcome.

Initial command vocabulary:

- `SetBoundaryValue`
- `SetMaterialProperty`
- `SetPhysicsModel`
- `SetNumericalOption`
- `ValidateCase`
- `InitializeCase`
- `RunSolver`, `PauseSolver`, `ResumeSolver`, `StopSolver`
- `CreateProbe`, `CreateContour`, `CreateSlice`
- `SaveCase`, `SaveCheckpoint`, `LoadCheckpoint`

Commands that mutate case configuration must record the affected revision and
change impact (`HOT`, `RESTART`, or `REBUILD`). Execution commands are not
pretended to be reversible: their history entry records intent and outcome,
while undo/redo applies only to supported configuration commands.

Undo/redo is deliberately a later implementation step. The command boundary
comes first so that adding it does not require another GUI rewrite.

## 6. Capabilities

The GUI must derive actions and pages from capabilities reported by the
backend/application layer. It must never display an action whose backend
contract is unavailable.

Capability names use stable dotted identifiers, for example:

```text
solver.steady
solver.transient
solver.incompressible
solver.energy
solver.turbulence
solver.mpi
initialization.uniform
initialization.field
post.contour
post.vector
post.streamline
post.slice
post.probe
post.derived_field
```

A capability has an identifier, availability, optional reason when unavailable,
and source/version metadata. Unknown capabilities are ignored by older clients;
known unavailable capabilities are rendered disabled with an explanatory
message rather than silently simulated.

## 7. Workbench navigation and layout

The target Workbench is a dockable shell:

```text
+--------------------------------------------------------------------------+
| CFDX | Project | Save | Check | Run | Stop | Search                      |
+--------------+-------------------------------------------+----------------+
| WORKFLOW     |                                           | PROPERTIES     |
| Project      |              3D VIEW / VIEWPORT          | selection      |
| Setup        |              renderer adapter             | contextual     |
| Run          |                                           | fields         |
| Results      +-------------------------------------------+----------------+
|              | MONITORS / CONSOLE / DIAGNOSTICS         |                |
+--------------+-------------------------------------------+----------------+
| Status: Ready | cells | MPI | iteration | residual                    |
+--------------------------------------------------------------------------+
```

The initial shell uses Qt dock widgets and named object IDs so that layouts can
be persisted safely. The layout persistence layer may store panel visibility,
size, and dock placement, but never stores CFD data or solver state. If a saved
layout is incompatible with a future version, CFDX falls back to the default
layout.

The existing conceptual Setup / Run / Results separation remains available as
navigation groups, but the user is not forced through three isolated pages.

## 8. Incremental PR plan

Each PR must be reviewable independently, based on `cfdx/gui`, and keep the
optional GUI importable without PySide6.

| PR | Branch | Deliverable | Verification |
| --- | --- | --- | --- |
| A | `cfdx/gui-architecture` | This architecture and migration contract | Documentation review |
| B | `cfdx/gui-workbench-shell` | Dockable Workbench shell with placeholder panels | Qt smoke test and screenshot/manual review |
| C | `cfdx/gui-application-state` | Read-only state snapshot and event adapter | Headless unit tests |
| D | `cfdx/gui-commands` | Command boundary for setup/run actions | Command and session integration tests |
| E | `cfdx/gui-navigation` | Workflow tree, status badges, stable selection | Qt/model tests |
| F | `cfdx/gui-properties` | Contextual properties with validation routing | Setup/validation integration tests |
| G | `cfdx/gui-run-center` | Monitor/log/metrics dock backed by controller | Execution/Qt integration tests |
| H | `cfdx/gui-results` | Results tree and display-object contracts | Renderer-contract tests |
| I | `cfdx/gui-migration` | Migrate existing widgets and retire tab-only shell | Full GUI regression suite |

The branch `cfdx/gui` is the integration line for these PRs. It is intentionally
kept at a known-good baseline and is not used as a long-lived feature branch
containing unreviewed work.

## 9. Acceptance criteria for the first shell

The next implementation PR (B) is complete when:

1. `cfdx-gui` still starts with the optional GUI dependencies installed;
2. the shell exposes toolbar, workflow tree, viewport, properties, monitor/
   console, and status regions;
3. panels can be hidden, shown, docked, and restored;
4. placeholder actions do not claim unsupported backend capabilities;
5. no existing case/run/results regression is introduced;
6. headless installations can still import `cfdx` without Qt/PyVista.

This keeps the visual Workbench work observable without prematurely moving
physics or persistence logic into Qt.
