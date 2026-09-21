# CFDX Fluent-like application workflow

This document defines the first vertical slice of the application architecture.

## User workflow

1. Create/open a case.
2. Inspect geometry and mesh in a common scene model.
3. Configure physics, materials, boundaries and numerics through typed parameters.
4. Validate the case.
5. Run a steady iteration target or transient time target.
6. Pause or stop and inspect monitors/checkpoints.
7. Edit hot parameters directly, or edit restart/rebuild parameters with an explicit warning.
8. Resume/restart from the checkpoint.
9. Register derived fields, monitors and reports.
10. Render the same case through GUI, TUI or automation.

## Current implementation

The application layer currently contains:

- `CaseModel`: shared typed case state and revisions.
- `SimulationController`: lifecycle, run/pause/stop, transient time, edit gating and checkpoints.
- `MonitorManager`: persistent monitor samples.
- `SceneModel`: stable geometry/mesh object identities and selection.
- `PostProcessor`: derived-field and report registration.
- `TuiRenderer`: a dependency-free status view suitable for a terminal frontend.

The GUI/VTK renderer and the concrete CFD solver adapter intentionally remain separate. This keeps the core usable on HPC nodes without a display server and allows a future Qt/VTK frontend to use the same controller API.

## Validation

`tests/unit/test_application_layer.cpp` is a vertical integration test covering:

- transient time advancement;
- pause and stop;
- hot versus rebuild-required edits;
- checkpoint metadata;
- TUI state rendering;
- monitor samples;
- scene selection;
- derived fields and reports.

Build and execute with the normal CMake/CTest workflow:

    cmake -S . -B build -DCFDX_BUILD_TESTS=ON
    cmake --build build --parallel
    ctest --test-dir build --output-on-failure

The test is deliberately solver-independent. Numerical solver validation remains in the existing M1-M4 benchmark/validation suites.


## Complete workflow slice

The application layer covers the full UI-neutral workflow contract:

- case tree with Geometry, Mesh, Physics, Materials, Boundaries, Numerics, Solver, Monitors and Reports nodes;
- property editing routed through the simulation controller and classified as hot/restart/rebuild;
- steady and transient run targets, including run-to-time;
- explicit run/pause/stop state transitions;
- solver adapter lifecycle (validate, begin, iterate, end);
- checkpoint storage through an interface, with memory and portable text implementations;
- command bus for TUI/automation commands (validate, run, run-until, pause, stop, set, checkpoint, restore);
- common scene selection and toolbar state for a future GUI;
- monitor, report and derived-field models shared by GUI/TUI/automation.

This is intentionally the application contract rather than a display toolkit. A Qt/VTK frontend can bind to Workbench, while HPC/headless execution can use the same objects through CommandBus.

### Workflow acceptance matrix

| Case | Covered by |
|---|---|
| Create/open-like case state | Workbench + CaseModel |
| Navigate setup tree | CaseTree |
| Edit setup | PropertyEditor |
| Validate | SimulationController + SolverAdapter |
| Steady run | RunTarget |
| Transient run | RunTarget + TimeControl |
| Run to physical time | run-until |
| Pause/stop | SimulationController |
| Hot/restart/rebuild edit | ChangeImpact |
| Checkpoint | CheckpointStore |
| Restore | CommandBus + CheckpointStore |
| Geometry/patch selection | SceneModel + Workbench |
| Live monitors | MonitorManager |
| Reports/derived fields | PostProcessor |
| TUI/automation | CommandBus |
| Concrete solver integration | SolverAdapter |
|