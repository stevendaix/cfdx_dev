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
8. Save/load the HDF5 case, optionally with its paired DAT restart artifact, and resume/restart from the validated state.
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

The GUI/VTK renderer, case persistence service and concrete CFD solver adapter intentionally remain separate. This keeps the core usable on HPC nodes without a display server and allows a future Qt/VTK frontend to use the same controller API.

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
