# CFDX GUI workflow validation

The GUI now separates setup, execution and results state while preserving one Case/Session source of truth.

| Stage | Contract | GUI | Real solver validation |
|---|---|---|---|
| Mesh import | HDF5 mesh catalog | Yes | Real mesh metadata |
| Setup | Case + typed contracts | Yes | Headless validated |
| Check | ValidationReport | Yes | Machine-readable |
| Initialize | uniform/field | Yes | Backend contract only |
| Run | ExecutionController | Yes | Existing subprocess integration |
| Monitor | SolverMetrics + MonitorSeries | Yes | Iteration/time keyed |
| Pause | process controller | Yes | Tested controller path |
| Save DAT | case_io | Yes | DAT copied/reloaded |
| Restart | restart option | Yes | Requires solver implementing configured option |
| Time series | ResultSeries | Yes | VTK-family discovery |
| 3D | PyVista adapter | Optional | contour/slice/glyph available |

A GUI integration test must not label a restart as numerically validated unless a solver executable actually consumes the DAT. The current application layer can carry the DAT and construct the restart command, but a solver-specific restart acceptance test remains required for a fully numerical end-to-end claim.
