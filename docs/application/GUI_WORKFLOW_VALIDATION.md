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
| Restart | restart option + DAT | Yes | **Validated for native steady incompressible U/p;** temperature/turbulence and MPI N→M remain separate |
| Time series | ResultSeries | Yes | VTK-family discovery |
| 3D | PyVista adapter | Optional | contour/slice/glyph available |

The native steady incompressible solver now has a solver-specific DAT acceptance test covering U/p restart consumption and direct-continuation equivalence. This does not imply restart support for temperature/turbulence or MPI N→M. The GUI must therefore report those capabilities separately rather than treating every DAT field as solver-consumable.
