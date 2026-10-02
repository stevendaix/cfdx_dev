# CFDX Developer Guide

**Status: PARTIAL — architecture exists; detailed migration is ongoing.**

The Developer Guide explains how CFDX is organised, extended, tested and integrated. It complements the Theory Guide (mathematics) and the V&V Guide (evidence).

## Development contract

- C++20 is the numerical core.
- CMake is the build-system boundary.
- Python/pybind11 is the Python-facing boundary where applicable.
- Numerical contracts are explicit: no silent algorithm fallback.
- A failing verification or qualification test is evidence to investigate, not a reason to relax a tolerance.
- Source paths and API names must be checked against the repository before publication.

## Developer domains

| # | Chapter | Purpose |
|---|---|---|
| 00 | Development philosophy | Engineering rules and scientific integrity |
| 01 | Repository architecture | Source/test/documentation boundaries |
| 02 | Build system | CMake, configurations and dependencies |
| 03 | Code architecture | Core/physics/API layering |
| 04 | C++ guidelines | Interfaces, ownership and numerical code |
| 05 | Python bindings | pybind11 boundary and API design |
| 06 | Data model and I/O | Case/state/output semantics |
| 07 | Numerical method registry | Schemes and method selection |
| 08 | Physics development | Adding equations and closures |
| 09 | Mesh and geometry development | Geometry/topology contracts |
| 10 | Linear algebra development | Operators, solvers and preconditioners |
| 11 | Solver development | Nonlinear and coupled algorithms |
| 12 | Testing | Unit, numerical and integration tests |
| 13 | Debugging | Diagnostics and failure isolation |
| 14 | Performance | Profiling and reproducible optimisation |
| 15 | Parallelism | Domain decomposition and communication |
| 16 | GPU | Explicit accelerator paths and fallback policy |
| 17 | CI/CD | Automated build and verification gates |
| 18 | PR workflow | Review, evidence and merge discipline |
| 19 | Release and maintenance | Reproducibility and compatibility |

## Cross-domain rule

For every numerical feature:

Theory equation → Developer implementation contract → V&V test/evidence → qualification decision.

The three guides must use the same terminology and status vocabulary.

Existing material under docs/development/ and the existing developer-domain trees remains authoritative during migration until reconciled.
