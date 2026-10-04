# CFDX Developer Guide

**Status: PARTIAL — the 20-chapter structure and core contracts exist; detailed migration of legacy developer material is ongoing.**

The Developer Guide explains how CFDX is organised, extended, tested and integrated. It complements Theory (mathematics) and V&V (evidence).

## Development contract

- C++20 is the numerical core.
- CMake is the build-system boundary.
- Python/pybind11 is the Python-facing boundary where applicable.
- Numerical contracts are explicit: no silent algorithm fallback.
- Failed verification/qualification is evidence to investigate, not a reason to relax tolerances.
- Source paths and API names must be checked against the repository before publication.

```{toctree}
:maxdepth: 2

00-development-philosophy/chapter
01-repository-architecture/chapter
02-build-system/chapter
03-code-architecture/chapter
04-cpp-guidelines/chapter
05-python-bindings/chapter
06-data-model-and-io/chapter
07-numerical-method-registry/chapter
08-physics-development/chapter
09-mesh-and-geometry-development/chapter
10-linear-algebra-development/chapter
11-solver-development/chapter
12-testing/chapter
13-debugging/chapter
14-performance/chapter
15-parallelism/chapter
16-gpu/chapter
17-ci-cd/chapter
18-pr-workflow/chapter
19-release-and-maintenance/chapter
```

## Cross-domain rule

**Theory equation → Developer implementation contract → V&V test/evidence → qualification decision.**

Existing material under `docs/development/` remains a migration source until reconciled. It must not silently become a second authority for claims already covered by this guide.
