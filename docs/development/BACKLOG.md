# CFDX — Backlog

Dernière mise à jour : 2026-09-21

## M0 — Stabilisation et validation

| ID | Title | Status |
|----|-------|--------|
| M0.1 | Mesh topology | DONE |
| M0.2 | Geometry | DONE |
| M0.3 | Mesh quality | DONE |
| M0.4 | Fields | DONE |
| M0.5 | Boundary fields | DONE |
| M0.6 | Interpolation | DONE |
| M0.7 | FVM operators | DONE |
| M0.8 | Linear algebra | DONE* |
| M0.9 | HDF5 / export | DONE |
| M0.10 | Import/export | **IN VALIDATION** |
| M0.11 | MPI | DONE* |
| M0.12 | Execution abstraction | DONE* |
| M0.13 | Temporal discretization | DONE |
| M0.14 | Thermodynamics & transport | DONE |
| M0.15 | Diagnostics & monitoring | DONE |

\* DONE au niveau fonctionnel/documentaire ; la validation finale CI reste la référence avant de considérer le socle totalement figé.

### M0.10 — Import/export

- [x] Native OpenFOAM reader
- [x] Universal mesh dispatcher
- [x] meshio bridge
- [x] Gmsh through meshio
- [x] Higher-order common topologies reduced to corner topology
- [x] Polyhedral cells when exposed by meshio
- [x] Physical-group / boundary patch mapping
- [x] Topology validation
- [x] OpenFOAM regression test
- [x] meshio bridge regression test
- [x] Synthetic geometry matrix
- [ ] Release CI
- [ ] DebugSanitizers CI

## IN_PROGRESS

- M0.10 — final CI validation of PR #3

## BLOCKED

- None

## NEXT

- M1.1 — Diffusion / convection discretization
- M1.2 — Continuity equation
- M1.3 — Pressure-velocity coupling
- M1.4 — SIMPLE / PISO
- M1.5 — Cavity / Poiseuille validation
