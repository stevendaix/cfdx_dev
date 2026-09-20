# CFDX — Backlog

Dernière mise à jour : 2026-09-19

## PLANNED

| ID | Title | Dependencies | Status |
|----|-------|-------------|--------|
| M0.1-T01 | Point storage | — | DONE |
| M0.1-T02 | Face CSR connectivity | M0.1-T01 | DONE |
| M0.1-T03 | Owner/neighbour storage | M0.1-T01, M0.1-T02 | DONE |
| M0.1-T04 | Cell CSR connectivity | M0.1-T01, M0.1-T02 | DONE |
| M0.1-T05 | Boundary patches | M0.1-T01, M0.1-T02 | DONE |
| M0.1-T06 | Topology validation | M0.1-T01..T05 | DONE |
| M0.2-T01 | Face geometry | M0.1-T01..T06 | DONE |
| M0.2-T02 | Cell geometry | M0.2-T01 | DONE |
| M0.2-T03 | Skewness | M0.2-T01 | DONE |
| M0.2-T04 | Non-orthogonality | M0.2-T01 | DONE |
| M0.3-T01 | Mesh validator | M0.2-T01..T04 | DONE |
| M0.4-T01 | Field<T, Location> template | M0.1-T01..T06 | DONE |
| M0.4-T02 | Field metadata | M0.4-T01 | DONE |
| M0.4-T03 | StorageHandle | M0.4-T01 | DONE |
| M0.5-T01 | BoundaryField / PatchField | M0.4-T01 | DONE |
| M0.6-T01 | Cell→Face interpolation | M0.4-T01 | DONE |
| M0.7-T01 | Gauss gradient | M0.2-T01, M0.4-T01 | DONE |
| M0.7-T02 | Divergence | M0.4-T01 | DONE |
| M0.7-T03 | Laplacian | M0.7-T01 | DONE |
| M0.7-T04 | Flux | M0.4-T01 | DONE |
| M0.7-T05 | Surface/Volume integrate | M0.4-T01 | DONE |
| M0.8-T01 | SparseMatrix CSR | — | DONE |
| M0.8-T02 | Vector | M0.8-T01 | DONE |
| M0.8-T03 | LinearSystem | M0.8-T01, M0.8-T02 | DONE |
| M0.8-T04 | CG solver | M0.8-T03 | DONE |
| M0.8-T05 | BiCGStab solver | M0.8-T03 | DONE |
| M0.8-T06 | GMRES solver | M0.8-T03 | PLANNED |
| M0.8-T07 | Preconditioners | M0.8-T03 | PLANNED |
| M0.9-T01 | HDF5 writer | M0.1-T01..T06, M0.4-T01 | DONE |
| M0.9-T02 | HDF5 reader | M0.9-T01 | DONE |
| M0.9-T03 | Round-trip validation | M0.9-T01, M0.9-T02 | DONE |
| M0.9-T04 | Hash computation | M0.1-T01..T06 | PLANNED |
| M0.10-T01 | OpenFOAM importer | M0.1-T01..T06 | PLANNED |
| M0.10-T02 | Gmsh importer | M0.1-T01..T06 | PLANNED |
| M0.10-T03 | meshio importer | M0.1-T01..T06 | PLANNED |
| M0.10-T04 | CFDX case writer | M0.9-T01 | PLANNED |
| M0.11-T01 | Domain decomposition | M0.1-T01..T06 | PLANNED |
| M0.11-T02 | Ghost cells | M0.11-T01 | PLANNED |
| M0.11-T03 | Halo exchange | M0.11-T02 | PLANNED |
| M0.12-T01 | CPU backend | M0.4-T01, M0.7-T01..T05 | PLANNED |
| M0.12-T02 | ExecutionPolicy enum | — | PLANNED |
| M0.12-T03 | Memory planner | M0.4-T01 | PLANNED |

## IN_PROGRESS

- M0.9-T02 (HDF5 reader) — bug de mémoire dans le lecteur (buffer d'attribut trop petit). En cours de débogage.

## BLOCKED

- Aucune

## COMMITTED

- M0.1 (Mesh topology) — commit `790d007`
- M0.2 (Geometry) — commit `8314db0`
- M0.3 (Mesh quality) — commit `f53c9cd`
- M0.4 (Fields) — commit `d4e24f9`
- M0.5 (Boundary fields) — commit `662b678`
- M0.6 (Interpolation) — commit `2240cd6`
- M0.7 (FVM operators) — commit `e68d3b8`
- M0.8 (Linear algebra) — commit `cdf76d4`
- M0.9-T01 (HDF5 writer) — non commité