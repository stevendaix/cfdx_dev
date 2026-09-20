# CFDX — État de l'agent

## Dernière mise à jour : 2026-09-20 (Memory Planner v4 applied)

## Phase actuelle

Phase 0 — Repository Bootstrap (M0 completed)

## Tâche actuelle

M0.10-T02 — Gmsh importer (XML .msh via meshio) + Memory Planner v4 intégration (budget/reordering)

## Tâches complétées

- [x] Repository GitHub créé et pushé
- [x] Spécification v0.7 commitée (dc376cb)
- [x] SPEC_GUIDE.md, CHECKLIST.md, TASKS.md créés
- [x] prompt.md (agent master) commité
- [x] docs/development/ ROADMAP.md, BACKLOG.md, DECISIONS.md créés
- [x] Arborescence source créée (src/cfdx/...)
- [x] pybind11 installé (3.1.0)
- [x] HDF5 1.10.7 en cours de build (background)
- [x] M0.9-T02: Fix HDF5 reader memory bug (H5Tget_size instead of H5Aget_storage_size)
- [x] M0.9-T03: Round-trip validation (HDF5 case round-trip, all 7 tests pass)
- [x] M0.9-T04: Hash computation (SHA256 checksums for case.cfdx.h5 self-contained integrity)
- [x] M0.10-A1: Memory Planner v4 (memory_planner.md → src/cfdx/core/memory/)
- [x] Dernier commit validé: 6e69c23 — M0.10-A1 Memory Planner v4 skeleton

## Tâches bloquées

Aucune

## Problèmes connus

- HDF5 dev headers absents du système → HDF5 construit localement en statique (third_party/hdf5)
- CUDA toolkit incomplet (nvcc 11.5, pas d'headers) → GPU désactivé pour M0
- Aucun sudo → pas d'apt install
- Catch2 non disponible → tests C++ natifs (pas de framework tiers)

## Dernier commit validé

ec7763d — "M0.9-T03/T04: Update agent state - round-trip validation complete, hash computation next; AGENT_STATE and BACKLOG synced"

## Prochaine tâche

M0.10-T02 — Gmsh importer + Memory Planner v4 budget/reordering tests

## Décisions architecturales

Voir docs/development/DECISIONS.md (D-001 à D-008)

## Décisions humaines en attente

Aucune

## Dépendances externes

| Dépendance | Statut |
|-----------|--------|
| C++17 (gcc 11.4) | OK |
| CMake 4.3.2 + Ninja | OK |
| Python 3.10 + numpy/scipy/h5py/meshio/gmsh/pytest | OK |
| pybind11 3.1.0 | OK (installé) |
| HDF5 1.10.7 | En cours de build (third_party/hdf5) |
| MPI (mpicc/mpirun) | OK |
| CUDA toolkit | Incomplet — GPU reporté |
| OpenFOAM 13 (/opt/openfoam13) | OK (référence comparaison) |