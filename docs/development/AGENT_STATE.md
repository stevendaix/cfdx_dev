# CFDX — État de l'agent

Dernière mise à jour : 2026-09-19T10:15+02:00

## Phase actuelle

Phase 0 — Repository Bootstrap

## Tâche actuelle

Aucune (démarrage)

## Tâches complétées

- [x] Repository GitHub créé et pushé
- [x] Spécification v0.7 commitée (dc376cb)
- [x] SPEC_GUIDE.md, CHECKLIST.md, TASKS.md créés
- [x] prompt.md (agent master) commité
- [x] docs/development/ ROADMAP.md, BACKLOG.md, DECISIONS.md créés
- [x] Arborescence source créée (src/cfdx/...)
- [x] pybind11 installé (3.1.0)
- [x] HDF5 1.10.7 en cours de build (background)

## Tâches bloquées

Aucune

## Problèmes connus

- HDF5 dev headers absents du système → HDF5 construit localement en statique (third_party/hdf5)
- CUDA toolkit incomplet (nvcc 11.5, pas d'headers) → GPU désactivé pour M0
- Aucun sudo → pas d'apt install
- Catch2 non disponible → tests C++ natifs (pas de framework tiers)

## Dernier commit validé

da93ef9 — "Add spec guide, task tracker and full checklist"

## Prochaine tâche

M0.1-T01 — Point storage (src/cfdx/core/mesh/point.h)

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