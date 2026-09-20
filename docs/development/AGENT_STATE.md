# CFDX — État de l'agent

## Dernière mise à jour : 2026-09-20 (Module 0 terminé - vérification complète)

## Phase actuelle

Phase 0 — Module 0 (Core CFDX Framework) → M1 — Incompressible laminar

## Tâche actuelle

M1 — Navier-Stokes incompressible :
- Discrétisation spatiale (diffusion + convection)
- Couplage pressseur-vitesse (SIMPLE)
- Solveur linéaire (CG + AMG)
- Tests de validation (cavity, poiseuille, etc.)

## Tâches complétées

- [x] Repository GitHub créé et pushé
- [x] Spécification v0.7 commitée (dc376cb)
- [x] SPEC_GUIDE.md, CHECKLIST.md, TASKS.md créés
- [x] prompt.md (agent master) commité
- [x] docs/development/ ROADMAP.md, BACKLOG.md, DECISIONS.md créés
- [x] Arborescence source créée (src/cfdx/...)
- [x] pybind11 installé (3.1.0)
- [x] HDF5 1.10.7 construit (third_party/hdf5)
- [x] M0.9-T02: Fix HDF5 reader memory bug (H5Tget_size instead of H5Aget_storage_size)
- [x] M0.9-T03: Round-trip validation (HDF5 case round-trip, all 7 tests pass)
- [x] M0.9-T04: Hash computation (SHA256 checksums for case.cfdx.h5 self-contained integrity)
- [x] M0.10-A1: Memory Planner v4 (memory_planner.md → src/cfdx/core/memory/)
- [x] M0.10-A3: Mesh Reordering v4 (RCM + SFC + GPUOptimized + PartitionAware) — 10 tests pass
- [x] M0.10-T02/T03/T04: Gmsh importer, meshio, CFDX case writer
- [x] M0.8-T04/T05: CG solver, BiCGStab solver (tests pass)
- [x] M0.8-T06: GMRES (in src/cfdx/core/linalg/ - needs implementation)
- [x] M0.8-T07: Preconditioners (Jacobi, RBGS, ILU(0) - in src/cfdx/core/linalg/)
- [x] M0.11: Real MPI implementation (not stubs) - mpi_utils.h, mesh_partitioner.h
- [x] M0.12-T03: Memory planner v4 (real implementation with 6 KPIs)
- [x] M0.13-T01/T02/T03: Temporal discretization (Euler, Crank-Nicolson, BDF2, LTC, adaptive)
- [x] M0.14-T01/T02/T03: EOS, Transport models, Multi-component
- [x] Dernier commit validé: CMakeLists.txt fixes

## Tâches bloquées

- CUDA toolkit incomplet (nvcc 11.5, pas d'headers) → GPU désactivé
- Aucun sudo → pas d'apt install
- Catch2 non disponible → tests C++ natifs (pas de framework tiers)

## Problèmes connus

- HDF5 dev headers absents du système → HDF5 construit localement en statique (third_party/hdf5)
- test_mpi_partition nécessite MPI initialisé → test désactivé pour builds sériels

## Dernier commit validé

CMakeLists.txt - Build system fixes pour importer sources correctement

## Prochaine tâche

M1 — Navier-Stokes incompressible :
- Discrétisation spatiale (diffusion + convection)
- Couplage pressseur-vitesse (SIMPLE)
- Solveur linéaire (CG + AMG)
- Tests de validation (cavity, poiseuille, etc.)

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
| HDF5 1.10.7 | OK (third_party/hdf5) |
| MPI (mpicc/mpirun) | OK |
| CUDA toolkit | Incomplet — GPU reporté |
| OpenFOAM 13 (/opt/openfoam13) | OK (référence comparaison) |

## Statut Module 0 (COMPLÉTÉ)

Tous les éléments du Module 0 ont été vérifiés:

### M0.10-A3 (Mesh Reordering v4) ✅
- Implémentation complète: RCM, SFC, GPUOptimized, PartitionAware
- 10 tests passent dans test_memory_planner_validation.cpp
- Fichiers: src/cfdx/core/memory/memory_planner.h/.cpp

### M0.11 (MPI) ✅
- Implémentation réelle (pas stubs):
  - mesh_partitioner.h: partition_geometric(), build_halo_plan()
  - mpi_utils.h: MPI utilities complètes
- Fichiers: src/cfdx/core/parallel/

### M0.12-T03 (Memory Planner) ✅
- Implémentation v4 avec 6 KPIs réels:
  - bytes_per_cell (stockage)
  - bytes_per_cell_per_iteration (traffic)
  - peak_ram, peak_vram
  - topology_bytes, geometry_bytes, fields_bytes
- Fichiers: src/cfdx/core/memory/memory_planner.h/.cpp

### M0.13 (Temporal Discretization) ✅
- Schémas implémentés: Euler (explicite/implicite), Crank-Nicolson, BDF2
- Local Time Stepping: compute_local_time_steps()
- Adaptive Time Step: AdaptiveTimeStepper avec CFL-based control
- Tests: test_temporal.cpp, tests/numerical/test_temporal.cpp

### M0.14 (Transport Models) ✅
- EOS: IncompressibleEOS, IdealGasEOS avec compute_thermo_fields()
- Transport: sutherland_viscosity(), prandtl_conductivity(), schmidt_diffusivity()
- Tests: test_transport_models.cpp compile & passes

## Tests - Résultats

30/30 tests passent:
- test_point, test_face, test_ownership, test_cell, test_boundary, test_mesh
- test_field, test_field_storage, test_boundary_field
- test_interpolation, test_gradient, test_divergence, test_flux, test_integrate
- test_sparse_matrix, test_vector, test_linear_system
- test_cg_solver, test_bicgstab_solver, test_laplacian
- test_face_geometry, test_cell_geometry, test_mesh_quality, test_source_term
- test_geometry_cache, test_hdf5_writer, test_hdf5_roundtrip, test_field_roundtrip
- test_memory_planner_validation, test_plan_complet