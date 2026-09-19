# CFDX — Décisions d'architecture

## D-001 — Langage et outils

- **Décision :** C++17 pour le cœur numérique, Python 3.10 pour l'orchestration (pybind11 pour l'API Python).
- **Justification :** Conformité à la spécification §63-64 (Python orchestre, C++ calcule). Disponible sur le système.
- **Date :** 2026-09-19

## D-002 — Système de build

- **Décision :** CMake avec Ninja si disponible, Make en fallback.
- **Justification :** Standard HPC, bon support MPI/CUDA/pybind11. CMake 4.3.2 et Ninja disponibles.
- **Date :** 2026-09-19

## D-003 — HDF5

- **Décision :** HDF5 1.10.7 construit localement en statique dans `third_party/hdf5` (pas de dev headers système).
- **Justification :** Aucun paquet `libhdf5-dev` installable sans root. HDF5 1.10.7 est la version compatible avec l'environnement (libhdf5-103 = 1.10.7).
- **Date :** 2026-09-19
- **Impact :** Aucune dépendance externe requise pour le build CFDX.

## D-004 — Tests

- **Décision :** Tests en C++ natifs (pas de framework tiers pour l'instant), exécutés via CMake `add_test` + `ctest`. Tests numériques en Python (numpy/scipy) pour les références analytiques.
- **Justification :** Aucun framework de test C++ installé (Catch2 non disponible). pytest disponible pour la couche Python.
- **Date :** 2026-09-19

## D-005 — Mémoire

- **Décision :** `std::vector` pour le stockage CPU de base (M0). Layout SoA (Structure of Arrays) privilégié pour les opérateurs FVM.
- **Justification :** Spécification §38 — cache, SIMD, vectorisation, accès séquentiel. Le SoA est naturellement compatible avec les opérateurs FVM.
- **Date :** 2026-09-19

## D-006 — GPU

- **Décision :** Le GPU est désactivé pour le Module 0 (pas de toolkit CUDA complet — nvcc 11.5 présent mais pas d'headers). Le `ExecutionPolicy::GPU` et `GPU_OUT_OF_CORE` sont réservés pour les phases ultérieures.
- **Justification :** Spécification §39 — "GPU development occurs as an execution capability and must not contaminate the physics architecture." M0 doit être CPU-only.
- **Date :** 2026-09-19

## D-007 — API publique

- **Décision :** Les headers CFDX exposent des interfaces conceptuelles (Mesh, Field, BoundaryField, FVM operators). Les détails HDF5, MPI, CUDA sont cachés derrière des PIMPL ou namespaces internes.
- **Justification :** Spécification §46.
- **Date :** 2026-09-19

## D-008 — OpenFOAM

- **Décision :** OpenFOAM 13 à `/opt/openfoam13` est utilisé comme référence de comparaison (§80), pas comme dépendance de build.
- **Justification :** Spécification §80 — "OpenFOAM est une référence de comparaison, pas la définition mathématique de la vérité."
- **Date :** 2026-09-19