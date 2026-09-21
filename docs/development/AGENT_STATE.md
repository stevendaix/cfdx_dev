# CFDX — État de l'agent

## Dernière mise à jour : 2026-09-21

## Phase actuelle

**M0 — Stabilisation finale / validation des imports** → M1 — Navier-Stokes incompressible

La branche `audit/mesh-importers` complète l'audit de l'import maillage. M1 ne doit démarrer qu'après validation CI de cette PR.

## Tâches complétées

- [x] Socle Mesh, géométrie, qualité et champs
- [x] HDF5 writer/reader + round-trip + hash
- [x] Algèbre linéaire : CG, BiCGStab, GMRES, préconditionneurs et interfaces
- [x] MPI : décomposition, ghost cells et halo exchange
- [x] Temporal : Euler, Crank-Nicolson, BDF2, LTS/adaptation
- [x] EOS et modèles de transport
- [x] Diagnostics, résidus et profiling
- [x] Memory Planner / Mesh Reordering
- [x] Audit et durcissement des imports OpenFOAM / meshio / Gmsh
- [x] Tests géométriques synthétiques pour les imports
- [x] CI Release + DebugSanitizers + tests Python meshio

## Validation encore requise

- [ ] CI de la branche `audit/mesh-importers` verte en Release
- [ ] CI de la branche `audit/mesh-importers` verte en DebugSanitizers
- [ ] CTest complet vert
- [ ] Tests meshio/Gmsh/VTU de la matrice géométrique verts
- [ ] Merge de la PR #3 après validation CI

## Points techniques connus

- CUDA reste désactivé tant qu'un toolkit complet n'est pas disponible.
- HDF5 est fourni via `third_party/hdf5` dans l'environnement actuel.
- OpenFOAM 13 est une référence de comparaison et non une dépendance de build.
- La couverture des formats meshio dépend des readers disponibles dans la version installée ; la topologie doit également être représentable par le modèle volumes/faces/owner-neighbour de CFDX.

## Prochaine étape

Après merge et CI verte de la PR #3 :

1. figer M0 dans la documentation ;
2. démarrer M1 avec diffusion/convection ;
3. construire les équations de quantité de mouvement et de continuité ;
4. implémenter le couplage pression-vitesse ;
5. ajouter les cas de validation cavity et Poiseuille.

## Décisions architecturales

Voir `docs/development/DECISIONS.md`.
