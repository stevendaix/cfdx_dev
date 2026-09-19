# CFDX — Guide de la spécification

## Présentation

Ce répertoire contient la spécification technique complète de **CFDX**, un framework CFD généraliste destiné au calcul scientifique haute performance (CPU / GPU / GPU out-of-core / MPI / HPC).

- **Fichier principal :** [`spec.md`](./spec.md)
- **Version :** 0.7
- **Statut :** Architecture cible / spécification de développement

## Structure des documents

| Fichier | Rôle |
|---------|------|
| `spec.md` | Spécification technique complète (sections 1 à 99) |
| `CHECKLIST.md` | Checklist exhaustive de toutes les tâches, organisées par phase |
| `TASKS.md` | Suivi dynamique des tâches — à faire / en cours / fait |

## Comment utiliser ces fichiers

1. **Lire `spec.md`** pour comprendre l'architecture, les principes et les décisions figées.
2. **Consulter `CHECKLIST.md`** pour voir l'état global de chaque tâche (case à cocher).
3. **Mettre à jour `TASKS.md`** au fur et à mesure de l'avancement : déplacer une tâche de *à faire* → *en cours* → *fait*.

## Liens utiles

- **OpenFOAM** : [openfoam.org](https://openfoam.org/) — framework de référence pour les tests de comparaison (§80)
- **Installation locale** : `/opt/openfoam13` (OpenFOAM v13, GNU GPL v3)
  - Pour l'activer dans un shell : `source /opt/openfoam13/etc/bashrc`
  - Binaires : `/opt/openfoam13/bin/` (ex: `foam`, `simpleFoam`, `icoFoam`, ...)
  - Tutorials : `/opt/openfoam13/tutorials/`
  - Code source : `/opt/openfoam13/src/`
- **meshio** : [github.com/nschloe/meshio](https://github.com/nschloe/meshio) — adaptateur d'import (§23)
- **HDF5** : [www.hdfgroup.org](https://www.hdfgroup.org/) — format natif de stockage (§6)
- **Gmsh** : [geonumlab.gmsh.info](http://geonumlab.gmsh.info/) — générateur de maillage (§22)

## Sommaire rapide de la spécification

- **§1-4** : Vision, principes architecturaux, objectifs
- **§5-6** : Format de cas CFDX, structure HDF5
- **§7-19** : `/meta`, `/case`, `/mesh`, points, faces, owner/neighbour, connectivité cellules, patches, géométrie, faces non planaires, calcul géométrique, invariants, validateur
- **§20-23** : Import/export (OpenFOAM, Gmsh, meshio)
- **§24-32** : Fields, BoundaryField, opérateurs FVM (gradient, interpolation, divergence, laplacien)
- **§33-36** : Algèbre linéaire (SparseMatrix, solveurs, préconditionneurs)
- **§37** : Module 0 — contenu définitif (M0.1 à M0.12)
- **§38-53** : Modèle mémoire CPU, abstraction device, storage states, execution policy, runtime, memory planner, GPU full, GPU OOC, halo management, async transfers, pinned memory, performance policy, memory-aware execution
- **§54-58** : MPI, domain decomposition, checkpoint/restart
- §59-61 : Output, results, `data.h5`
- **§62-64** : Python API, C++ API, architecture cible
- **§65-71** : Modules physiques M1 à M7
- **§72-75** : Configuration, reproductibilité, déterminisme, précision
- **§76-81** : Validation analytique, tests, conservation, comparaison OpenFOAM, niveaux de validation
- **§82-86** : Benchmarks, profiling, logging
- **§87-90** : CLI, inspection HDF5, robustesse, compatibilité
- **§91-94** : Séparation case/runtime, execution graph, principes GPU/OOC
- **§95** : Roadmap (Phase 0 à Phase 10)
- **§96** : Critères de sortie du Module 0
- **§97-99** : Architecture finale, décisions figées, principe directeur

## Décisions d'architecture figées (v0.7)

1. HDF5 est le format natif du cas.
2. `case.cfdx.h5` est autoportant.
3. La topologie est la source de vérité du maillage.
4. La géométrie est dérivée de la topologie.
5. OpenFOAM est une référence, pas la définition de CFDX.
6. meshio est un adaptateur, pas une dépendance du cœur.
7. Les Fields sont indépendants du backend matériel.
8. L'algèbre linéaire est indépendante de la physique.
9. Les opérateurs FVM sont indépendants du backend.
10. SIMPLE/PISO/PIMPLE/Rhie-Chow appartiennent au Module 1.
11. CPU/GPU/GPU-OOC sont des Execution Policies.
12. Aucun fallback GPU→CPU silencieux pendant un calcul.
13. Le choix CPU/GPU/GPU-OOC est fait par le Runtime avant exécution.
14. Le GPU-OOC utilise domain decomposition + tiles + halos.
15. Les transferts CPU↔GPU ne doivent pas être présents dans la boucle GPU normale.
16. Pinned memory est limitée à des buffers de staging.
17. La physique ne connaît jamais CUDA/MPI/mémoire GPU directement.
18. Python orchestre ; C++ calcule.
19. Le parallélisme est une propriété du Runtime, pas de la physique.
20. Les benchmarks sont mesurés et non définis arbitrairement à l'avance.

## Principes fondamentaux

> **La physique ne doit pas connaître le backend d'exécution.**

Une équation physique doit pouvoir être exécutée sur CPU, GPU, mode GPU out-of-core, en MPI, sans modifier sa formulation.

> **La topologie reste la source de vérité.**

Si la géométrie stockée est absente ou invalide, CFDX doit pouvoir la reconstruire.

> **Module 0 avant la physique.**

Aucun solveur physique complet ne doit être développé avant validation du noyau de calcul.