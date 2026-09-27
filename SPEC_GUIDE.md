# CFDX — Guide de la spécification

## Présentation

Ce répertoire contient la spécification technique complète de **CFDX**, un framework CFD généraliste destiné au calcul scientifique haute performance (CPU / GPU / GPU out-of-core / MPI / HPC).

- **Fichier principal :** [`spec.md`](./spec.md)
- **Version :** 0.8
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

## Contrat des artefacts de calcul figé (v0.8)

CFDX utilise trois classes d'artefacts distinctes :

| Artefact | Rôle | Source de vérité |
|---|---|---|
| `<case>.cfdx.h5` | définition complète du cas : mesh, setup, physique, BC, IC, numerics, solver | **oui** |
| `<case>.dat.h5` | état numérique : champs calculés, itération, temps, IDs globaux, restart | non, dépend du case |
| `<case>_<time>.vtu` | visualisation / post-traitement | non |

Règles impératives :

- `case.cfdx.h5` est autonome et valide sans DAT.
- `case.cfdx.h5` n'est jamais un restart numérique.
- `case.cfdx.h5` ne stocke pas l'itération courante ni les champs de continuation.
- `<case>.dat.h5` ne contient pas le setup source et ne remplace pas le case.
- Un restart consomme explicitement un case compatible et un DAT compatible.
- VTU est une sortie de visualisation et ne doit jamais être utilisé pour reconstruire le setup.

## Décisions d'architecture figées (v0.8)

1. HDF5 est le format natif.
2. `case.cfdx.h5` est autoportant et constitue la source de vérité de la définition complète du cas.
3. `<case>.dat.h5` est l'état numérique/checkpoint/restart séparé.
4. `<case>_<time>.vtu` est un artefact de visualisation/post-traitement.
5. `case.cfdx.h5` n'est pas un restart numérique.
6. La topologie est la source de vérité du maillage.
7. La géométrie est dérivée de la topologie.
8. OpenFOAM est une référence, pas la définition de CFDX.
9. meshio est un adaptateur, pas une dépendance du cœur.
10. Les Fields sont indépendants du backend matériel.
11. L'algèbre linéaire est indépendante de la physique.
12. Les opérateurs FVM sont indépendants du backend.
13. SIMPLE/PISO/PIMPLE/Rhie-Chow appartiennent au Module 1.
14. CPU/GPU/GPU-OOC sont des Execution Policies.
15. Aucun fallback GPU→CPU silencieux pendant un calcul.
16. Le choix CPU/GPU/GPU-OOC est fait par le Runtime avant exécution.
17. Le GPU-OOC utilise domain decomposition + tiles + halos.
18. Les transferts CPU↔GPU ne doivent pas être présents dans la boucle GPU normale.
19. Pinned memory est limitée à des buffers de staging.
20. La physique ne connaît jamais CUDA/MPI/mémoire GPU directement.
21. Python orchestre ; C++ calcule.
22. Le parallélisme est une propriété du Runtime, pas de la physique.
23. Les benchmarks sont mesurés et non définis arbitrairement à l'avance.

## Principes fondamentaux

> **La physique ne doit pas connaître le backend d'exécution.**

Une équation physique doit pouvoir être exécutée sur CPU, GPU, mode GPU out-of-core, en MPI, sans modifier sa formulation.

> **La topologie reste la source de vérité.**

Si la géométrie stockée est absente ou invalide, CFDX doit pouvoir la reconstruire.

> **Module 0 avant la physique.**

Aucun solveur physique complet ne doit être développé avant validation du noyau de calcul.