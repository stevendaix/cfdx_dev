# CFDX Workbench moderne : audit, bibliographie, code et feuille de route

**État : lot incrémentiel implémenté dans le dépôt fourni, le 8 novembre 2026.**

Le document de conception transmis par l'utilisateur fixe la direction : une interface unifiée pour explorer, préparer, exécuter, surveiller, analyser et qualifier des simulations CFDX, avec une seule source de vérité scientifique et des budgets de mémoire explicites. Le présent lot apporte des changements concrets dans le code sans prétendre achever l'ensemble du Workbench industriel.

## 1. État du dépôt et changements réellement livrés

Dans l'archive de départ, on retrouve déjà un Workbench Qt avec docks, espaces de travail, palette de commandes, explorateur HDF5, rendu PyVista/VTK, diagnostics de maillage, pilotage du solveur, pipeline de post-traitement et tests GUI. Il faut donc étendre ces composants, non les remplacer.

| Composant | Situation dans la version fournie | Modification de ce lot |
| --- | --- | --- |
| Navigation et docks Qt | Existants | Nouvel espace `Validate` et dock `Qualification Center` |
| Registre de qualification | `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` | Lecteur unique, borné et en lecture seule ; filtres et recherche transversale |
| CLI de suivi | Scripts de campagne déjà présents | `scripts/cfdx_qualification_status.py` avec filtres et JSON exploitable en CI |
| Viewport | Qt + PyVista, contour/slice/glyph | Ajout `Clip` et `Threshold`, édition des paramètres, garde-fous numériques |
| Pipeline de filtres | Opérations ordonnées non destructives | État publié après évaluation/rendu ; rollback d'affichage sur erreur de rendu |
| CI graphique | `.github/workflows/cfdx-gui.yml` | Ajout des nouveaux tests au profil GUI |
| Validation physique NACA/SST/etc. | Registre déclaré, preuves hétérogènes | **Pas de nouvelle convergence physique démontrée** : aucune promotion automatique de statut |

### 1.1 Contrat de provenance pour le suivi transversal

Le registre compte **25 cas** dans l'archive : `PASS=3`, `DIAGNOSTIC=6`, `READY=10`, `BLOCKED=4`, `PLANNED=2`. Ces statuts ne certifient pas le commit courant. Les valeurs sont lues depuis le fichier source, sans duplication d'une base parallèle. Le lecteur rejette les doublons de clés JSON, les identifiants dupliqués, les statuts ou schémas non pris en charge, les constantes non finies et les fichiers de plus de 4 MiB. Une lecture modifiée en cours d'opération est rejetée. L'interface conserve l'ancien instantané valide en cas d'erreur de rechargement.

La vue affiche la recherche texte, les filtres de domaine et de statut, une synthèse et les critères enregistrés. Elle avertit explicitement que **PASS déclaré != PASS confirmé sur HEAD**. Cette distinction est indispensable dans CFDX, notamment pour comparer NACA, Ghia, SST, CHT, convection naturelle, rayonnement, multiphase et tests de solveurs sans focalisation sur une seule famille.

```sh
python scripts/cfdx_qualification_status.py
python scripts/cfdx_qualification_status.py --status DIAGNOSTIC
python scripts/cfdx_qualification_status.py --format json --query SST
```

Le JSON expose `head_evidence: "unverified"`, sans inventer une exécution scientifique.

### 1.2 Pipeline moderne non destructif

Les opérations `contour`, `slice`, `glyph`, `clip`, `threshold` sont des données immuables. Chaque évaluation repart de la source originale. La validation rejette NaN, infini, normale nulle, bornes inversées, champ absent et opération inconnue. `clip(..., inplace=False)` conserve les données sources ; `threshold` utilise une préférence explicite pour les données cellulaires, cohérente avec la documentation PyVista. Avant publication d'un changement, le résultat est calculé ; en cas d'échec du rendu, le système tente de restaurer la dernière scène valide et ne publie pas le nouveau pipeline.

Ces changements ne constituent **pas** un moteur de streaming de maillages industriels. Les filtres fonctionnent encore sur des datasets chargés en mémoire. Une projection de performance de 10 ou 50 millions de cellules ne peut donc être déclarée validée sur cette base.

## 2. Architecture cible et principes de programmation

```text
CFDX Workbench (Qt: navigation, selection, proprietes, statut)
  |-- Project Explorer + HDF5 Inspector (metadata first)
  |-- Simulation Setup (conditions limites et modeles contextuels)
  |-- Scientific Viewport (pipelines VTK/PyVista, budgets memoire)
  |-- Simulation Cockpit (processus, logs, residus, ressources)
  |-- Results & Analysis (profils, flux, forces, comparaisons)
  `-- Validation Center (cas, preuves, rapports, seuils)
                  |
      Services / modele scientifique partage
         |-- lecteur HDF5 / CGNS / formats de resultat
         |-- operateurs de visualisation immuables
         |-- superviseur des calculs / checkpoints
         `-- registre de preuves versionnees
                  |
               CFDX core (C++ et solveurs, source de verite)
```

Regles : pas de calcul physique dans les widgets, pas de valeur de démonstration affichée comme un résultat réel, aucune solution déclarée convergée sur le seul code de retour du processus, et aucune modification silencieuse des seuils. Les cas sources et checkpoints sont de nature différente ; l'interface doit afficher cette différence.

## 3. Recherche bibliographique et décisions techniques

| Bibliographie officielle | Idée pertinente pour CFDX | Choix |
| --- | --- | --- |
| [ParaView - Remote and parallel visualization](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html) | Séparation logique client / données / rendu ; traitement distant | **P1/P2 : adopter une architecture similaire**, sans dupliquer immédiatement un serveur ParaView |
| [VTK - vtkStreamingDemandDrivenPipeline](https://vtk.org/doc/nightly/html/classvtkStreamingDemandDrivenPipeline.html) | Demander des pièces ou sous-étendues du maillage | **P1 : lecteur par blocs**, à qualifier sur données réelles |
| [PyVista - Using Common Filters](https://docs.pyvista.org/examples/01-filter/using_filters.html) | `slice`, `clip`, `threshold`, `extract_surface`, etc. | **Implémenté pour Clip/Threshold** dans le pipeline actuel |
| [h5py - Datasets, chunking](https://docs.h5py.org/en/stable/high/dataset.html) | Lectures par sélection, chunking, stockage compressé | **P1 : conserver un aperçu borné**, ajouter l'extraction bornée et un index spatial |
| [Qt for Python - Threading Technologies](https://doc.qt.io/qtforpython-6/overviews/qtdoc-threads-technologies.html) | Workers, pools, signals ; ne jamais bloquer le thread UI | **P1 : étendre aux lectures importantes**, avec annulation et invalidation |
| [OpenFOAM - Boundary conditions](https://doc.openfoam.com/2212/tools/processing/boundary-conditions/) | Groupes inlet/outlet/wall/coupled/constraint | **P1 : fiche contextuelle par type**, validation de la combinaison physique et des dimensions |
| [SU2_GUI (source)](https://github.com/su2code/su2gui) | Interface de configuration, lancement et post-traitement | **P2 : s'inspirer des parcours utilisateur**, sans réimplémenter SU2GUI |
| [PETSc - SNES nonlinear solvers](https://petsc.org/release/manual/snes/) | Newton, Newton-Krylov, line search, trust region, Anderson, NGMRES | **P1/P2 : cockpit des opérateurs**, montrer le vrai résidu et les méthodes disponibles ; aucune méthode numérique ajoutée ici |
| [NASA - Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/) | Cas de vérification/validation turbulence et données de référence | **P1/P2 : lier les cas SST aux références**, avec convergence spatiale et erreurs physiques indépendantes |

### 3.1 Bonnes idées à intégrer, sans copier les logiciels concurrents

**Ergonomie façon Fluent.** Une navigation `General > Models > Materials > Cell Zones > Boundary Conditions > Methods > Initialization > Run`, tandis que chaque sélection d'une face présente uniquement les entrées utiles. Le schéma de condition limite doit déclarer champs requis, unité, type, compatibilité de modèle, validation de jauge de pression, détection d'une sortie sans référence et traitement des retours de flux. Les validations restent dans le cœur scientifique ou un service métier, pas dans Qt.

**Diagnostic plus fort que le simple résidu.** Le cockpit doit juxtaposer résidus linéaires vrais, résidus non linéaires sur l'état retourné, bilan de masse local/global, flux Rhie-Chow, erreurs énergétiques, forces et moniteurs spécifiques. Les normalisations et normes doivent être exposées. Sur SIMPLE/SIMPLEC/COUPLED/Newton-Krylov, un affichage commun n'implique pas que les algorithmes soient équivalents.

**Rapports de qualification cross-cases.** Chaque preuve doit contenir au minimum `git_sha`, empreinte du maillage et du cas, version du solveur, date et commande de lancement, résultat du processus, critères numériques, résidus vrais, bilans, comparaisons indépendantes, artefacts et politique de validité. **Le registre de planification seul ne suffit pas.** L'affichage doit indiquer `NOT_RUN`, `RUNNING`, `BLOCKED`, `FAILED`, `PASSED`, `STALE` au niveau d'une preuve sur HEAD, indépendamment du statut de travail déclaratif.

**Streaming et rendu industriel.** L'aperçu doit être adapté au budget RAM/VRAM et au nombre de cellules ; lecteur des metadonnees, extraction de surface à partir des faces visibles, partitions et cache avec clé d'invalidation. Pour 1M/10M/50M cellules, mesurer ouverture, premier affichage, p95 de navigation, consommation de mémoire, temps de filtre, annulation, déchargement ; un benchmark de HDF5 seul ne valide pas le rendu 3D.

**Reproductibilité.** Les filtres et les vues doivent pouvoir être sérialisés par une spécification versionnée, rejouée et comparée entre révisions du solveur ; l'export doit conserver provenance, unités et paramètres.

## 4. Plan priorisé, critères de fin

| Priorité | Travail à effectuer | Critère d'acceptation |
| --- | --- | --- |
| **P0 (ce lot)** | Vue Validate, recherches transversales, Clip/Threshold, tests | Registre unique/borné, pipeline non destructif, échec explicite, tests ciblés verts |
| **P0 suivant** | Vérification native GUI et OpenGL en CI réelle | Ouverture Qt, actions utilisateur et captures, aucune regression ; 0 skip obligatoire |
| **P1** | Ingestion de preuves de qualification versionnées | Aucune promotion sans SHA/artefacts ; tests de preuves périmées et campagnes partielles |
| **P1** | Wizard BC type Fluent, validateurs physiques, preflight | Messages actionnables, impossibilité d'exécuter un cas physiquement incomplet |
| **P1** | Lecture blocs HDF5, extraction de surface, workers interruptibles | Budget RAM testé avec 1M et 10M, erreur explicite si 50M non pris en charge |
| **P1** | Cockpit de convergence multi-méthodes | Separation résidu linéaire/non linéaire/conservation, suivi des cas en parallèle |
| **P2** | Visualisation distante, comparaison de campagnes, CGNS, sessions rejouables | Benchmarks et invariants physiques, export reproductible |
| **P2** | Bibliothèque de cas physiques NASA/NACA/Ghia/CHT/radiation/VOF | Preuves sur le HEAD exact et comparaisons indépendantes, sans assouplir les seuils |

## 5. Vérification effectivement exécutée dans cet environnement

- Groupe registre/pipeline/CLI/Workbench/inspecteur HDF5 : **81 tests réussis, 20 ignorés**, environ 6,19 secondes.
- Groupe qualité de maillage, topologie, navigateur et tests du moteur 3D : **48 tests réussis, 10 ignorés**, environ 1,62 seconde.
- Total de ces deux groupes disjoints : **129 réussis et 30 ignorés**. Aucun échec dans ces groupes.
- Les dépendances PySide6, PyVista et PyVistaQt ne sont pas disponibles localement ; tests Qt/OpenGL réellement executes : **aucun**. Les tests ignorés ne sont pas des reussites.
- La campagne Python complète n'est **pas certifiée** : `meshio` manque pour la collecte d'un module et le parcours élargi sans ce module a dépassé le temps local d'exécution. Aucune CI GitHub distante ni tests C++ complets exécutés pour cette livraison.
- Aucun résultat NACA, SST, thermique ou radiatif n'a été reconvergé dans ce lot : **aucune qualification physique nouvelle**.

Reproduire depuis la racine du dépôt :

```sh
PYTHONPATH=src/cfdx/python python -m pytest -q \
  tests/python/test_qualification_center.py \
  tests/python/test_visualization_pipeline.py \
  tests/python/test_visualization_renderer_transaction.py \
  tests/python/test_capabilities.py \
  tests/python/test_workbench.py \
  tests/python/test_case_inspector.py \
  tests/python/test_hdf5_inspector_benchmark.py

PYTHONPATH=src/cfdx/python python -m pytest -q \
  tests/python/test_mesh_quality.py \
  tests/python/test_volume_mesh_quality.py \
  tests/python/test_mesh_topology_validation.py \
  tests/python/test_gui_browser_widgets.py \
  tests/python/test_renderer.py \
  tests/python/test_renderer_headless_contract.py \
  tests/python/test_renderer_pyvista.py \
  tests/python/test_gui_3d.py \
  tests/python/test_specialized_setup.py
```

Pour une validation GUI complète : `python -m pip install -e '.[gui]'` et utiliser la CI `cfdx-gui.yml` sur le HEAD exact de la branche, avec les dependances système OpenGL disponibles. Il est indispensable de contrôler les tests ignorés et de conserver les screenshots/JUnit en artefacts de CI.

## 6. Mode d'emploi et integration

1. Partir de la base exacte correspondant à l'archive fournie : ne pas rebaser le patch en aveugle sur une PR modifiée.
2. Ouvrir le Workbench et choisir l'espace `Validate` pour consulter le registre déclaratif et ses critères.
3. Dans `Analyze`, utiliser `Clip` ou `Threshold` sur un champ compatible puis consulter la pile des opérations. L'opération ne touche pas au champ source.
4. Utiliser le CLI de statut en CI ou pendant une revue technique ; ne jamais confondre le JSON produit avec un certificat numérique.
5. Exécuter les tests GUI en environnement complet avant la promotion du code.

**Livraison :** ZIP du dépôt mis à jour, patch Git cumulatif, ce rapport et un Markdown autonome contenant le patch complet. Ce lot est une évolution logicielle partiellement vérifiée, non une certification industrielle et numérique complète.
