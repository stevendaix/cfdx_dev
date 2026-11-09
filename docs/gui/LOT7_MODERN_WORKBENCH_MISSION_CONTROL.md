# CFDX Workbench — Lot 7 : Mission Control, parcours moderne et sûreté des états

## Objectif

Compléter les panneaux Qt du Workbench par une vue d'ensemble **branchée sur les événements réels** de l'application. Le thème par défaut reste **OLED Black** (`#000000`). Les vues `Setup`, `Mesh`, `Run`, `Analyze`, `Validate` demeurent ; `Overview` devient un sixième espace (Ctrl+0).

## Fonctionnalités livrées

1. **Mission Control** (`Overview`) affiche le nom du projet, les sections de configuration réellement renseignées, le statut d'un contrôle du cas, l'état du solveur, les itérations/temps et les résultats chargés. La liste d'anomalies donne accès au workspace concerné par double-clic ou par le bouton « Open selected ».
2. **Preuve à trois niveaux** : une section configurée **n'est pas** un contrôle réussi ; un contrôle du cas **n'est pas** une convergence du solveur ; les certificats natifs `execution.json` et `convergence.json` **ne remplacent pas** une qualification scientifique indépendante.
3. **Vérification périmée** : un contrôle réussi est suivi par un marqueur de révision et invalidé par les modifications apportées via l'application. Un tableau vide de diagnostics, sans appel à `validate_case`, ne valide rien. Les contrôles ne couvrant pas de maillage sont libellés « case only ».
4. **Navigation modernisée** : `Ctrl+0` pour Overview, `Ctrl+1`–`Ctrl+5` conservés, palette de commandes `Ctrl+K`, accès rapide au contrôle du cas et au maillage.
5. **Projets récents** : menu File → Open Recent et entrées dans la palette de commandes. Les chemins sont dédupliqués, bornés à 8 et validés ; seuls les projets ouverts ou sauvegardés avec succès sont enregistrés. Un maillage d'échange sans configuration ne remplace pas le cas actif.
6. **Cycle de vie** : nouveau projet = libération de l'ancien viewport, des panneaux de fichiers/maillage, des résultats sélectionnés et des traces cellulaires. Les fenêtres fermées se désabonnent des événements Qt.
7. **Carte cellules** : la superposition du lot 6 préfère uniquement le volume DAT non filtré, pour conserver les identifiants de cellules d'origine. Les opérations 3D utilisent une copie superficielle du maillage et un budget de 1 million de cellules pour la cartographie ; au-delà, seuls les centres des hotspots échantillonnés sont affichés. Pas de valeur inventée pour les cellules absentes de l'échantillon.
8. **Export** : snapshot JSON atomique du tableau de bord, avec provenance et empreintes des certificats déjà vérifiés lorsqu'ils existent. Il ne constitue pas une preuve indépendante de convergence physique.

## Contrat de vérité

- `state.workflow`: décrit la configuration et la disponibilité déclarative des fonctionnalités ; pas une certification physique.
- `Application.validate_case()`: produit des diagnostics sur le cas. Le Workbench conserve si cette vérification est actuelle ; `mesh` reste explicitement hors portée quand aucun maillage n'a été passé au validateur.
- `RunCenterState.convergence_evidence`: seul le lecteur natif **run-scoped** qui a validé les fichiers du solveur peut fournir un verdict de convergence.
- `CellTrace`: série limitée de cellules échantillonnées (top-K). Une cellule non mesurée n'est jamais étiquetée « convergée ».
- Qualification indépendante (références NACA/Ghia/SST/CHT/rayonnement) : reste une campagne scientifique distincte.

## Tests

```bash
PYTHONPATH=src/cfdx/python pytest -q \
  tests/python/test_workbench_mission_control.py \
  tests/python/test_workbench_mission_control_qt.py \
  tests/python/test_application_wiring.py \
  tests/python/test_capabilities.py
```

Tests de régression plus larges : groupes `application_*`, `run_center`, `workbench_*`, `cell_convergence`, `mesh_*`, `boundary_*`, `qualification_*`, `case_inspector`, `visualization_*`, `theme`.

Les tests dépendant de PySide6, de PyVista et d'OpenGL doivent être exécutés dans la CI graphique ; ils ne doivent pas être comptabilisés comme réussis lorsqu'ils sont ignorés.

## Références techniques

- Qt 6, Model/View et filtrage d'éléments sans dupliquer le modèle : https://doc.qt.io/qt-6/qsortfilterproxymodel.html
- OpenFOAM, function objects et journaux de surveillance à l'exécution : https://doc.openfoam.com/2306/tools/post-processing/function-objects/
- OpenFOAM, `fieldMinMax` pour la localisation d'anomalies et de leurs cellules : https://doc.openfoam.com/2312/tools/post-processing/function-objects/field/fieldMinMax/
- OpenFOAM, `runTimeControl` : https://doc.openfoam.com/2306/tools/post-processing/function-objects/utilities/runTimeControl/

Ces sources motivent l'architecture et l'ergonomie ; le lot 7 n'implémente pas leurs algorithmes numériques.

## Limites restantes

- Pas de certification GUI native Qt/OpenGL dans l'environnement local dépourvu de PySide6/PyVista.
- Le contrôle du cas actuel ne remplace pas la validation des BC natives C++ multi-champs ni la qualification sur maillage.
- Le champ spatial est encore limité au top-K jusqu'à l'ajout d'un format volumique complet par cellule.
- Les calculs longs NACA/BFS/SST et les campagnes physiques ne sont pas rejoués par ces modifications de l'interface.
