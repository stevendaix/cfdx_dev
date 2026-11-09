# CFDX Workbench — Lot 4 : boundaries, mesh health et preuves d'exécution

Date : 2026-11-08. Base : **ZIP exact `CFDX_Workbench_lot3_noir_moderne.zip`**.

## Principes invariables

- Le fond noir OLED du lot 3 reste le thème par défaut.
- Pas de nombres physiques fabriqués ou de seuils du solveur assouplis.
- Un processus avec code de sortie 0 n'est **pas** une preuve de convergence.
- Les événements issus des threads ne doivent pas toucher directement les widgets Qt.
- Un résultat sans donnée, non calculable ou non mesuré ne doit pas être présenté comme zéro.

## 1. Éditeur de conditions limites à champs multiples

Implémentation : `boundary_editor.py`, `boundary_matrix.py`, `boundary_setup.py`, `setup_panel.py`, `execution.py`.

**Ce qui fonctionne :** édition d'un champ par patch, mise en scène de plusieurs champs sans perte des anciens, valeurs finies, unités visibles, guide de choix, navigation par nom de patch, tableau des champs, suppression d'un champ et retour exact au contrat Python mono-champ.

**Sécurité au lancement :** la même validation de définition est appelée dans `ExecutionController.start()` et `restart()`, même si le bouton « Validate » n'a pas été utilisé. Une définition multi-champs non prise en charge est rejetée avant création d'un processus. Changement destructif de type et suppression demandent confirmation côté Qt.

**Limite native à ne pas masquer :** le modèle de configuration Python est sérialisé dans le dataset `/case/config`; le solveur C++ de production lit l'attribut HDF5 **`case_setup_json`**, dont les entrées sont `boundary_conditions` avec `patch_name`, `type`, `value_type`, `velocity_vector`, `pressure`, etc. Le bridge bidirectionnel et les cas multi-champs pour vitesse/pression/énergie/SST ne sont **pas encore certifiés**. « Contrat mono-champ accepté dans l'éditeur » ne signifie pas « champ effectivement appliqué par le C++ ». On ne prétend pas avoir complété cette intégration.

**Priorité d'intégration suivante :** établir une table de correspondance exhaustive Python/`CaseSetup`, ajouter un export natif transactionnel qui préserve les conditions non modifiées, puis vérifier dans le C++ que chaque condition de chaque patch influence l'opérateur réel (dont SST, chaleur, gradients et reflux). Interdire le lancement si la configuration sauvegardée et lue par le C++ diffère de celle affichée.

## 2. Mesh Health Dashboard

Implémentation : `mesh_health.py`, `case_file_panel.py`, `volume_skewness_view.py`.

- Ouverture HDF5 : lecture des métadonnées et contrôles existants, sans forcer le chargement de tous les datasets.
- Calcul qualité : appels à `analyze_cfdx_volume_skewness` (mêmes conventions que le code existant), dans un QRunnable du pool Qt.
- Résumé : nombre de cellules effectivement évaluées, valeurs indisponibles, min/médiane/p95/max, histogrammes, seuils d'alerte réglables, liste des pires cellules.
- Viewport : coupes selon X/Y/Z, non-orthogonalité ou skewness, surbrillance des cellules au-dessus du seuil, sélection des cellules critiques et de leurs voisines.
- Exports : JSON structuré avec identité du fichier, contrôle de source inchangée. L'alerte graphique n'est pas un critère de rejet ou de convergence du solveur.

**Limite** : le chargement des données de maillage volumique est borné par le budget actuel du lecteur existant (limites par taille/cellules) ; ne pas annoncer de rendu 10–50 millions de cellules. Pour les très gros maillages, conserver l'inspection de métadonnées et ajouter un futur pipeline streaming/partitionnement natif.

## 3. Cockpit de convergence avec preuves natives

Implémentation : `convergence_evidence.py`, `execution.py`, `application/run_center.py`, `run_center_panel.py`.

- Capture des identités de `execution.json` et `convergence.json` **avant** le calcul, puis refus des fichiers inchangés, manquants, liens symboliques, JSON corrompus, clés dupliquées, NaN/infinités et dépassement de budget (64 MiB/fichier).
- Confrontation du code de sortie du processus, du verdict natif global (`execution.json`), du verdict écoulement (`convergence.json`), des versions de schéma, du compteur d'itérations et du dernier point d'historique.
- Le résultat du solveur multi-physique (énergie/turbulence et étapes éventuelles) est distinct du seul résultat momentum/pression.
- Affichage des métriques finales effectivement présentes, regroupées par équation/conservation/linéaire/mise à jour, provenance et SHA-256 de chaque document.
- Export JSON du certificat lu, en indiquant explicitement que la qualification physique indépendante n'a pas été effectuée.
- En l'absence de preuve valide, `CONVERGED` n'est jamais déduit de la disparition du processus. Un retour 0 sans preuve est `STOPPED` non vérifié ; un code 1 certifié non convergent est `STOPPED` ; un échec de processus est `FAILED`.

**Limite de preuve :** les deux fichiers natifs n'embarquent pas un ID de campagne généré avant lancement ni un lien cryptographique au commit, à la géométrie, aux critères ou aux références expérimentales. Les SHA-256 assurent l'empreinte des fichiers observés, pas une attestation externe. Le verdict natif n'est ni un GCI ni une validation physique NASA/DNS.

## 4. Vérification et commande

Depuis la racine du dépôt :

```bash
PYTHONPATH=src/cfdx/python python -m pytest -q \
  tests/python/test_workbench_lot4_evidence.py \
  tests/python/test_workbench_lot4_boundary_matrix.py \
  tests/python/test_workbench_lot4_mesh_health.py \
  tests/python/test_workbench_lot4_gui.py \
  tests/python/test_execution.py \
  tests/python/test_application_e2e.py
```

La CI `.github/workflows/cfdx-gui.yml` reprend les quatre nouveaux fichiers de test dans sa campagne xvfb + Qt/PyVista.

Tests de sûreté : certificat manquant/périmé/non fini, divergence entre étapes, sauvegarde multi-champs et suppression sans pertes, HDF5 réel, source maillage modifiée, exception pré-lancement et régression E2E.

Les tests Qt directs nécessitent PySide6, pyvistaqt et un environnement OpenGL/Xvfb ; en leur absence, ils sont comptés *ignorés*, jamais validés. Le solveur C++ de production n'est pas modifié par ce lot et ses qualifications longues n'ont pas été exécutées.

## 5. Bibliographie utile

- Qt — QThreadPool / QRunnable : https://doc.qt.io/qt-6/qthreadpool.html
- OpenFOAM — Boundary conditions (conditions par champ et familles) : https://www.openfoam.com/documentation/user-guide/5-models-and-physical-properties/5.1-boundary-conditions
- PETSc — SNES monitoring et convergence : https://petsc.org/main/manual/snes/
- CFDX natif dans ce dépôt : `apps/production_solver/main.cpp`, `src/cfdx/io/hdf5/case_hdf5_io.cpp`, `src/cfdx/io/runtime/convergence_history.h`.

Le lot ne copie pas les méthodes PETSc/OpenFOAM, il adopte leurs **contrats et distinctions scientifiques vérifiables**.
