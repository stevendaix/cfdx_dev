# CFDX Workbench — Lot 5 : diagnostic spatial de convergence par cellule

Date : 8 novembre 2026. Base du lot : `CFDX_Workbench_LOT4_COMPLET_TESTE.zip`.

## But et limites scientifiques

L'outil répond à la question **« quelle cellule empêche la convergence ? »** sans déduire un résultat local d'un simple résidu global. Le solveur de production instrumente, en option, son opérateur de momentum réassemblé et le déséquilibre conservatif des flux de masse par cellule. Chaque cellule est identifiée par son **ID d'origine** et le centre géométrique du maillage utilisé. L'état *difficile* n'est pas une preuve causale que la cellule *empêche* seule la convergence : il signale une zone à investiguer.

Les trois dimensions de la mesure sont volontairement distinctes :

1. `momentum_abs[d] = |sum_j A_d[i,j] U_d[j] - b_d[i]|` évalué sur l'état non linéaire réassemblé ; **pas** un résidu du solveur linéaire de pression.
2. `momentum_scaled = max_d(momentum_abs[d] / S_d)`, avec `S_d = max_i |b_d[i]|` borné par le plancher numérique existant et **figé lors de la première itération acceptée**. Le schéma ne change aucun seuil de convergence ; la normalisation n'est qu'un indicateur de classement spatial.
3. `mass_imbalance_abs = |sum_faces signe(owner/neighbour) mdot_face|`, flux natif Rhie–Chow / conservation du solveur, **en unités de débit massique**, sans normalisation locale implicite.

`nondecrease_streak` compte **pour une même cellule** les transitions entre itérations NON LINÉAIRES **acceptées** où le score reste au-dessus du seuil de référence et ne diminue pas de plus de 2 %. La troisième transition consécutive déclenche le signal indicatif `>=3`. Le compteur est réinitialisé si le schéma convectif change pendant une continuation documentée. Des changements d'autres paramètres de l'opérateur peuvent néanmoins réduire la comparabilité : **ne pas interpréter le compteur comme une preuve mathématique de stagnation**. La hausse d'un résidu compte aussi comme « non-décroissance » ; ce n'est pas nécessairement un plateau.

Le traceur **n'altère ni les équations, ni les critères physiques, ni les paramètres des solveurs linéaires**. Les tentatives rejetées par la transaction non linéaire ne sont pas comptabilisées.

## Utilisation native

```bash
cfdx_production_solver --mesh /absolute/path/case.cfdx.h5 \
    --output-dir /absolute/path/results --iterations 500 \
    --cell-hotspots --cell-hotspots-every 5
```

Ouvrir **Run → Cell hotspots** dans le Workbench noir et cocher « Cell hotspots (opt-in) » **avant** de lancer un nouveau calcul avec `cfdx_production_solver`. Le panneau lit le fichier `results/cell_convergence.jsonl` en tâche Qt séparée, propose trois classements (*momentum / RHS*, *mass imbalance*, *non-decrease streak*) et une période de recalage de 2,5 s. Sélectionner un ID dans le tableau montre ses coordonnées, son historique réellement échantillonné et le centre de cette cellule dans le viewport ; la caméra se déplace vers ce point. Le maillage source HDF5 (ou DAT renvoyant à cette géométrie) doit d'abord être ouvert dans la même interface.

Si un jeu de données volumique avec le même nombre de cellules est affiché (par exemple un checkpoint DAT ou un résultat volumique), le Workbench peint maintenant les cellules échantillonnées elles-mêmes avec une échelle **bleu → rouge** et une barre de couleurs. Les cellules sans mesure directe restent en gris neutre ; le Workbench ne fabrique jamais de champ complet à partir d'un top-K partiel. En l'absence de géométrie volumique compatible, l'overlay retombe sur les points aux centres des cellules.

Pour examiner les traces sans GUI :

```bash
PYTHONPATH=src/cfdx/python python -m cfdx.cell_convergence \
    /absolute/path/results/cell_convergence.jsonl --metric persistence --top 24
PYTHONPATH=src/cfdx/python python -m cfdx.cell_convergence \
    /absolute/path/results/cell_convergence.jsonl --json > hotspots_summary.json
```

La trace native contient un en-tête, des snapshots des **24 pires cellules** dans chacune de trois catégories (union dédupliquée, au plus 72 cellules par snapshot par défaut), puis un marqueur de fin si l'exécution a terminé l'écriture. Le suivi de persistance s'effectue **sur toutes les cellules à chaque itération acceptée** ; les valeurs détaillées ne sont stockées que pour les échantillons exportés. Les trous de la série temporelle d'une cellule signifient « **non enregistré** », pas « convergence locale prouvée ». Les fichiers n'offrent donc **pas** une carte complète du champ résiduel à chaque itération ; le viewer ne colorise **que** les points observés, sans interpolation trompeuse.

## Protection mémoire, identité et sécurité

- Diagnostic désactivé par défaut : aucune matrice ni historique spatial supplémentaire sans option explicite.
- Budget C++ : maximum **5 millions de cellules** pour l'activation de ce mode et au plus 64 échantillons par groupe (24 retenus en production). Les allocations de travail sont `O(n_cells)` ; le flux disque est borné par la fréquence et le top-K. **5 millions de cellules ne signifient pas 5 millions de points affichables en mémoire Qt.**
- Le lecteur Python accepte au plus **32 Mio** et **128 Kio par ligne**, refuse les objets non réguliers, symlinks Linux et JSON non fini / dupliqué / incohérent ; une ligne partielle en fin d'écriture n'est jamais présentée comme valide. Une nouvelle exécution efface le fichier précédent, en refusant les liens symboliques pour la sortie native.
- Le viewport compare le chemin absolu du maillage source et le nombre de cellules lorsqu'il dispose d'un objet DAT. **Un chemin identique n'est pas une preuve cryptographique d'un maillage inchangé** ; la vérification complète d'empreinte géométrique est à ajouter avant d'affirmer un contrat d'identité fort.
- `execution.json` et `convergence.json` restent les seuls artefacts utilisés pour le verdict de convergence natif ; ce diagnostic ne fournit **jamais** un certificat de convergence ou de qualification physique.
- Une exécution interrompue avant la première itération acceptée peut contenir un marqueur de fin à zéro itération : il n'implique aucun succès physique.
- Les champs turbulents `k`, `ω`, énergie et radiation ne sont **pas encore** instrumentés spatialement dans cette livraison : aucune affirmation de diagnostic SST/thermique cellulaire complet.

## Tests du lot

1. `test_cell_convergence_tracker.cpp` : IDs, scores re-normalisés sur base fixe, persistance de la même cellule sur états acceptés, changement d'opérateur, fréquence d'échantillonnage, arrêt sans itération et arrêt entre snapshots.
2. `test_cell_convergence.py` : lecteur avec bornes, conflits entre en-tête/frames, NaN, cellules dupliquées, limites de taille, JSON partiel, classement, historique troué et absence de certificat.
3. `test_cell_hotspots_gui.py` : création Qt, lecture asynchrone, sélection et historique, opt-in CLI avec refus du mauvais exécutable. Exécuter dans la CI Xvfb/PySide6 ; **ignoré localement** quand Qt absent.
4. `test_cell_hotspots_native_e2e.py` : CI de bout en bout du solveur vers le lecteur, sur le maillage CTest ; nécessite `meshio` dans l'environnement de test.
5. Smoke manuel indépendant : génération d'un cube HDF5 1 cellule via le writer C++, `cfdx_production_solver --cell-hotspots --cell-hotspots-every 2`, lecture du JSONL par Python. Observé : 2 itérations, 2 snapshots, cellule #0 au centre (0,5;0,5;0,5), résidu nul. Ce smoke **ne démontre pas** une stagnation réelle ni la précision aérodynamique.

## Prochaines extensions justifiées

- Fournir la reprojection d'une trace entre maillages avec une correspondance géométrique **validée** ; ne jamais réutiliser des IDs sur un maillage adapté/remanié.
- Étendre le résidu spatial natif aux équations de pression, `k`, `ω`, énergie, rayonnement et VOF, avec unités, normalisations, bilans et certifications indépendantes par équation.
- Ajouter la localisation des sources locales du défaut : convectif, diffusif, gradient de pression, condition limite, correction non orthogonale ; cartes et histogrammes pour les voisinages des cellules défaillantes.
- Campagne physique N64/N128/N256 et SST/thermique : corréler les hotspots aux critères et à la sensibilité au maillage, puis comparer les temps jusqu'à convergence sur les mêmes critères.

## Références

- OpenFOAM documentation, residual controls et définition algébrique des résidus : https://doc.openfoam.com/2306/tools/processing/numerics/solvers/residuals/monitoring/
- PETSc SNES, évaluation du résidu non linéaire réel `F(x)` : https://petsc.org/main/manual/snes/
- ParaView, pipelines et filtres Threshold : https://docs.paraview.org/en/latest/UsersGuide/filteringData.html
- Qt, exécution concurrente avec QThreadPool : https://doc.qt.io/qt-6/qthreadpool.html

Les références fournissent une logique de diagnostic et d'architecture, **pas** une preuve de conformité physique de la nouvelle instrumentation CFDX.
