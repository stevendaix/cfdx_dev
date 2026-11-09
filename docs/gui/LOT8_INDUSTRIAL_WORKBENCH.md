# Lot 8 — Cartographie complète, BC natives et outils industriels

Ce lot étend le Workbench noir du lot 7 sans modifier les tolérances de convergence du solveur.

## 1. Cartographie spatiale native (`Run`)

```bash
cfdx_production_solver --mesh case.cfdx.h5 --output-dir results \
    --iterations 500 --cell-hotspots --cell-hotspots-every 5 --cell-map-full
```

- Le traceur C++ existant calcule les défauts de momentum reconstitués après **acceptation** de chaque itération non linéaire et le déséquilibre local des flux massiques.
- `--cell-map-full` crée `results/cell_convergence_map.bin`, en conservant l'ancien `cell_convergence.jsonl` top-K. Le fichier binaire (`CFDXMAP1`) comprend toutes les cellules, un bit de validité par cellule et une borne stricte de **512 Mio**.
- Les valeurs affichées du momentum sont `log10(local momentum/RHS / référence fixe)`; plage de couleurs fixe **[-2,4]** pour comparer les itérations. La référence est un indicateur de résidu, **pas** une erreur locale estimée de solution. Valeur manquante = gris neutre, jamais bleu arbitraire.
- Choisir `Run → Cell hotspots` puis itération/métrique. Pour peindre les cellules (et non les points), charger d'abord un DAT volumique qui conserve l'ordre original des cellules, dans la limite du rendu de la vue. Le fichier d'origine est vérifié.
- La carte binaire se lit en streaming/mmap avec `cfdx.cell_field_map.CellFieldMap`. La couleur de l'historique du volume n'est pas un verdict de convergence globale.

## 2. Conditions limites distinctes par champ (`Setup`)

- Nouveau contrat natif `CaseSetup.boundary_conditions[].field_conditions` pour `velocity`, `pressure`, `temperature`; `k` et `omega` restent dans le contrat SST existant `sst_scalar_conditions`.
- Le C++ applique les paramètres natifs vitesse/pression et thermique. Les modes non supportés sont refusés ; aucune valeur ne se propage implicitement d'un champ à l'autre.
- Dans `Setup`, l'action d'export écrit **un nouveau fichier HDF5 natif** avec `case_setup_json`, en conservant les paramètres physiques et numériques existants. Elle refuse d'écraser l'original, de déduire une physique depuis `/case/config`, ou de faire passer un export pour un lancement qualifié.
- Les modèles SST restent expérimentaux là où indiqués par le solveur ; ce pont n'établit pas une qualification des lois de paroi.

## 3. Exploration HDF5 progressive (`Analyze`)

- `Large-field inspector` : catalogue de champs natifs, statistiques numériques et aperçus bornés par **8 Mio** de valeurs décompressées et limites de métadonnées.
- Lecture dans un worker Qt séparé ; refus des liens HDF5 non internes, champs virtuels, chunks surdimensionnés, formes inconnues.
- Ce service ne remplace **pas** encore un renderer VTK partitionné de 10/50 millions de cellules. Le volume 3D non partitionné reste limité par les budgets du Workbench.

## 4. Comparaison contrôlée (`Analyze`)

- `Compare runs` : pour un champ scalaire HDF5 et des géométries avec hash de topologie strictement identiques, même temps, même unité et même nombre de cellules, calcul des différences L2/RMS/L∞/L2 relatif par blocs sans interpolation.
- `Solver metrics` : comparaison de deux paires archivées `execution.json` + `convergence.json`, avec deltas sur les indicateurs communs. Ajoute les **forces dimensionnelles** de `wall_forces.json` lorsque les conventions et les références coïncident.
- Les métriques ne certifient **pas** que les deux calculs utilisent la même physique ou le même maillage ; les forces ne sont pas transformées en coefficients sans références. `elapsed_wall_seconds` reste non disponible, car absent du format de sortie natif actuel.

## 5. Centre de qualification et CI (`Validate`)

- Un nouveau profil `bfs` est ajouté au sélecteur CI, en plus des profils NACA, Ghia, SST, thermique, etc.
- Le workflow ciblé produit `qualification/cfdx-qualification-evidence.json` lié au SHA du commit.
- Le Workbench peut ouvrir ce manifeste ; chaque cas garde deux états distincts : statut **déclaré** du registre et **preuve physique** propre au commit.
- Une prétendue preuve PASS sans critères numériques physiques, fichiers et SHA256, ou portant sur un autre commit reste **UNVERIFIED**. Un CTest vert seul ne prouve pas la qualification physique.

## Tests et limites

Les tests du lot doivent couvrir : lectures HDF5 sûres et comparaisons incompatibles, roundtrip C++ des conditions multi-champs, émission binaire sur le solveur natif, absence de certification fictive, régression des interfaces et du registre. Les cas NACA/BFS/Ghia/SST/thermiques/radiatifs longs et le rendu natif Qt/OpenGL ne sont pas implicitement validés par ces tests de contrat.
