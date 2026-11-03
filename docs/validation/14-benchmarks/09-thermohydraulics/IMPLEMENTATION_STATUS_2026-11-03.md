# CFDX — TH-01/03/07/08 — état réel d'implémentation et de qualification (2026-11-03)

**Aucun des quatre benchmarks TH n'est physiquement qualifié par ce patch.**
Une validation "parfaite" ne peut pas résulter de la seule convergence itérative.
Ce document sépare les tests de contrat, les exécutions binaires réelles, et
les comparaisons indépendantes que les entrées actuelles ne permettent pas.

## Modifications réellement intégrées

- Pilote de production C++ : le runtime SST expérimental fournit `nu_t` à l'énergie;
  l'énergie calcule `lambda_eff=lambda_mol+rho*Cp*nu_t/Pr_t`. `Pr_t` doit être
  déclaré, fini et strictement positif. Il n'est pas supposé.
- Couplage expérimental SST+Boussinesq, avec itérations de point fixe en
  température/vitesse/pression/**viscosité turbulente**, bilan d'énergie natif,
  et arrêt au rouge si le solveur vitesse n'est pas convergé.
- Transport SST : les faces `EMPTY` des maillages 2D extrudés sont traitées
  à gradient nul pour `k` et `omega`; la sérialisation HDF5 conserve ces BC.
  Budget intérieur SST réglable par `sst_transport_max_iterations` entre
  10 et 2000 sans modification du seuil `1e-8`.
- Export CSV/VTU cellule et face, pressions correctement interprétées en SST,
  bilan de masse indépendant, gradients/cisaillement et sommes des flux.
- Générateur BFS de **débogage seulement**; lanceur de campagnes TH-03/TH-07
  réelles avec trois maillages et preuves SHA-256 obligatoires.
- Diagnostics réattachement à cisaillement signé, comparateur linéaire contre
  quadratique au mur (ne s'applique qu'aux parois horizontales alignées),
  processeur des séries temporelles s'il existe des données physiques réelles;
  il ne fabrique ni séries ni fluctuations turbulentes.
- Qualification fail-closed : un manque d'entrée réelle, d'historique,
  de convergence, de données de référence ou de profils sur trois grilles
  donne `INCOMPLETE_OR_FAIL`. Le mode `--require-pass` renvoie exit 2.
- CI additionnelle : validations des métadonnées, C++/Python, vrais smoke
  HDF5 SST+énergie, et vrai SST sur marche de débogage.

## Mesures CFDX réellement exécutées

| Essai | État | Résultats vérifiés | Interprétation |
|---|---|---|---|
| BFS 192 cellules, Re_h 36 000, entrée uniforme | `converged=true`, 145 itérations | Résidu momentum `5.64976e-12`, continuité `8.72791e-13`, déséquilibre massique global `4.04903e-14`, local `4.32298e-10`; x_r/h au premier ordre `1.76882`, reconstruction quadratique `2.26478` | Physiquement faux vs Driver–Seegmiller `6.26 +/-0.10` |
| BFS 402 cellules, même hypothèse d'entrée | `converged=true`, 220 itérations | Résidu momentum `1.56581e-11`, déséquilibre massique global `7.21960e-14`; x_r/h premier ordre `2.22988`, quadratique `2.81068` | Sensible au maillage et au gradient pariétal; toujours faux |
| BFS 688 cellules | Arrêt d'exécution, 45 itérations seulement | Résidu momentum `6.35348e-4`, continuité `1.80855e-4`, déséquilibre global `2.96148e-5` | **Non convergé**; pas de x_r/h fin ni GCI démontrable |
| Conduction stationnaire SST, 15 cellules, `Pr_t=0.40` | `converged=true` | Flux pariétal froid `~2.76432 W`, bilan énergétique relatif `1.35069e-13` | La fermeture de conductivité est active ; ne valide pas une convection turbulente |
| Même cas, `Pr_t=1.80` | `converged=true` | Flux pariétal froid `~0.93706 W`, bilan relatif `1.43457e-14` | Sensibilité réelle au Pr_t ; validation constitutive encore requise |
| Boussinesq/SST, température uniforme 300 K | `converged=true` | Température uniforme, flux nul, bilan de masse nul | Équilibre trivial, pas de convection mixte |
| Boussinesq/SST, parois 300/310 K | `converged=false`, 120 itérations | Déséquilibre massique relatif `0.566609` | Porte de convergence correctement rouge; couplage non qualifié |

**Incompatibilité explicite :** le maillage de débogage BFS utilise `Re_h=36 000`
alors que le contrat TH-03 du dépôt exige `Re_h=37 500`; ils ne
constituent pas un scénario physiquement identique. Même la source publique
NASA/OpenFOAM utilisant `Re_h=36 000` requiert les profils turbulents
mesurés à l'entrée, des parois résolues, et des stations de validation.
Le générateur ne crée volontairement pas de preuves d'identité expérimentale.

## Couverture du cahier des charges et blocages réels

| Benchmark | Acquisition de sorties | Capacité de production | Blocage qualification |
|---|---|---|---|
| TH-01 DNS channel | Exports natifs et comparateur sont présents | **NON** : transport DNS/LES instationnaire 3D périodique et statistiques moyennées absents | Résoudre le 3D instationnaire et comparer 180/395/590 avec DNS épinglées |
| TH-03 BFS RANS SST | Solveur SST expérimental + réattachement | **PARTIELLE** : HDF5 réel et BC SST + transport de production | Profil de couche limite mesuré, paroi `y+`, 3 maillages convergés, `Cp`, `Cf`, U aux 5 stations, GCI |
| TH-07 convection mixte | Conductivité SST/Pr_t, énergie, flottabilité, bilans indépendants | **PARTIELLE** : échec encore observable pour un gradient 10 K | Convergence flottabilité à Gr non nul, profils DNS Re_b4328/Pr0.71, 3 maillages chacun |
| TH-08 jonction T | Exports et post-traitement de vrais historiques d'échantillonnage | **NON** : transport thermique et turbulence URANS/LES instationnaires non branchés | Séries stationnaires en temps physique, RMS, PSD, CI temporelles, études dt et 3 maillages |

Le comparateur donne au mieux un **pré-contrôle** de preuves, et non la
validation expérimentale certifiée. Une précision parfaite ne doit pas
être promise : l'erreur de discrétisation, de modèle turbulent et
l'incertitude expérimentale sont distinctes.

## Tests de non-régression

- C++ : `test_thermohydraulic_observables`, `test_production_sst_thermal_coupling`,
  `test_thermal_driver_vv`, `test_production_buoyancy_outer_loop_vv` : **4/4 PASS**.
- Python ciblé : **39/39 PASS**, plus un E2E BFS natif optionnel
  `tests/python/test_th03_native_sst_e2e.py` : **1/1 PASS**.
- Qualificateur strict sans preuves : statut `INCOMPLETE_OR_FAIL`, **exit 2 attendu**.

### Reproduire l'essai BFS natif (débogage, 192 cellules)

```bash
export PYTHONPATH=src/cfdx/python
python scripts/generate_th03_bfs_case.py --output build/th03_debug.cfdx.h5
build/th-preflight/cfdx_production_solver \
  --mesh build/th03_debug.cfdx.h5 --output-dir build/th03_debug_out \
  --iterations 250 --experimental-sst --benchmark-export
python scripts/compare_th03_wall_gradients.py \
  --face-csv build/th03_debug_out/benchmark_observables/face_fluxes.csv \
  --cell-csv build/th03_debug_out/benchmark_observables/cell_fields.csv \
  --density 1.2 --kinematic-viscosity 1.5592777777777776e-05 \
  --bulk-speed 44.2 --step-height 0.0127 \
  --output build/th03_debug_wall_fit.json
```

Pour la qualification, remplacer d'abord le maillage de débogage par un
jeu traçable conforme au cas expérimental, puis exécuter `scripts/run_thermohydraulic_qualification.sh`.
