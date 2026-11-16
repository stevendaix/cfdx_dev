# CFDX V36 — réparation du solveur BFS V6 et qualification de la trace native

Date : 11 novembre 2026. Préimage **exacte** : `cfdx_dev-physics-bfs-v6-hilo-rhie-chow-production.zip`, SHA-256 `c29fc9c4aa92057c99aecd8fc8bc087124eb12174ba388c88063c34b3a32bff6`.

## Portée

La fusion V35 contenait des fragments contradictoires de `steady_incompressible_solver.h` et de `apps/production_solver/main.cpp` : structure `IncompressibleSolverControls` amputée, redéclarations PCD et arguments Krylov incohérents, portions du cycle de rollback/acceptation SST absentes, et redéclarations d'options COUPLED. Ce n'était pas un échec de convergence mesuré mais un **échec de compilation**, antérieur à toute validation BFS.

V36 remplace les **deux unités de production concernées** par leur préimage complète fonctionnelle V34, maintient les modules et fichiers spécifiques à la base BFS V6, puis réintègre sélectivement les ajouts réellement utiles absents de V34. `pressure_velocity_algorithms.h` perd seulement des déclarations d'options répétées.

### Réparations conservatrices

- Restauration de la structure de contrôle SST et des cycles PTC V31–V32, du choix JFNK V22–V34, des protections de restauration d'état et des opérateurs PCD/Schur de production.
- Conservation du mode pressure high/low (`pressure_high_low_production`), du flux Rhie–Chow autoritaire, des checkpoints DAT6 et de la capture CSR, sans les promouvoir.
- Raccordement du diagnostic cellulaire (JSONL et carte binaire) depuis les matrices momentum et les flux de masse FV effectivement réassemblés, *uniquement après* acceptation de l'itération extérieure. Les échelles de normalisation sont gelées au premier état accepté ; les essais rejetés ne créent pas d'amélioration artificielle.
- Correction du pied de trace : compte exactement le dernier état observé et accepté, et non la dernière tentative d'itération retournée par le solveur.
- Conservation des conditions limites explicites, distinctes pour `velocity` et `pressure`, avec refus des vecteurs mal dimensionnés/non finis, des gradients Dirichlet contradictoires et des rôles physiques incompatibles. Les champs `temperature` et SST restent gérés par leurs opérateurs distincts.
- Ajout du contrôle `test_cell_hotspots_native_hdf5_v36` sans `meshio` : exécute le même solveur avec diagnostic OFF/ON sur un maillage natif HDF5, compare **à l'identique** les fichiers de convergence, vérifie le pied des cartes, le caractère non certifiant de la trace et le rejet d'une fréquence d'échantillonnage nulle.

## Mathématiques et limites

Les cartes de cellules exposent les défauts d'équations réassemblées `|A_u(U,p,k,omega) U-b_u|` et les défauts locaux `|sum_f mdot_f|`. Elles ne sont **pas** la norme globale indépendante du système complet

`F(U,p,k,omega) = (R_U, R_continuite, R_k, R_omega)`.

La direction d'un solveur de Newton–Krylov pour ce système doit vérifier séparément `J(x)v ≈ [F(x+eps*v)-F(x)]/eps` avec réassemblage simultané de `nu_t`, `G`, contributions de paroi et Rhie–Chow. Les six variantes PCD, LSC, BFBt, BlockLocal, SIMPLE-Schur et SIMPLEC-Schur doivent partager le **même** opérateur `F` et la **même** norme physique, et jouer le rôle de préconditionneurs uniquement. V36 restaure l'architecture existante mais ne prétend **pas** avoir réalisé cette campagne V35/V36 globale.

Seuil de convergence SST conservé : `1e-8`. Rien n'est activé par défaut. Ni le plafond SST-2003, ni les schémas de transport, ni les critères de conservation ne sont assouplis.

## Reproduction (après application du patch cumulatif à l'archive V6 uniquement)

```bash
cmake -S . -B build-v36 -DCMAKE_BUILD_TYPE=Release -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON -DCFDX_BUILD_TESTS=ON
cmake --build build-v36 --target cfdx_production_solver test_production_boundary_contract \
  test_cell_convergence_tracker test_cell_field_map_writer \
  test_sst_jfnk_v22 test_sst_spatial_block_preconditioner_v24 \
  test_sst_branch_preconditioner_v34 -j 2
ctest --test-dir build-v36 -R 'test_cell_hotspots_native_hdf5_v36|test_production_boundary_contract|test_cell_convergence_tracker|test_cell_field_map_writer|test_sst_jfnk_v22|test_sst_spatial_block_preconditioner_v24|test_sst_branch_preconditioner_v34' --output-on-failure
```

Le drapeau `CFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON` compile `main.cpp` en `-O0` pour raccourcir les essais : **ce n'est pas une qualification de performance Release optimisée**. Ne jamais l'utiliser pour comparer des temps de convergence industriels.

Pour activer explicitement les cartes : `cfdx_production_solver --mesh case.h5 --output-dir out --iterations 20 --cell-hotspots --cell-map-full --cell-hotspots-every 1`.

## Preuves de cette livraison

- `cfdx_production_solver` compilé et lié.
- **17/17 CTest ciblés** (incluant le maillage HDF5 quick, les tests SST V22/V24/V31–V34, Rhie–Chow, la BC contract, le suivi cellulaire et le contrat de paroi V6) réussis dans le build local avec MPI OFF.
- Un cas HDF5 quick BFS à champ uniforme nul, deux itérations, sortie complète : `cell_convergence.jsonl` et `cell_convergence_map.bin` présents, `accepted_iterations=2`, deux frames, `convergence_certificate=false`.
- Deux tests historiques `.msh` de production tentés et **échoués avant résolution**, faute du module Python `meshio` dans l'environnement. Ne pas les déclarer PASS ou les ignorer dans le bilan. `test_bfs_quick` n'a pas été exécuté : cible lourde non construite. Un test C++ COUPLED lourd a atteint le budget de compilation local.
- **Non démontrés** : convergence SST BFS 712/20 540 cellules, six Schur exécutés sur cas réel, réattachement NASA `x_r/h=6.26 ± 0.10`, certification OpenFOAM native, CI intégrale, DebugSanitizers intégral, CPU multicoeur/MPI et GPU, performances.

Aucune modification du dépôt GitHub et aucun push. Statut : **fusion V36 techniquement compilable, tests ciblés validés, qualification physique BFS encore ouverte**.
