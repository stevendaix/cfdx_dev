# CFDX BFS V11 - Pression, Schur/PCD, SST et statut physique

Date de travail : 9 novembre 2026. Source : V10 complète, elle-même issue de l'archive utilisateur `cfdx_dev-master (2)(1).zip` (SHA-256 `c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c`).

## 1. Décisions de qualification

Cette livraison ne déclare **aucune convergence physique complète** sur Driver-Seegmiller. Les trois niveaux demeurent distincts : exécution, convergence de chaque équation sur état retourné, validation des observables contre les données expérimentales. Aucune tolérance physique ou linéaire n'a été abaissée, aucun code de retour d'échec n'est transformé en PASS. Les sorties `TIMEOUT_INCOMPLETE` ne sont pas des preuves de divergence numérique.

Références :
- OpenFOAM v2606 : https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/ ; `Re_h=36000`, `U=44.2 m/s`, `h=0.0127 m`, et `x_r/h=6.26 +/- 0.10`.
- NASA TMR : https://tmbwg.github.io/turbmodels/backstep_val.html ; profils x/h=1,4,6,10, et mesures Cf/Cp/contraintes turbulentes.
- OpenFOAM natif : https://gitlab.com/openfoam/core/openfoam/-/tree/OpenFOAM-v2606/tutorials/incompressible/simpleFoam/backwardFacingStep2D . Le `blockMesh` natif n'est PAS installé ici. Les HDF5 utilisés sont une reconstruction six blocs, pas une sortie bit-à-bit de blockMesh.

## 2. Défauts et correctifs V11

### Pression / Rhie-Chow

A/B réel en production sur 712 cellules et 2907 faces, une itération externe :
- SIMPLE V10 par défaut : correcteur arrêté avec residu relatif `4.24078912e-11`, seuil inchangé `1e-11`.
- `v18_sparse_flux_pressure_jacobian=true` + `v19_full_pressure_flux_step=true` + `pressure_high_low_production=true` : résidu CSR direct `2.390865895e-11`, défaut vrai aux faces initial `5.692669085e-11`, puis vrai défaut final `1.079505469e-16` après un raffinement compensé. Cela prouve la fermeture **du correcteur en flux aux faces seulement**.
- Avec le SST nominal de 800, le cas SIMPLE reste non convergé : la résolution turbulente ne termine pas dans le budget court de 18 s. Il est interdit d'en déduire le succès U,p,k,omega.

### COUPLED Schur/PCD

L'ancienne fonction `nested_schur_pressure_solve_usable` acceptait `CONVERGED` sans contrôler le tableau retourné. Le correctif impose, pour CONVERGED **et** MAX_ITER_REACHED, taille attendue, composantes finies, et residus rapportés finis/non négatifs. Aucun fallback silencieux.

`CFDX_BFS_PCD_INNER_TRACE=1` trace chaque application Fp et, lorsque utilisé, Kp. Les mesures incluent le défaut vrai sur la matrice **d'origine**, non seulement le résidu rapporté par un solveur de matrice équilibrée. Le pilote Python agrège nombre d'applications, pire défaut réel et actions rejetées ; toute trace malformée bloque sa qualification tout en conservant le journal brut.

Quatre commandes de recherche expérimentales sont maintenant exposées dans les `raw_settings` HDF5 du solveur de production ; chacune est rejetée hors de COUPLED/PCD :
- `pcd_fp_row_equilibration_experimental`
- `pcd_fp_quality_guard_experimental`
- `pcd_fp_second_pass_experimental`
- `pcd_lumped_mass_commutator_experimental`

La première matrice PCD, six blocs 712 cellules : 208 Fp, 208 défauts vrais relatifs Fp au-dessus de `1e-3`, médiane de `1`, pire `1`, et 0 appel Kp en `LegacySingleConvectionDiffusionInverse` (normal). FGMRES extérieur : `52` itérations, residu `1.05585e-2`, statut `NOT_CONVERGED`.

Expériences de sélection indépendantes, limite temporelle 12 secondes pour chacune :
- Fp `quality_guard` : 4/4 actions Fp refusées (défaut > 1e-3) ; échec rapide au premier FGMRES, 0 itération. Ce comportement est **fail-closed**, pas une convergence.
- Fp `row_equilibration` : 40 actions Fp, pire défaut vrai relatif 7.69829 (amplification !) ; FGMRES extérieur arrêté à 10 itérations, résidu ~1.11502e-2. L'équilibrage ne corrige donc pas la cause sur cette matrice.
- Fp `second_pass` : action très coûteuse ; expiration du budget avec défaut extérieur encore ~1.06119e-2 lors de la tentative observée. Pas de succès certifié.
- `lumped_mass_commutator` : 971 appels Kp observés, avec défauts indépendants typiquement de l'ordre de 1e-11 ; pas de Fp inverse en ce mode. Le budget temporel s'épuise avant la fin de FGMRES ; **pas de qualification**.

Ces ordres de grandeur concernent la première itération externe d'une **seule grille de diagnostic**, ils ne constituent ni classement définitif des Schur ni recommandation universelle.

### SST, 800 itérations internes

A/B sur COUPLED BlockLocal, mêmes données et mêmes tolérances :
- Méthode héritée : `k_residual=0.00555927`, `omega_residual=1.243515` après 800 itérations (échec).
- Méthode optionnelle `sst_omega_newton_sink=true`, `sst_sequential_omega_first=true` : certification interne retournée à itération 213 ; `k_residual=8.198997e-13`, `omega_residual=3.025602e-10`. Après réassemblage, momentum global `U=0.527537`, continuité `0.00176003`. **La solution globale reste NON CONVERGÉE.**

`CFDX_SST_DIAGNOSTICS=1` ajoute un diagnostic indépendant de 5 conditions, sans modifier la logique de rejet : residuals internes, acceptation du pas, finitude, variation de champs et défauts des équations *réévaluées sur l'état retourné*. Une mauvaise fermeture SST sur SIMPLE avec flux compensés ne doit pas être masquée par un bon résidu k ou omega de l'ancienne matrice.

### Fonctions de paroi

Le code de référence calcule déjà la viscosité Spalding, `omegaWallFunction` et un proxy de production G ; un audit a trouvé des différences importantes, et le composant V8 conserve explicitement `openfoam_equivalent=false`. Il manque la pondération commune des cellules aux coins, la mise à jour couplée du champ de patch nut et de G et la certification du budget SST sous la même physique. **Aucune équivalence OpenFOAM en production n'est revendiquée.**

## 3. Qualification à trois maillages

La grille reconstruit la géométrie six blocs d'OpenFOAM mais n'est pas le polyMesh natif. Le cas reconstruit complet compte 20 540 cellules. Une tentative V11 à 1 itération, arrêt externe après 20 s, s'est arrêtée pendant `PRESSURE_LINEAR_BEGIN` : état **INCOMPLETE**, sans observables valides. Les limites des ressources locales empêchent la campagne N128/N256 complète ici.

Pour réaliser l'étude requise : importer les vrais `constant/polyMesh` natifs avec les mêmes BC SST/OpenFOAM ; vérifier les définitions et dimensions NASA du `U_ref` (défini près de x/h=-4, pas nécessairement la vitesse inlet), puis exécuter au moins trois maillages de même famille. La qualification exige vrais résidus des quatre champs, conservation, fonctions de paroi, Cf réattaché, profils x/h={1,4,6,10}, contraintes de Reynolds, et GCI sur une famille de grilles comparables. L'outil `qualify_bfs_driver_seegmiller_v10.py` reste fail-closed et ne calcule rien si le certificat est absent. **Nouveau code V11 : préparation A/B avec preuve SHA-256 du champ physique de chaque HDF5 et du choix numérique.**

## 4. Tests locaux

- Construction de `cfdx_production_solver` réussie, GCC14.2, HDF5, MPI OFF, sans GPU. Le gros fichier `main.cpp` est compilé `-O0` via `CFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON`, option explicitement DIAGNOSTIQUE ; ce n'est ni une campagne Release optimisée de performance, ni une qualification MPI/GPU.
- Huit cibles CTest ciblées passent (pression Rhie-Chow, guard Schur, PCD/équilibration/sélection, SST, mur).
- 81 tests Python BFS/PCD passent et 6 sont ignorés sur le périmètre élargi. Ces tests ne sont pas une CI complète.
- Six profils BFS et cinq profils PCD spécifiques sont effectivement tentés sur la même source HDF5, chaque trace conserve statut d'exécution, schéma, SHA-256, preuves et première étape bloquante. Échec et timeout doivent rester non qualifiés.
- Le garde de corruption de vectorisation Schur est exercé en C++ sur `CONVERGED`, `MAX_ITER_REACHED`, taille incorrecte, non-fini et statut DIVERGED.

## 5. Reproduction

```
cmake -S . -B build/bfs-v11 -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_ENABLE_MPI=OFF -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build/bfs-v11 --target cfdx_production_solver test_coupled_schur_action_guard test_pressure_rhie_chow_compensated -j2
ctest --test-dir build/bfs-v11 -R 'test_coupled_schur_action_guard|test_pressure_rhie_chow_compensated' --output-on-failure
python -m pytest -q tests/python/test_bfs* tests/python/test_pcd_fp_replay_campaign.py
python scripts/qualify_bfs_v11_numerical_matrix.py \
  --reference-case /chemin/vers/bfs_reference.cfdx.h5 \
  --binary build/bfs-v11/cfdx_production_solver \
  --out /chemin/vers/preuves_v11 --sst-inner-budget 800 --outer-iterations 1 --wall-timeout 60
```

Le code et les preuves exacts sont dans les fichiers livrés. Aucun push/merge GitHub n'a été effectué.
