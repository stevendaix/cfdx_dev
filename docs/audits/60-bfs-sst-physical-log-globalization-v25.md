# CFDX V25 — globalisation logarithmique SST et dérivées réassemblées aux parois

**Date :** 10 novembre 2026. **Statut : expérimental, désactivé par défaut. Qualification physique : NON acquise.**

## 1. Cause examinée

Sur le BFS six blocs réduit (712 cellules), la V24 obtient à l'itération interne SST 10 une direction spatiale Newton–Krylov dont le vrai résidu directionnel GMRES vaut `0.16613194317`, inférieur au seuil **linéaire** expérimental `0.25`, mais le pas ne passe pas le mérite utilisé par le solveur. Le meilleur candidat V22 Newton–Jacobi est alors accepté avec mérite SST réassemblé `0.88990745136`, l'état extérieur restant rejeté car les équations SST ne sont pas fermées.

Deux évaluations du mérite coexistent intentionnellement : V22 agit sur un système en inconnues logarithmiques, mis à l'échelle avec des coefficients **figés**, alors que les portes du solveur utilisent le résidu FV **réassemblé** et les normalisations de lignes recomputées. Le pas testé dans V22 ne peut donc pas être déduit du seul résidu directionnel GMRES. Cela n'autorise pas à remplacer le critère physique par le critère Krylov.

## 2. Modifications de production, toutes opt-in

- `sst_log_globalization_v25_experimental=true` exige **aussi** `sst_jfnk_v22_experimental=true` et `sst_spatial_block_preconditioner_v24_experimental=true` dans `numerics.raw_settings`; tous restent désactivés par défaut.
- Le noyau partagé `coupled_newton_krylov_step` expose maintenant la direction *après satisfaction du vrai critère GMRES*, même si sa recherche de pas avec mérite gelé est rejetée. Cette exposition est diagnostique et **ne change aucune décision d'acceptation des appels existants**.
- `sst_globalize_log_direction_v25()` parcourt des mises à jour en `exp(alpha*direction_log)` : rayon de confiance maximal 1,5 en log, puis au plus 12 subdivisions par 2. Chaque candidat est fini, positif, re-projette la contrainte propriétaire lorsque demandée, recalcule `G`, les faces et le résidu SST FV avec normalisations courantes. Il n'est admissible que s'il diminue **strictement** ce vrai mérite ; la tolérance SST `1e-8` n'est jamais modifiée.
- La proposition V25 entre en compétition avec Picard, V22 Newton–Jacobi, V23 blocs et V24 spatial. La sélection finale est faite à nouveau par le solveur avec sa porte physique habituelle. Échec => ancienne proposition, ou rollback si rien n'est admissible.
- `sst_correct_wall_blocks_v25()` remplace au plus **4** blocs locaux de paroi du préconditionneur par les dérivées numériques réassemblées de `[R_k,R_omega]` en log(k), log(omega). La sélection par défaut porte sur les propriétaires de paroi où le défaut `R_k` normalisé est le plus important. `CFDX_SST_V25_WALL_CELL=357` priorise une cellule propriétaire **explicitement désignée pour l'audit**, sans hypothèse sur sa numérotation dans le code général. La correction n'est appliquée que si son écart dépasse 1 %, que le facteur reste inversible et que toutes les évaluations restent admissibles. Les autres blocs restent inchangés.
- `SST_JFNK_V22` journalise `v25_direction`, `v25_attempted`, `v25_valid`, `v25_alpha`, `v25_first_merit`, `v25_best_trial_merit`, `v25_wall_sampled_ids`, `v25_wall_corrected`, `v25_wall_relative_error`, `selected`.

Ces mécanismes n'installent pas PETSc et n'affirment pas être `SNESLINESEARCHBT`. La recherche V25 est une politique bornée de diagnostic et de choix ; elle n'introduit **aucune modification** du modèle SST ni de ses conditions aux limites.

## 3. Résultats de production : même HDF5 réduit, 10 itérations SST

| Critère à l'itération interne 10, première tentative extérieure | V24 | V25, sélection automatique des parois | V25, audit paroi 357 |
|:--|--:|--:|--:|
| Vrai résidu GMRES, direction spatiale | 0.16613194317 | 0.16613194317 | 0.23209984454 |
| Pas V25 logarithmique physique trouvé | Sans objet | Oui | Non |
| Alpha logarithmique | Sans objet | 0.00087718879 | 0 |
| Mérite V25 candidat, avant sélection | Sans objet | 0.92595972137 | Aucun pas admissible |
| Mérite physique avant sélection | 0.92596107454 | 0.92596107454 | 0.92596107454 |
| Mérite **effectivement** sélectionné | **0.88990745136** | **0.88990745136** | **0.88990745136** |
| Méthode sélectionnée | Newton–Jacobi V22 | Newton–Jacobi V22 | Newton–Jacobi V22 |
| Blocs paroi corrigés | 0 | 0 | **1 (cellule 357)** |
| Écart jacobien propriétaire corrigé | Non applicable | Non appliqué | **11.542184024 %** |
| État extérieur convergé | Non | Non | Non |

La V25 démontre qu'il **existe** un petit pas descendant sur une direction spatiale auparavant rejetée par la recherche avec échelles gelées. Ce pas est moins efficace que Newton–Jacobi, et **n'est pas choisi**. La correction exacte du bloc 357 ne suffit pas ; elle modifie la direction et supprime même ce petit pas admissible, alors que le critère linéaire de 0,25 est toujours satisfait. Le solveur conserve le résultat de référence V22 dans les deux cas : pas de régression du mérite à l'itération étudiée, mais pas de solution globale obtenue non plus.

### Réserves strictes

- Le budget `sst_transport_max_iterations=10` est **diagnostique**, la tolérance `1e-8` reste inchangée. Le calcul échoue à la porte des équations retournées et **annule l'état externe**. Les retry extérieurs doivent être ignorés pour cette comparaison.
- La réévaluation du champ et de la cible de paroi en V25 n'est pas une démonstration d'équivalence numérique aux pondérations `updateCoeffs()` et `manipulateMatrix()` d'OpenFOAM Foundation.
- Le balayage spatial V24 reste CPU **séquentiel**. Aucune qualification multicœur CPU, GPU, MPI ou multi-GPU de cette voie.
- Les six Schur/PCD/AMG/MGR sur le grand BFS, la convergence à 800 itérations, les trois états extérieurs, NACA, NASA et les autres modèles sont **non testés** dans ce lot. Aucun seuil physique assoupli, aucun push GitHub.

## 4. Tests et reproduction

Tests ciblés : 8/8 CTest, incluant 7 tests C++ V25 et 4 tests Python de protocole dans le projet de test. Compilation du solveur de production avec **Clang -O0 diagnostic** et liaison avec la bibliothèque native CFDX : réussies (aucun benchmark Release optimisé).

```sh
cmake -S . -B build-v25 -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v25 --target cfdx_production_solver test_sst_globalization_v25
ctest --test-dir build-v25 -R 'test_sst_globalization_v25|test_bfs_sst_v25_protocol' --output-on-failure
# Copier le vrai cas V24 712-cellules, ajouter l'option ci-dessus dans
# numerics.raw_settings, conserver tolérance et BC HDF5 strictement identiques.
CFDX_SST_DIAGNOSTICS=1 build-v25/cfdx_production_solver \
    --mesh bfs_v25_on_10.cfdx.h5 --output-dir output_v25 \
    --iterations 1 --experimental-sst
# Diagnostic exact pour un owner :
CFDX_SST_DIAGNOSTICS=1 CFDX_SST_V25_WALL_CELL=357 \
    build-v25/cfdx_production_solver --mesh bfs_v25_on_10.cfdx.h5 \
    --output-dir output_owner357 --iterations 1 --experimental-sst
python3 scripts/diagnose_bfs_sst_v25.py v25_first_attempt.log
```

Les journaux et HDF5 utilisés sont dans les preuves V25 ; les exécutions ont été interrompues volontairement après avoir capturé le premier rejet extérieur : **ce ne sont pas des calculs longs terminés**. L'arrêt contrôlé conserve le statut non qualifié.

## 5. Étape technique suivante

1. Vérifier indépendamment le produit directionnel `Jv` à deux amplitudes de perturbation (effet des contraintes propriétaires et précision soustractive), et la qualité de la direction contre la norme `L∞` du résidu dynamique. V25 montre que la direction spatiale peut être descent en pratique mais trop faiblement.
2. Améliorer la stratégie de construction du préconditionneur paroi : choisir les cellules où les dérivées de `G`, `nu_t`, diffusion et cibles de paroi sont importantes **sans index de maillage codé en dur**, puis tester la régularité de ces blocs ; ne pas promouvoir un bloc vrai mal conditionné.
3. Tester une recherche de pas fondée sur un objectif lisse supplémentaire, par exemple une norme L2 équilibrée, **uniquement comme stratégie pour proposer des candidats**, avec acceptation finale sur le même mérite physique L∞ et les mêmes portes SST.
4. Obtenir au plus 800 itérations SST sur 712 cellules avec plusieurs états extérieurs effectivement acceptés, puis élargir à 20 540 cellules et aux six variantes COUPLED sous un protocole identique ; aucune validation NASA avant convergence et comparabilité expérimentale.

### Bibliographie de méthode

- PETSc SNES : https://petsc.org/release/manual/snes/ — plusieurs recherches de pas et séparation Jacobien / préconditionnement / résidu non linéaire.
- PETSc SNESLINESEARCHBT : https://petsc.org/main/manualpages/SNES/SNESLINESEARCHBT/ — contexte théorique ; **non implémenté à l'identique**.
- OpenFOAM Foundation v11 : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html — `updateCoeffs()`, poids par face et manipulation matricielle, **équivalence non démontrée**.
