# CFDX — BFS V15 : action PCD de rang plein, preuve du résidu extérieur et triplet de paroi atomique

**État :** développement expérimental ; aucune qualification aérodynamique ou RANS Driver–Seegmiller.

## 1. Résumé et provenance

Base : livraison CFDX BFS V14 du 9 novembre 2026, reconstruite depuis l'archive utilisateur `cfdx_dev-master (2)(1).zip` via le patch cumulatif. V15 n'est pas un HEAD distant ni une PR intégrée.

Les tests sont des contrôles logiciels et numériques ciblés. Ils ne prouvent ni la convergence de la campagne longue des six Schur, ni le raccordement des trois conditions de paroi OpenFOAM à l'opérateur volumes finis.

## 2. PCD : noyaux gauche/droit et rang de l'action

Pour une matrice auxiliaire potentiellement singulière `Fp`, les bases indépendamment vérifiées satisfont `Fp Z ≈ 0`, `Fpᵀ W ≈ 0`, `ZᵀZ ≈ I`, `WᵀW ≈ I`. Le résidu relatif original inclut le plancher de compatibilité `||Wᵀ b||/||b||` ; il n'est jamais remplacé par le résidu projeté.

L'action historique de préconditionnement supprimait `Wᵀb` puis `Zᵀx`. Cela faisait perdre des directions au préconditionneur PCD, même lorsque l'opérateur auxiliaire déflaté était inversible. Avec `pcd_fp_full_rank_action_v15_experimental=true`, V15 résout `[(Fp + WZᵀ)] x = b` **avec le second membre complet et sans projection finale de x**. Il calcule séparément le résidu de cet opérateur déflaté et celui de `Fp` original. `Fp` et la matrice COUPLED physique ne sont pas modifiés.

Le test C++ sur une matrice singulière 2×2 non symétrique reproduit le défaut de rang, puis valide un résidu vrai déflaté de `2,48253e-16` et la conservation de la direction du noyau droit.

## 3. Actualisation lorsque Fp change

Les bases sont identifiées par l'empreinte FNV64 des **coefficients ET de la structure** CSR de Fp. L'option `pcd_fp_bisided_basis_directory_v15_experimental` choisit une base pour la matrice numérique courante. Changement de coefficients → nouvelle empreinte, nouveau fichier `.cfpb` exigé. Les bases sont recertifiées en C++ sur la matrice courante. Base absente, périmée ou non certifiable → échec explicite.

`python scripts/build_pcd_fp_basis_registry_v15.py --capture X.cfp --capture Y.cfp --directory registry` prépare ce registre par SVD dense bornée. Le backend fourni refuse les matrices de plus de 1024 lignes : ce n'est pas une actualisation dense viable pour 20 540 cellules et ce n'est pas une SVD incrémentale de production. Un calcul robuste scalable des noyaux restera nécessaire.

## 4. Résidu physique COUPLED et contraction

Le résidu vrai de la matrice COUPLED `||b-Ax||/||b||` est recalculé en précision étendue même lorsque FGMRES échoue, sans changer ni les seuils ni le statut Krylov. Un statut `CONVERGED` dont le résidu recalculé dépasse la tolérance est refusé.

Sur le BFS 712 cellules, les deux modes projeté et rang plein ferment les applications Fp internes à ~`1e-8`, mais **le système COUPLED extérieur ne converge pas**. Le diagnostic du mode rang plein aboutit à un résidu COUPLED réel `0,0369139` avec statut d'échec après 1124 itérations internes, et aucun état physique qualifié. Aucun troisième état accepté n'a été obtenu.

## 5. Triplet SST de paroi préparé atomiquement

Le nouveau module `sst_openfoam_wall_triplet_v15.h` calcule `nu_t`, `G` et `omega` sur la **même** vitesse et le **même** état SST. Il applique une moyenne de face pondérée par aire pour la production des cellules touchant plusieurs faces. Il vérifie la qualité géométrique, les valeurs physiques et la complétude des patches. Le résultat n'est publié qu'après `transaction.commit()` de l'itération externe.

BFS réduit : `70` faces, `69` cellules propriétaires, `1` cellule multi-faces. Cas SST rejeté à 10 itérations internes : aucun triplet publié. Dans les traces, `momentum_applied=0`, `sst_applied=0`, `openfoam_equivalent=0` restent explicites. **La condition aux limites OpenFOAM est encore à raccorder réellement aux trois assemblages de production.** Aucun défaut de cette implémentation de référence ne doit être confondu avec une preuve d'équivalence.

## 6. Qualification et portes restant fermées

- FGMRES extérieur des six Schur : non convergé / non qualifié ; les tests de l'inverse Fp ne suffisent pas.
- Trois états externes acceptés et contraction stricte du résidu global `U,p,k,omega,nu_t` : non acquis.
- Fonctions de paroi complètes et indépendamment validées : non acquises (triplet V15 lecture seule).
- Maillages OpenFOAM natifs et comparaison Driver–Seegmiller (`x_r/h`, `Cf`, vitesse, contraintes de Reynolds) : non exécutés.
- Tests MPI/GPU, DebugSanitizers et CI complète : non exécutés.
- Aucun seuil physique n'a été abaissé ; aucun push GitHub.

## 7. Reproduction minimale

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_GPU=OFF -DCFDX_ENABLE_MPI=OFF -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build --target cfdx_production_solver test_pcd_full_rank_action_v15 test_sst_wall_triplet_v15 -j 2
ctest --test-dir build -R 'test_(pcd_full_rank_action_v15|sst_wall_triplet_v15)' --output-on-failure
python -m pytest -q tests/python/test_build_pcd_fp_basis_registry_v15.py
```

Pour BFS, activer `v15_stage_openfoam_wall_triplet_reference=true` et `pcd_fp_full_rank_action_v15_experimental=true` dans les paramètres bruts du cas ; conserver les tolérances d'origine. La valeur `pcd_fp_bisided_basis_directory_v15_experimental` ne doit viser qu'un répertoire de certificats dérivés des captures exactes.

## 8. Sources méthodologiques

- PETSc KSP : https://petsc.org/main/manual/ksp/ (systèmes singuliers, modes gauches/droits, résidus préconditionnés et vrais, limites de FGMRES à droite).
- OpenFOAM : `omegaWallFunction`, `nutUSpaldingWallFunction`, `kqRWallFunction` ; vérifier leur couplage dans l'implémentation native et les conditions de Driver–Seegmiller avant de prétendre à l'équivalence.
