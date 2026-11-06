# CFDX — adaptation native ciblée des méthodes PETSc/HYPRE pour l'anisotropie

**Date :** 2026-11-05
**Source de départ :** `CFDX_AMG_anisotrope_archive_corrigee_2026-11-04.zip` fourni dans la conversation.
**Objectif :** exploiter des recettes algorithmiques utiles de PETSc et HYPRE **sans installer ni lier** ces bibliothèques ; préserver l'interface CFDX, les critères de résidu vrai et les tests existants.

## Portée et attribution

- PETSc et HYPRE **ne sont pas intégrés comme bibliothèques**. Aucun nom de méthode CFDX ne doit être interprété comme une implémentation exacte de `HYPRE_BoomerAMG` ou de `PETSc.PCGAMG`.
- Les nouveaux algorithmes sont des **réimplémentations C++17 de CFDX**, pas une copie textuelle de leurs fichiers sources. Les fonctions PETSc/HYPRE ne sont pas appelées. Si du code source upstream est ultérieurement recopié directement, il faudra importer l'en-tête de copyright correspondant et inclure les notices de licence.
- Licence HYPRE : **MIT OU Apache-2.0**, https://github.com/hypre-space/hypre . Licence PETSc : **BSD-2-Clause**, https://petsc.org/main/install/license/ . Dans l'état livré, aucune copie verbatim de code upstream n'exige de redistribuer un fichier `third_party`.
- Références de conception : HYPRE BoomerAMG (graphe de connexions fortes, HMIS/PMIS, lissage et interpolation) https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html ; PETSc PCGAMGAGG (filtrage du prolongateur) https://petsc.org/main/manualpages/PC/PCGAMGAGG/ ; PETSc KSPFGMRES https://petsc.org/main/manualpages/KSP/KSPFGMRES/ .

## Implémenté

### 1. Lignes algébriques anisotropes (CFDX indépendant)

**Nouveau :** `src/cfdx/core/linalg/algebraic_line_preconditioner.h`.

1. Parcours du CSR et identification des arêtes à coefficient négatif dominant dans chaque ligne ; exigence de paire transposée compatible pour éviter de construire de fausses directions sur des matrices convection-diffusion non symétriques.
2. Tri déterministe, filtre de direction (rapport 32 entre forts et faibles coefficients), sélection des arêtes par union-find, sans cycle et avec degré maximal de deux : chaque bloc est un chemin tridiagonal.
3. Factorisation tridiagonale en préparation, application en O(n) ; pivots strictement positifs, tests de valeurs finies, rejet des mises à jour numériques invalides ou qui changent la direction forte sans modifier la factorisation précédente.
4. Sur un opérateur isotrope, les lignes détectées ont longueur un : comportement block-Jacobi réduit au Jacobi scalaire.

**Limites :** ce détecteur est algébrique, non géométrique ; la qualité dépend du signe des coefficients et de la forte hiérarchie des connexions. Les opérateurs non-M-matrix ou non symétriques ne sont pas qualifiés pour CG avec cette option. Il n'a pas vocation à remplacer les lisseurs des problèmes généraux.

### 2. V-cycle AMG avec lissage directionnel (référence méthodologique HYPRE)

`src/cfdx/core/linalg/amg_preconditioner.h` : politique opt-in `AMGSmootherPolicy::AlgebraicLines` sur chaque niveau, pré/post-correction par blocs tridiagonaux si la couverture est suffisante ; maintien du Jacobi d'origine ailleurs. Si le bloc grossier ne satisfait pas les pivots, le niveau continue au Jacobi ; le nombre de ces replis est exposé via `line_factor_fallback_levels()` et n'est pas interprété comme une validation d'un solveur direct.

Le choix FGMRES évite d'affirmer à tort la SPD/symétrie d'un V-cycle qui comporte des lisseurs directionnels et des repliements. Un CG explicitement forcé avec les nouvelles variantes AMG est rejeté par le dispatch.

### 3. Filtre de prolongateur à agrégation lissée (référence méthodologique PETSc)

Ajout d'un seuil relatif aux poids de chaque ligne de la prolongation lissée (par défaut zéro, **0,05 dans la variante opt-in anisotrope**). Après filtrage, chaque ligne est renormalisée : sur le mode constant, \(P\\mathbf 1 = \\mathbf 1\). Les entrées nulles ou non finies sont rejetées. Le Galerkin \(A_c=P^T A P\) reste celui de CFDX.

Sur la matrice 128×128, forte anisotropie x, \(10^3\), les entrées du premier prolongateur passent de **59 648 à 27 136**, avec 5→10 niveaux dans cette configuration. L'impact est mesuré sur le coût total, pas seulement le nombre d'itérations.

### 4. Trois préconditionneurs sélectionnables

Dans `LinearSolverRequest.preconditioner` :

- `PreconditionerModel::AnisotropicLine` — blocs lignes seuls, avec Krylov externe ;
- `PreconditionerModel::AnisotropicBoomerAMG` — C/F de CFDX et lissage par lignes ;
- `PreconditionerModel::AnisotropicSmoothedAggregationAMG` — agrégation lissée de CFDX + lignes + filtrage 5 %.

Les variantes existantes `NativeAMG` et `SmoothedAggregationAMG` restent **inchangées comme choix par défaut**, car les résultats de vitesse ne sont pas uniformes. Lorsqu'un AMG à lignes est explicitement choisi avec `krylov=Auto`, le dispatch impose FGMRES ; le CG forcé est refusé. Les nouveaux AMG exigent un opérateur elliptique symétrique à diagonale positive.

**Exemple :**

```cpp
LinearSolverRequest request;
request.krylov = KrylovModel::FGMRES;
request.preconditioner = PreconditionerModel::AnisotropicSmoothedAggregationAMG;
request.gmres_restart = 35;
Vector solution(A.n_rows(), 0.0);
auto result = solve_linear_system(A, rhs, solution,
                                  LinearProblemKind::PressurePoisson,
                                  request, 250, 1e-9);
// Qualification uniquement si status == CONVERGED ET
// result.physical_residual_relative <= 1e-9.
```

## Résultats mesurés — matrices discrétisées indépendantes des solveurs

Matrice SPD de type volumes-finis, diffusion à cinq points + terme de réaction 1, solution analytique manufacturée (même CSR, RHS, initialisation et seuil pour les anciens/nouveaux AMG). Directions x/y et ratios 1000/1 000 000 ; taille 128×128, **16 384 inconnues**. Benchmark couplé **3 répétitions** de chaque outil, médiane du temps préparation+résolution. Matériel/charge de cette session seulement : ne pas traiter comme une comparaison de performance avec les logiciels PETSc/HYPRE officiels.

| Direction | Ratio | AMG originel | Variante native optimisee | Iterations ancien → nouveau | Mediane ancien → nouveau (s) | Rapport ancien/nouveau |
|---|---:|---|---|---:|---:|---:|
| x | 10³ | C/F | C/F + lignes | 5 → 2 | 0.0484 → 0.0453 | **1.07×** |
| x | 10⁶ | C/F | C/F + lignes | 5 → 2 | 0.0515 → 0.0446 | **1.15×** |
| y | 10³ | C/F | C/F + lignes | 5 → 2 | 0.0513 → 0.0580 | **0.88× (ralentit)** |
| y | 10⁶ | C/F | C/F + lignes | 5 → 2 | 0.0479 → 0.0600 | **0.80× (ralentit)** |
| x | 10³ | Agrégation lissée | Lignes + filtre | 5 → 3 | 0.4573 → 0.0434 | **10.54×** |
| x | 10⁶ | Agrégation lissée | Lignes + filtre | 6 → 2 | 0.1509 → 0.0346 | **4.36×** |
| y | 10³ | Agrégation lissée | Lignes + filtre | 6 → 3 | 0.4179 → 0.0474 | **8.81×** |
| y | 10⁶ | Agrégation lissée | Lignes + filtre | 6 → 2 | 0.1467 → 0.0368 | **3.99×** |

**Ablation contrôlée de l'agrégation lissée, direction x / ratio 1000 / 16 384 inconnues** (une exécution instrumentée avec 4 configurations) :

| Filtrage | Lignes | Itérations | Setup (s) | Solve (s) | Nb. NNZ du premier prolongateur |
|---|---|---:|---:|---:|---:|
| Non | Non | 5 | 0.421 | 0.417 | 59 648 |
| 5 % | Non | 6 | 0.022 | 0.042 | 27 136 |
| Non | Oui | 2 | 0.484 | 0.479 | 59 648 |
| 5 % | Oui | 3 | 0.025 | 0.050 | 27 136 |

**Conclusion de l'ablation :** le filtrage de la prolongation produit le principal gain de vitesse mesuré ; le lissage par lignes réduit les itérations mais n'est pas toujours plus rapide. L'analyse en setup/solve confirme qu'une diminution du nombre d'itérations ne suffit pas à valider une optimisation.

**Réserves numériques :** pour le préconditionneur par lignes seul, le défaut sur la solution manufacturée peut être de l'ordre de 10⁻⁶ lorsque le ratio atteint 10⁶ malgré un résidu relatif inférieur à 10⁻⁹ ; le conditionnement explique cette différence. Les variantes AMG donnent généralement des erreurs de solution beaucoup plus petites sur ces cas. Il faut conserver **deux critères indépendants** : résidu vrai et erreur de solution lorsque celle-ci est connue.

## Tests et preuves

- `test_petsc_hypre_native_port_anisotropy` : huit groupes incluant axes x/y, ratio 10³/10⁶, anisotropie diagonale à 45°, invariance du mode constant de P, rejet d'un couplage CG non qualifié, invalidation sans mutation et mise à jour numérique du préconditionneur.
- Campagne Release ciblée (2026-11-05) : **10/10 tests CTest** passant, incluant les gardes AMG antérieurs, la qualification anisotrope et la prévention de résolution dense d'un niveau grossier trop grand.
- Le test ciblé a également été exécuté sous ASan+UBSan : **8/8 groupes réussis**, sur les chemins concernés (voir journal fourni).
- Campagnes détaillées : `benchmark_native_petsc_hypre_port.py`, CSV de répétitions, médianes et test ablation.

### Reproduire

```bash
cmake -S . -B build-native -DCMAKE_BUILD_TYPE=Release -DCFDX_ENABLE_MPI=OFF
cmake --build build-native --target test_petsc_hypre_native_port_anisotropy test_anisotropic_krylov_matrix -j 4
ctest --test-dir build-native --output-on-failure -R '(test_petsc_hypre_native_port_anisotropy|test_anisotropic_krylov_matrix)'
./build-native/test_petsc_hypre_native_port_anisotropy --benchmark 128 128
./build-native/test_petsc_hypre_native_port_anisotropy --ablation 128 128
python3 tools/benchmarks/benchmark_native_petsc_hypre_port.py \
  --base-binary build-native/test_anisotropic_krylov_matrix \
  --port-binary build-native/test_petsc_hypre_native_port_anisotropy \
  --nx 128 --ny 128 --runs 3 --output medians.csv --raw raw.csv
```

## Étapes avant promotion en production

1. Exécuter les matrices de pression réelles CFDX (N128/N256, cellules très étirées/non orthogonales) ; vérifier la conservation locale, la réduction du résidu vrai et la qualité finale des forces. Ces campagnes physiques **n'ont pas été exécutées dans ce lot**.
2. Tester les systèmes Neumann purs et les quasi-noyaux de pression, y compris les préconditionneurs à lignes dans les formulations qui nécessitent une projection explicite. Le benchmark actuel utilise un terme de réaction et n'exerce pas de nullspace singulier.
3. Tester l'interpolation filtrée sur des maillages non structurés, des graphes anisotropes non-M-matrix, des couplages signés et des modes lents autres que la constante ; requalifier la robustesse en cas de modification du graphe de force au cours des itérations non linéaires.
4. Qualifier les backends MPI/GPU et la thread-safety du contexte de diagnostic (`NativeAMGPreconditioner::last_error` n'est pas synchronisé) avant usage concurrent.
5. Si une comparaison à HYPRE BoomerAMG / PETSc PCGAMG **officiels** est demandée, exécuter séparément les versions officielles sur les mêmes CSR avec le même RHS, car **elles ne sont pas installées ni utilisées** ici.
6. Ne pas présenter ces développements comme un port complet de HYPRE/PETSc ni comme une PR distante fusionnée ; un patch Git et une archive sont fournis pour révision et publication.

## Provenance exacte du code

Les nouveaux fichiers et modifications sont du code C++17 développé pour CFDX en utilisant les publications et documentations PETSc/HYPRE comme références. **Aucune portion textuelle amont de PETSc ou HYPRE n’est vendorizée dans ce lot.** Un portage de fichiers amont originaux, même partiel, nécessiterait une analyse de leurs contrats et dépendances puis l’attribution de chaque fichier effectivement repris. La solution livrée est autonome, sélective et testée, mais ne constitue ni PETSc complet ni BoomerAMG officiel.