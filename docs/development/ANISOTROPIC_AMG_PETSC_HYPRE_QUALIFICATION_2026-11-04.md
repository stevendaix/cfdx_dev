# CFDX — qualification des solveurs anisotropes (AMG natif, agrégation lissée, Krylov)

**Date :** 4 novembre 2026
**Base :** archive utilisateur `cfdx_dev-master(1).zip` (le SHA `0529a79` de la branche locale est un commit d'import, **pas** un SHA GitHub amont).
**Statut :** correctifs natifs codés et campagnes ciblées exécutées ; références PETSc/HYPRE facultatives **non exécutées** faute de bibliothèques ; PR GitHub **non publiée** faute d'accès réseau/écriture.

## 1. Objectif et frontière de vérification

Améliorer le coût total `préparation + résolution` sur des matrices elliptiques de diffusion/pression **anisotropes** sans changer la tolérance d'acceptation du **résidu physique vrai**. Comparer CG, GMRES, FGMRES, ILU(0), AMG de type Boomer natif, AMG à agrégation lissée ; prévoir une référence indépendante avec les vraies bibliothèques HYPRE/PETSc lorsque celles-ci sont disponibles.

**Ne pas confondre les produits :** `NativeBoomerAMGPreconditioner` est une implémentation CFDX inspirée des principes de BoomerAMG, **pas** `HYPRE_BoomerAMG`. `NativeSmoothedAggregationAMGPreconditioner` est une implémentation CFDX d'agrégation lissée, **pas** `PETSc.PCGAMG`. Le projet est autonome ; aucune API HYPRE/PETSc in-process n'est déclarée ni présentée comme intégrée par cette PR.

## 2. Défauts avérés et correctifs livrés

### A — Chemin de résolution cubique sur le dernier niveau AMG

**Défaut confirmé.** Une interpolation peut arrêter le grossissement alors que la dernière matrice compte des milliers d'inconnues. `smooth_coarsest()` utilisait malgré tout une factorisation dense en `O(n³)` / `O(n²)` mémoire. Sur le cas **8 192 inconnues, isotrope, réaction=1**, les niveaux retenus étaient notamment `8192 → 2692` pour l'agrégation lissée. Le benchmark était alors interrompu après 18 s, même avec une limite de 35 itérations Krylov. Le temps de *setup* était rapide : le défaut était dans l'**application** du V-cycle.

**Correction :** uniquement pour la matrice terminale de taille `n > 16`, employer une correction creuse linéaire **bornée à 24 balayages Jacobi amortis**, initialisée à zéro et vérifiée en valeurs finies. L'élimination dense pivotée reste réservée aux systèmes terminaux d'au plus 16 inconnues. La tolérance de qualification et le test de résidu final n'ont pas été modifiés. La correction grossière reste volontairement approximative : le solveur Krylov extérieur doit toujours vérifier l'équation initiale.

**Résultat :** la campagne à 8 192 inconnues termine désormais ; l'agrégation lissée isotrope affiche environ **0,017–0,019 s** pour 3 itérations FGMRES dans cette session, contre un dépassement de **18 s** avant correction. C'est une comparaison locale à compilation/plateforme similaires, **pas** une estimation de gain universel.

**Régression anti-réapparition :** `test_amg_sparse_coarse_regression`, taille 128 × 64, teste le dernier niveau à plus de 16 inconnues, le solveur FGMRES+AMG et la solution manufacturée ; timeout CTest de 10 s.

### B — Coût de sélection des agrégats et points C/F

**Défaut :** les sélections de graines parcouraient de façon répétée tous les sommets, donnant un coût de sélection quadratique sur des graphes creux.

**Correction :** file de priorité maximale avec entrées périmées ignorées, et mise à jour incrémentale des influences/degrés ; maintien de la même politique déterministe de départage (plus petit indice pour influence égale). La complexité de **sélection** devient en ordre de grandeur `O((n + nnz) log(n + nnz))`, hors interpolation et produit de Galerkin. Le produit de Galerkin reste susceptible de grossir ; aucune complexité linéaire globale de l'AMG n'est annoncée.

### C — Dispatch du solveur à partir de la matrice réelle

**Défaut :** le routeur pouvait choisir un couple solveur/préconditionneur d'après le seul nombre d'inconnues alors que des indicateurs numériques étaient disponibles.

**Correction :** mesurer les caractéristiques CSR (symétrie, dominance diagonale, positivité diagonale, disparité locale de couplages négatifs). Pour un problème elliptique avec **couplage directionnel local ≥ 10³**, choix automatique **FGMRES + AMG natif** à partir de 24 inconnues, plutôt que de se fonder sur la taille seule. Le critère est un **proxy algébrique**, non une mesure géométrique exacte du rapport d'aspect d'une maille.

Le calcul rapide évite de créer le graphe topologique biparti `lignes + colonnes` lorsque seules les caractéristiques nécessaires au dispatch sont demandées ; les diagnostics complets conservent cette fonctionnalité.

### D — Contrats Krylov, mise à l'échelle et résidus

- **CG** : refus explicite de la mise à l'échelle ligne/colonne non symétrique ; seule la congruence symétrique préserve l'admissibilité de principe. Une matrice mesurée non symétrique ou à diagonale non positive est rejetée avec un statut `NOT_APPLICABLE` quand CG a été forcé ; le test de positivité des courbures durant CG reste nécessaire pour exclure des matrices indéfinies non détectées statiquement.
- **FGMRES** : utilisé comme chemin robuste lorsque le préconditionneur n'est pas garanti constant/SPD ; cela n'affirme pas que toutes les méthodes AMG soient inconditionnellement sûres avec CG.
- **Mise à l'échelle** : résolution interne demandée à `0,1 × tolerance` par prudence puis **certification sur `||Ax-b||₂/||b||₂` pour la matrice physique d'origine** ; aucun assouplissement de seuil n'est accepté. La solution et son résidu sont reportés dans les coordonnées physiques.
- **ILU(0)** : la combinaison CG+ILU(0) non symétrique est refusée ; GMRES/FGMRES+ILU(0) sont exercés.

## 3. Qualification exécutée

La matrice est une discrétisation volumes-finis à cinq points d'une diffusion réaction-diffusion 2D sur grille structurée, avec poids des voisins horizontaux/verticaux `a_x/a_y` contrastés, **rapport 1, 10³ ou 10⁶**, et terme de réaction diagonal `+1`. Le second membre est fabriqué `b = A x_exact` avec une solution trigonométrique non constante. Pour **toutes** les variantes qualifiées, le gate vérifie : statut convergé, `||Ax-b||₂/||b||₂ ≤ 10⁻⁹` et `||x-x_exact||₂/||x_exact||₂ ≤ 2 × 10⁻³`. Les itérations et le temps (`setup + solve`, pas uniquement le produit matrice-vecteur) sont enregistrés.

Les maillages testés : **32 × 24** (CTest), **64 × 64**, **128 × 64** et **128 × 128** (jusqu'à **16 384 inconnues**). La direction de forte anisotropie a été inversée (x et y). Le cas de 8 192 inconnues reproduit spécifiquement le défaut d'agrégation qui déclenchait auparavant le solveur dense.

### Résultats sélectionnés, exécution Release (environnement de cette session)

| Problème | Combinaison | Itérations | Résidu vrai relatif | Temps setup + solve |
|---|---|---:|---:|---:|
| 8 192 DOF, isotrope | FGMRES + AMG natif | 3 | 1,65 × 10⁻¹¹ | 0,032 s |
| 8 192 DOF, isotrope | FGMRES + agrégation lissée | 3 | 3,14 × 10⁻¹⁰ | 0,017 s |
| 8 192 DOF, anisotropie x=10⁶ | Auto → FGMRES + AMG natif | 5 | 4,18 × 10⁻¹¹ | 0,027 s |
| 8 192 DOF, anisotropie x=10⁶ | FGMRES + agrégation lissée | 6 | 2,89 × 10⁻¹¹ | 0,086 s |
| 16 384 DOF, anisotropie x=10³ | FGMRES + AMG natif | 5 | 2,68 × 10⁻¹¹ | 0,061 s |
| 16 384 DOF, anisotropie x=10³ | FGMRES + agrégation lissée | 5 | 9,96 × 10⁻¹⁰ | 0,559 s |
| 16 384 DOF, anisotropie x=10⁶ | FGMRES + AMG natif | 5 | 4,18 × 10⁻¹¹ | 0,053 s |
| 16 384 DOF, anisotropie x=10³ | GMRES + ILU0 | 18 | 6,48 × 10⁻¹⁰ | 0,050 s |
| 16 384 DOF, anisotropie y=10³ | CG + agrégation lissée | 6 | 2,56 × 10⁻¹¹ | 0,479 s |

**Interprétation :** le meilleur choix n'est pas identique partout. L'AMG natif est nettement plus rapide que l'agrégation lissée pour plusieurs cas fortement directionnels testés, mais il ne faut pas généraliser ce classement à d'autres opérateurs et architectures. L'agrégation lissée corrigée conserve son intérêt comme second algorithme indépendant de coarsening.

**CTest Release :** `test_matrix_diagnostics`, `test_amg_dispatch_certificate`, `test_null_space`, `test_native_external_methods`, `test_matrix_characteristics`, `test_amg_preconditioner_qualification`, `test_n11_scaling_equivalence`, `test_n8_amg_anisotropy_scaling`, `test_anisotropic_krylov_matrix`, et `test_amg_sparse_coarse_regression` : **10/10** exécutés et réussis dans cette campagne.

**Référence externe calculée :** SciPy `spsolve` sur les six matrices 64 × 64 (deux axes, trois ratios) donne **6/6 validations** avec résidu vrai de l'ordre de `10⁻¹⁵` à `10⁻¹⁴` et erreur relative sur la solution ≤ `2,2 × 10⁻¹¹`. Elle confirme les solutions numériques attendues et sert de comparaison indépendante **des solveurs**, non de preuve d'identité algébrique des méthodes AMG.

## 4. PETSc / HYPRE : infrastructure de comparaison, non validation imaginaire

L'outil `tools/benchmarks/compare_petsc_hypre_anisotropy.py` assemble **les mêmes coefficients et la même solution manufacturée** que le C++ et mesure :

- PETSc `KSPFGMRES` + vrai `PCHYPRE`/`boomeramg` ;
- PETSc `KSPFGMRES` + `PCGAMG` ;
- PETSc `KSPGMRES` ou `KSPFGMRES` + `PCILU` ;
- PETSc `KSPCG` + `PCGAMG`, sous réserve de convergence certifiée ;
- indépendant sans PETSc : SciPy solveur sparse direct.

Il produit un CSV avec nombre d'itérations, **résidu original** et **erreur sur l'exacte**, durée de préparation PETSc et durée de résolution. En cas d'absence de `petsc4py` ou de PETSc, le mode PETSc quitte avec **code 77 = SKIP** ; une combinaison non disponible / non convergée n'est jamais inscrite comme validée.

**État de cet environnement :** `petsc4py`, PETSc, HYPRE et MPI ne sont pas installés ; la tentative PETSc/HYPRE a réellement retourné `SKIP(77)`. Les performances HYPRE/PETSc ne sont donc **pas** renseignées ; il est impossible d'affirmer ici que CFDX bat ces logiciels.

Pour installer un banc de référence sur une machine équipée de PETSc construit avec HYPRE :

```bash
python3 tools/benchmarks/compare_petsc_hypre_anisotropy.py --backend petsc --nx 64 --ny 64 > petsc_hypre_4096.csv
```

## 5. Reproduction des tests natifs

```bash
cmake -S . -B build-aniso -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_USE_HDF5_SERIAL=OFF \
  -DCFDX_BUILD_RUNTIME_TESTS=OFF -DCFDX_BUILD_TESTS=ON
cmake --build build-aniso --target test_anisotropic_krylov_matrix \
  test_amg_preconditioner_qualification test_n8_amg_anisotropy_scaling \
  test_amg_dispatch_certificate -j 4
ctest --test-dir build-aniso --output-on-failure \
  -R 'test_anisotropic_krylov_matrix|test_amg_sparse_coarse_regression|test_amg_preconditioner_qualification|test_n8_amg_anisotropy_scaling|test_amg_dispatch_certificate'
./build-aniso/test_anisotropic_krylov_matrix --benchmark 128 128 > native_16384.log
python3 tools/benchmarks/compare_petsc_hypre_anisotropy.py --backend scipy --nx 64 --ny 64 > scipy_4096.csv
```

`--benchmark` fournit des chronos indicatifs, **sans gate de temps rigide** afin d'éviter les faux échecs CI. Le test dédié du solveur grossier dispose, lui, d'un timeout assez large pour détecter un retour au comportement catastrophique (solveur dense de plusieurs milliers d'inconnues).

## 6. Limites et travaux non clos

1. **PAS de qualification HYPRE/PETSc runtime** dans cet environnement ; le script de comparaison doit tourner sur une installation de référence vérifiée.
2. **Pas de MPI/GPU ni test CFD NACA N256 couplé** : les tests actuels sont des opérateurs elliptiques synthétiques représentant les difficultés de pression, pas une validation de conservation Navier–Stokes en 2D/3D et pas une preuve d'accélération bout-en-bout.
3. **Pas de garantie d'erreur avant pour les matrices quasi singulières.** Une expérimentation supplémentaire pour anisotropie `10⁹` + réaction `10⁻⁶` a montré une erreur relative en solution de l'ordre de `10⁻³` à `10⁻²` malgré des résidus vrais faibles. Dans ce régime le conditionnement se rapproche de la limite des doubles ; maintenir une gate *séparée* de précision de solution et de résidu, et envisager un changement d'échelle, un traitement du noyau et un raffinement itératif en précision élargie.
4. **Détecteur anisotrope local** ne remplace pas l'analyse géométrique (non-orthogonalité, transitions taille, coordonnées faces) nécessaire sur les vrais maillages CFDX. Les systèmes nonsymétriques et de point-selle couplés exigent des matrices représentatives dédiées.
5. **Extensions recommandées** avant de prétendre égaler HYPRE/PETSc : interpolation tronquée à complexité contrôlée, coarsening PMIS/HMIS et semi-coarsening pour fortes anisotropies ; contrôle de croissance Galerkin, lissage symétrique/ligne ou bloc, hiérarchie persistante entre itérations non linéaires ; opérateur near-null explicite, coarse solvers itératifs avec budget/contrat; moniteurs de résidu réel et erreur avant sur MMS.
6. **Sanitizers :** build DebugSanitizers ASan/UBSan terminé ; `test_amg_preconditioner_qualification`, `test_anisotropic_krylov_matrix` et `test_amg_sparse_coarse_regression` ont passé **3/3** en campagne dédiée. Cette couverture ciblée ne constitue pas une validation sanitizer de tous les modules CFDX.
7. **Publication PR** : la copie de cette archive n'a pas de remote GitHub et l'environnement ne résout pas `github.com`. Le livrable prêt à réappliquer est un commit/patch local, **pas une PR distante ouverte**.

## 7. Références méthodologiques (algorithmes, pas portage du code)

- HYPRE BoomerAMG (coarsening HMIS, interpolation et aggressive coarsening) : https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html
- PETSc PCGAMG (agrégation, complexité, interpolation, quasi-noyau) : https://petsc.org/main/manualpages/PC/PCGAMG/
- PETSc PCGAMGAGG (filtrage du prolongateur conservant le noyau) : https://petsc.gitlab.io/petsc/main/manualpages/PC/PCGAMGAGG/
- petsc4py `KSP.setGMRESRestart()` : https://petsc.org/main/petsc4py/reference/petsc4py.PETSc.KSP.html

## 8. Déploiement proposé de la PR

Titre : `perf(linalg): bound AMG coarse solves, accelerate coarsening, qualify anisotropic Krylov paths`.

Critères de merge : replay du patch sur HEAD exact; cibles CTest linéaires vertes Release et DebugSanitizers; campagne anisotropie x/y sur 4k/8k/16k ; comparaison PETSc/HYPRE sur hôte qui en dispose; tests de pression N128/N256 en circuit CFD complet; absence de régression en temps complet. Le lot présent constitue un premier correctif natif vérifié, **pas** la validation finale de la chaîne pression–vitesse CFDX.