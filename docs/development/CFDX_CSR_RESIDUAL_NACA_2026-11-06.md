# CFDX — matrices CSR, certification du résidu physique et convergence NACA

**Date :** 2026-11-06
**Base locale exacte :** `cfdx_1182_sources.zip` (1 509 fichiers avant ce lot).
**Statut :** correctifs algébriques qualifiés sur tests ciblés ; qualification CFD complète **non acquise**.
**CI GitHub :** non exécutée ni certifiée sur le HEAD de la PR #1182.

## Problèmes reproduits

### 1. Résidu relatif sur un système homogène

Dans `solve_reusable_scaled_cg`, le résidu relatif physique utilise auparavant `||b-Ax||/max(||b||,1e-300)`. Pour `b=0` et une approximation de haute précision, cette formule produit une erreur relative artificiellement gigantesque, tout en rapportant un résidu absolu nul à tort (`relative * ||b||`).

**Correction :** calculer le résidu initial vrai `r0=b-A*x0` ; choisir le dénominateur `||b||_2` si non nul, sinon `||r0||_2`. Après résolution, recalculer le résidu **sur la matrice et le vecteur originaux** : `r=b-A*x`, `residual=||r||_2`, `residual_relative=||r||_2/reference`. Cette règle ne réduit aucune tolérance.

**Contre-exemple SPD 3x3, `b=0` :** avant, relatif déclaré `1.69e284`, alors que la solution retournée vérifie `||Ax||=1.69e-16`; après, relatif vrai `3.14e-17` et statut `CONVERGED`. Le test dédié est dans `test_reusable_scaled_pressure_cg_v9.cpp`.

### 2. Doublons dans une matrice COO/CSR

`SparseMatrix::finalize()` triait auparavant les contributions mais ne fusionnait pas les entrées de même `(ligne,colonne)`. `A*x` cumulait tous les doublons ; des préconditionneurs comme Jacobi ne prenaient que le **premier** coefficient diagonal. Il en résultait des opérateurs algébriques incohérents et une dépendance à l'ordre d'assemblage.

**Correction :** CSR canonique avec colonnes strictement croissantes et une seule entrée par `(i,j)` ; somme en précision intermédiaire étendue ; refus des sommes hors plage `double` ; construction transactionnelle avant publication. Les valeurs nulles explicites sont conservées pour ne pas altérer silencieusement la structure. Les diagnostics et tests qui inspectent le nombre de doublons après `finalize()` sont adaptés au nouveau contrat.

**Test avant/après :** l'ancienne implémentation échoue sur les nouveaux contre-exemples (4 échecs parmi 21 cas) ; la version corrigée passe 21/21.

### 3. Réemploi Jacobi après une mise à jour invalide

`JacobiPreconditioner::setup()` effaçait son inverse diagonal avant de savoir si le nouveau système était acceptable. La préparation peut aussi produire `1/a_ii=Inf` lorsqu'un pivot est subnormal.

**Correction :** construire `next_inv_diag` hors état actif, vérifier dimensions, matrice et inverses finis, puis `swap()` atomique. Sur échec, conserver la précédente factorisation utilisable. Nouveau cas de régression pour inverse débordant et actualisation rejetée.

## Qualification ciblée exécutée sur cette archive

| Campagne | Résultat | Portée |
|---|---|---|
| Release, lot 1 | 12/12 CTest PASS | SST, SIMPLE, Schur, AMG, Neumann, CG homogène, transport, CSR |
| Release, lot 2 | 7/7 CTest PASS | diagnostics, Schur, anisotropie, port natif, géométrie NACA rapide |
| Pytest NACA/Schur ciblés | 18 PASS, 1 SKIP, 3 sous-tests PASS | outillage Python, sans CI globale |
| Pytest `tests/python` élargi, hors module Saturne | 946 PASS, 115 SKIP, 6 FAIL, 22 sous-tests PASS | 5 échecs OpenFOAM : `meshio` local absent ; 1 échec wheel : isolation PEP517 veut télécharger `setuptools` sans réseau |
| Wheel hors isolation (`--no-build-isolation`) | PASS | `cfdx-0.7.0-py3-none-any.whl` généré localement avec `setuptools 82.0.1` installé |
| Clang ASan/UBSan | cas CSR 21/21 et contre-exemple CG homogène PASS | chemins ciblés, pas suite sanitizers complète |
| `git diff --check` | PASS | aucune erreur d'espacement |

Pour reproduire les tests Python sans ombrage par `scripts/cfdx.py`, lancer `PYTHONPATH=src/cfdx/python python3 -m pytest -q tests/python`. Le module `test_saturne_adapter.py` nécessite `meshio` pour sa collecte. Les workflows GitHub prévoient son installation, mais elle est absente de cet environnement : ne pas modifier les oracles pour masquer le manque de dépendance. La commande `pip wheel --no-build-isolation --no-deps .` a réussi, mais le test `test_package_builds` utilise l’isolation de build qui demande un accès à PyPI.

Les deux lots CTest comprennent une répétition de `test_reusable_scaled_pressure_cg_v9` ; il s'agit de 18 tests distincts au total, dont la génération de maillage rapide. Une suite exhaustive de tous les tests du dépôt n'a pas été exécutée.

## Microbenchmark d'assemblage et de produit CSR

Programme autonome : `benchmarks/linear/bench_csr_coalescing.cpp`.

Matrice anisotrope synthétique à 36 864 inconnues, avec contributions COO volontairement répétées, même vecteur d'essai avant/après. Mesure de setup et du produit CSR, 3 répétitions séparées. Résultats médians locaux :

| Mesure | Avant | Après |
|---|---:|---:|
| Contributions mises en file | 440 832 | 440 832 |
| Coefficients CSR stockés | 440 832 | 183 552 |
| Temps de produit `A*x` | 433,36 µs | 171,19 µs |
| Gain sur ce noyau | — | ~2,53x |

Les sommes de contrôle valent `1.70604829689625` et `1.70604829689637` (écart d'arrondi dû au regroupement des sommes). La vitesse d'assemblage varie ; ce benchmark **ne démontre pas un gain de temps CFD global** ni un gain sur les matrices déjà sans doublons.

## Qualification physique NACA0012 — Re 1000, incidence 0°

Génération de trois maillages depuis la source fournie : 3 096, 12 336, 49 248 cellules (N64/N128/N256), sans cellules inversées signalées par le générateur. Solveur SIMPLE, `alpha_u=0.6`, `alpha_p=0.25`, `linear_tolerance=1e-10`, seuils de convergence non linéaire d'origine, convection UPWIND de référence.

| Indicateur | N64 | N128 | N256 |
|---|---:|---:|---|
| Itérations acceptées | 336 | 794 | 80, puis 240 (essais distincts) |
| Convergence itérative | PASS | PASS | **NON**, budgets épuisés |
| Temps solveur | 14.8082 s | 134.956 s | 60.1051 s (80 it.), 178.983 s (240 it.) |
| Cd | 0.0868919 | 0.128099 | non publiable |
| Cd pression | 0.0325429 | 0.0458963 | non publiable |
| Cd viscous | 0.0543490 | 0.0822031 | non publiable |
| Cl | -3.3032e-10 | -3.71379e-10 | non publiable |
| Continuité normalisée | 1.07666e-10 | 7.84596e-11 | 1.73322e-7 à 240 it. |

Le Cd varie d'environ 47 % entre N64 et N128 (valeur relative à N64) : **aucune indépendance au maillage démontrée**. Un résultat de pression linéaire certifié ne signifie pas que le résidu non linéaire de momentum a convergé. Sur N256, `momentum_residual_rhs_relative` descend de `1.8693e-3` à 80 itérations jusqu’à `5.89954e-4` à 240 itérations (non convergé), alors que le résidu linéaire de pression vaut `3.40449e-11` au dernier point. Il faut vérifier les erreurs de traction pression/visqueuse indépendamment, puis obtenir N256 convergé avant tout GCI. Le diagnostic `flux_velocity_mismatch_linf` doit être suivi sans l'assimiler automatiquement à une violation du bilan massique local.

## Essai de convection TVD bornée (opt-in, non promu)

Sur N64, la variante `CFDX_NACA_CONVECTION_TVD=1` a été exécutée avec les mêmes `alpha_u=0.6`, `alpha_p=0.25` et les mêmes critères. Après 335 itérations le contrôleur a correctement arrêté la solution au statut `STAGNATED` : `momentum_rhs_relative=7.25144e-5`, `continuity=3.114e-8`, `pressure_residual=4.13725e-12`. Ce cas ne valide donc **pas** la convection d’ordre supérieur ; la qualification N64/N128 ci-dessus est uniquement celle du schéma UPWIND de référence.

## Reprise et portes de qualification encore ouvertes

1. **P0 :** appliquer et qualifier le patch sur le HEAD exact visé, relancer Release, Python, DebugSanitizers et le contrôle d'audit avec bibliographie/fiches méthode ; ne pas confondre cet environnement local avec la CI GitHub.
2. **P0 :** N256 jusqu'à convergence non linéaire complète, sans modifier les tolérances ; chercher la première stagnation et examiner les contributions de quantité de mouvement, pression, flux Rhie–Chow et états rejetés. Un essai limité à 240 itérations n'est pas une preuve d'impossibilité.
3. **P0 :** convergence spatiale et traînée : confirmer pression/cisaillement par reconstruction indépendante des faces ; maillages N64/N128/N256, au besoin famille murale raffinée et extrapolation GCI lorsque les trois solutions sont physiquement convergées.
4. **P1 :** reprendre opérateurs skew/non orthogonaux avec Jacobien conservatif complet de COUPLED, vérification MMS quadratique et conservation locale face par face.
5. **P1 :** quantifier la qualité du Schur LSC/PCD/BFBt sur matrices **réelles**, comparer Vanka et Braess–Sarazin en itérations FGMRES, temps setup/apply/solve et coût mémoire.
6. **P0 :** éliminer la stagnation du NACA TVD d’ordre supérieur sans amoindrir les critères, en contrôlant la correction convective différée, les flux de face et l’assemblage du résidu physique sur l’état accepté.
7. **P1 :** sur anisotropie oblique et variable, comparer Line-GS + semi-coarsening, agrégation lissée et interpolation énergétique à tolérance de résidu vrai identique ; ne pas imposer systématiquement la variante la plus coûteuse.
8. **P1 :** SST, thermique, multiphase et Ghia/N9, avec campagnes de conservation et critères physiques indépendants.

**Important :** ce lot modifie le cœur CSR/CG et des tests, pas les seuils physiques. Les modifications peuvent améliorer la robustesse et la vitesse de certains opérateurs, mais ne représentent pas une validation générale de CFDX. Aucun push/merge GitHub n'a été exécuté.
