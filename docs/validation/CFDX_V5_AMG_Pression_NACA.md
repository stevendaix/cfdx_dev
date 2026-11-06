# CFDX V5 — AMG CSR et replay de pression NACA

Date : 5 novembre 2026. Base : archive CFDX V4. Portage **natif C++17** : aucune dépendance binaire PETSc ou HYPRE ; aucun fichier source de ces bibliothèques n'a été recopié.

## Objectif et modifications

La V4 sépare déjà le calcul symbolique et numérique de la projection de Galerkin : \(A_c=P^T A P\). La V5 ajoute :

1. Une adjacence \(P^T\) stockée en **trois tableaux CSR contigus** (offsets, indices de lignes fines, positions dans \(P\)), et non un tableau de petits vecteurs alloués séparément. L'ordre des contributions reste stable.
2. Une accumulation des poids d'interpolation **à marqueurs creux** pour les AMG de type Boomer et à agrégation lissée, sans `std::map` par ligne. Les poids sont triés par colonne avant l'écriture de la prolongation.
3. Une voie `GalerkinCSRProduct::numeric_into_existing()` : actualisation uniquement des valeurs d'un opérateur grossier dont la **structure CSR élaguée est inchangée**. Les changements de structure ou les données invalides provoquent un refus sans mutation. Une copie de l'ancienne hiérarchie est néanmoins maintenue pour l'actualisation transactionnelle globale.
4. Des scripts **à échec explicite** pour comparer V4 et V5 sur les captures NACA réelles, choisir expérimentalement un plafond d'interpolation et comparer les historiques physiques des premières itérations.

L'optimisation ne doit pas être décrite comme une copie du code PETSc/HYPRE. Elle s'inspire des principes publiés de séparation symbolique/numérique, réemploi de structure et maîtrise du grossissement, tout en conservant la logique CFDX.

## Résultats contrôlés

### Linéaire — 3 corrections successives N128 et N256

- N128 : 12 336 inconnues ; N256 : 49 248 inconnues. Six captures au total issues de l'assemblage de pression de production.
- Le solveur est **FGMRES avec les AMG natifs** ; le test vérifie \(\|b-Ax\|\) sur la matrice d'origine, tolérance relative maximale \(10^{-10}\).
- Boomer normal, Boomer avec 4 coefficients d'interpolation, et agrégation lissée : les nombres d'itérations et les résidus vrais observés sont **identiques V4/V5**.
- Les valeurs de convergence de la première matrice ont aussi été vérifiées sur une solution SciPy indépendante : **6 méthodes N128** et **5 méthodes N256**, toutes avec une erreur relative en pression inférieure à \(5,1\times10^{-11}\). GMRES+ILU0 et préconditionnement par lignes seul ne sont pas homologués sur N256 à 300 itérations.

### Temps médians de 3 répétitions appariées (ms)

| Configuration | V4 | V5 | Rapport V4/V5 |
| --- | ---: | ---: | ---: |
| N128 Boomer sans limite | 155,10 | 157,74 | 0,983 |
| N128 Boomer limite 4 | 152,77 | 144,60 | 1,056 |
| N128 agrégation lissée | 202,92 | 200,97 | 1,010 |
| N256 Boomer sans limite | 722,80 | 777,28 | 0,930 |
| N256 Boomer limite 4 | 695,57 | 707,33 | 0,983 |
| N256 agrégation lissée | 1077,16 | 995,03 | 1,083 |

Les mesures sont **mixtes et sensibles au bruit** ; la V5 n'apporte pas un gain universel. Il est prématuré de changer les choix AMG par défaut sur ces seules données.

Un seul relevé de mémoire maximale résidente sur N256 n'indique **pas de réduction nette** : Boomer environ 150 Mo dans les deux versions ; agrégation lissée environ 120 Mo. Il s'agit de mémoire de processus, non d'un profil d'allocations internes.

### Balayage expérimental du plafond Boomer

Sur trois captures par maillage et trois répétitions :

- **N128** : plafond 3, médiane 129,91 ms, contre 162,44 ms illimité (1,25 fois plus rapide dans cette série), complexité d'opérateur réduite de 4,56 à 3,23 ; 32 itérations contre 29 pour les trois corrections.
- **N256** : le meilleur plafond testé (4) donne 671,27 ms contre 702,10 ms (1,046 fois). Le seuil de recommandation du tuner est 1,10 fois : **aucun changement recommandé**.
- Il s'agit d'un **outil de sélection hors ligne** : il ne change pas les réglages du solveur de production. Une sélection fondée uniquement sur les résidus linéaires ne qualifie pas encore les forces et la convergence physique.

### Physique — début de trajectoire CFD NACA

Avec la même configuration `CFDX_NACA_PRESSURE_SOLVER=fgmres_boomer`, la même entrée N128/N256 et trois itérations externes :

- Nombre d'itérations pression N128 : **9, 8, 9** ; N256 : **9, 9, 9**, identiques V4/V5.
- Écarts maximaux **exactement nuls aux chiffres enregistrés** pour continuité, résidu de quantité de mouvement, déséquilibre de masse global, déséquilibre local et cohérence flux/vitesse.
- **La convergence NACA complète, les forces finales, la sensibilité au maillage, le GCI et les couplages SST/thermiques ne sont pas qualifiés par cette campagne.**

## Reproduction

Construire les tests V5 en Release et exécuter :

```bash
ctest --test-dir build --output-on-failure -R 'test_galerkin_csr_symbolic_numeric|test_amg_cap_sweep_fail_closed|test_naca_v5_benchmark_gates|test_amg_preconditioner_qualification|test_neumann_anisotropic_pressure_qualification'
```

Les captures requises doivent correspondre aux six manifests `N128_portable*.manifest` et `N256_portable*.manifest` et aux fichiers binaires référencés. On ne doit jamais comparer deux exports d'itérations ou de configurations différentes.

```bash
python3 tools/benchmarks/compare_naca_amg_v5.py \
  --base /path/to/v4/cfdx_replay_pressure_sequence \
  --candidate /path/to/v5/cfdx_replay_pressure_sequence \
  --captures /path/to/captures --repetitions 4 --output paired.csv

python3 tools/benchmarks/qualify_naca_amg_cap_sweep.py \
  --replay /path/to/v5/cfdx_replay_pressure_sequence \
  --captures /path/to/captures --repetitions 4 \
  --raw cap_raw.csv --summary cap_summary.csv

python3 tools/benchmarks/compare_naca_outer_histories.py \
  --baseline /path/to/V4/N128.outer_history.csv \
  --candidate /path/to/V5/N128.outer_history.csv \
  --expected-iterations 3 --max-absolute 1e-12
```

Les replays s'arrêtent avec un code d'erreur en cas de résidu insuffisant, fichier manquant, manifest incorrect, résultat non fini ou ordre de séquence invalide. Les tests Python incluent des contre-exemples qui doivent échouer.

## Critères d'intégration et travaux suivants

- Conserver les garde-fous de validité et la transactionnalité. La mise à jour `numeric_into_existing` échoue si le graphe CSR change, même de manière infime, afin de ne jamais recycler silencieusement un graphe périmé.
- Ne pas généraliser le plafond 3/4 sur la seule base de trois captures et de quelques répétitions ; le calibrer hors ligne sur davantage de matrices, de directions anisotropes et d'itérations réellement convergées.
- Étudier la réutilisation d'espaces de travail pour les produits Galerkin : des tableaux numériques `AP` et `P^TAP` sont encore alloués pendant les mises à jour.
- Définir et mesurer le **temps jusqu'à convergence CFD complète**, le déséquilibre de masse, les forces `Cd/Cl/Cm`, et la convergence spatiale N64/N128/N256 ; valider aussi N9, SST, thermique et rayonnement.
- Pour les cas parallèles, contrôler les coûts de communication du Krylov (réductions globales) ; aucun gain MPI/GPU n'est mesuré ici.

### Références algorithmiques (aucun code externe copié)

- PETSc, [MatPtAP](https://petsc.org/main/manualpages/Mat/MatPtAP/) et [MatProductSymbolic](https://petsc.org/main/manualpages/Mat/MatProductSymbolic/).
- HYPRE, [BoomerAMG documentation](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html).
