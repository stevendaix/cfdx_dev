# CFDX — correctifs anisotropes et audit des tests (6 novembre 2026)

## Périmètre et provenance

**Base de travail :** archive fournie `cfdx_dev-master(2).zip`, dont le commentaire
identifie la révision `9d4d132306b222fc3ccc0c7d4cbd0f5125801965`.
**Étude fournie :** `test.txt.txt` (933 lignes). Ses constats relatifs au commit
`fe9606a` ne sont pas automatiquement transposables à cette archive.

Les fichiers de cette archive ont été audités et modifiés localement. **Aucune
PR, aucun push, aucun merge ni aucun résultat GitHub Actions ne sont revendiqués.**
Les journaux récents de la CI GitHub n'étaient pas accessibles avec les moyens
d'accès disponibles ; la vérification externe sur le HEAD exact reste nécessaire.

## Défauts reproduits et correctifs

1. **Rafraîchissement AMG conservant des lignes devenues incorrectes.**
   `AlgebraicLinePreconditioner::update_values()` vérifiait le rapport des
   couplages négatifs maximal/minimal, mais un couplage extrêmement faible
   pouvait conserver ce rapport même lorsque la direction physique dominante
   avait tourné de 90 degrés. Le nouveau contrôle certifie, pour chaque
   cellule couverte par une ligne, la part initiale du couplage dominant
   représentée par cette ligne :

   \[
   d_i(A,L)=\frac{\max_{j\in L(i)}(-a_{ij})_+}
                  {\max_{j\ne i}(-a_{ij})_+},\qquad
   d_i(A_{\rm nouveau},L)\ge 0{,}90\,d_i(A_{\rm initial},L).
   \]

   Une rotation forte entraîne un refus explicite du *numeric refresh* :
   l'appelant doit reconstruire les lignes. L'ancien préconditionneur demeure
   utilisable après refus. Cette règle est prudente : elle peut augmenter le
   nombre de reconstructions en échange d'une meilleure protection contre le
   recyclage de structures inadaptées. Elle ne remplace pas un détecteur
   géométrique général des directions anisotropes.

2. **`setup()` destructif en cas d'échec.** Les objets
   `AlgebraicLinePreconditioner` et `AnisotropicLineAMGPreconditioner`
   conservent maintenant leur état précédent si la préparation du nouveau
   système échoue. La nouvelle hiérarchie Line-AMG est préparée puis échangée
   seulement après succès. Le coût de création n'est pas réduit : c'est une
   correction de contrat transactionnel.

3. **Préconditionnement en place (`rhs` et `out` aliasés).** Dans
   `AnisotropicLineAMGPreconditioner::apply`, le second membre était détruit
   par le premier lissage. La copie du second membre est désormais effectuée
   **uniquement en cas d'alias** ; les appels Krylov normaux n'en paient pas
   le coût.

4. **Allocations répétées en application du Line-AMG.** Les substitutions
   de Thomas sont calculées directement dans le vecteur destination, les
   produits matrice-vecteur utilisent un buffer préalloué local à l'appel et
   le buffer de correction AMG sert aussi au post-lissage. Cela préserve
   l'ordre de complexité et l'état numérique retourné, sans état mutable
   partagé entre appels concurrents.

5. **Mesure mémoire du benchmark N15.** Sur certains chemins `fork/exec`,
   `wait4.ru_maxrss` reflétait une borne héritée du parent, et le test
   comparait un tableau virtuel pas nécessairement résident. Sous Linux,
   le runner échantillonne désormais `VmHWM` de `/proc/<pid>/status`
   **après vérification de l'identité du processus exécuté**. Un échantillon
   manquant est `unavailable` et non une valeur inventée. Le test alloue et
   touche réellement 80 Mio. La métrique est un pic RSS post-exec relevé par
   échantillonnage, et non une mesure idéale pour les processus très courts.

6. **Test NACA instable sous charge Python complète.**
   `test_naca_controlled_campaign.py::test_nonzero_exit_is_failure` lançait
   un interpréteur Python avec seulement deux secondes de budget, alors que
   l'oracle devait simplement vérifier qu'un code de sortie 4 ne devient
   jamais PASS. Dans la campagne globale, le démarrage dépassait parfois le
   budget : l'enregistrement portait `reason=timeout` et ne possédait pas de
   champ `exit_code`. Il s'agissait d'une **fixture instable**, pas d'un
   défaut physique du solveur NACA démontré. Un petit script shell POSIX
   retourne désormais directement 4 ; l'oracle vérifie explicitement le
   motif `nonzero exit` et le code, tandis que la vérification des vrais
   timeouts reste un test indépendant. La CI de ce dépôt exécute ces tests
   sous Linux.

## Tests numériques ajoutés

`tests/unit/test_anisotropic_line_integrity.cpp` comporte quatre scénarios
indépendants :

- rotation de la direction dominante en présence de couplages résiduels très
  faibles, avec refus du rafraîchissement et préservation de l'ancien état ;
- échec de `setup()` conservant la factorisation par lignes ;
- échec de `setup()` conservant la hiérarchie Line-AMG ;
- égalité des applications Line-AMG en place et hors place, à 3e-12 près.

Les trois premières familles de défauts (rotation, setup destructif,
application aliasée) ont été constatées sur la version de départ à l'aide de
tests ciblés avant l'introduction des correctifs. Le quatrième test de setup
Line-AMG a été ajouté lors de l'extension transactionnelle.

## Bilan de compilation et de tests locaux

| Porte | Résultat | Portée |
|---|---|---|
| CTest Release étendu | **16/16 PASS** | Tests ciblés AMG, Krylov, Schur, Galerkin, Neumann anisotrope, solveur stationnaire, intégrité des lignes |
| Clang DebugSanitizers (ASan + UBSan) | **3/3 PASS** | Intégrité des lignes, Line-AMG, adaptation native anisotrope PETSc/HYPRE ; aucune erreur mémoire relevée dans ces exécutions |
| Tests Python N15 + NACA après correction | **8/8 PASS** | RSS N15, comportement `exit 4` NACA, refus des timeouts, métriques |
| Campagne Python initiale partielle | **945 PASS, 115 SKIP, 7 FAIL** | Les 7 échecs ont été diagnostiqués, mais cette campagne précède la correction N15 |
| Suite Python relancée après corrections, sans dépendances indisponibles | **945 PASS, 115 SKIP, 0 FAIL, 22 sous-tests PASS** | NACA et N15 passent dans la campagne complète exécutable, durée 178,97 s |
| Compilation Release totale | **NON TERMINÉE** | Plusieurs centaines de cibles ; aucune prétention de CI globale verte |
| CI GitHub du HEAD exact | **NON VÉRIFIÉE** | Journaux/jobs non accessibles depuis cette session |
| NACA N128/N256, SST, 3D et efforts | **NON QUALIFIÉS ICI** | La réussite des tests linéaires ne prouve pas la convergence CFD physique |

Répartition des 7 échecs initiaux Python : cinq tests OpenFOAM manquaient
`meshio`, le test packaging requérait une dépendance distante inaccessible
(réseau de cet environnement), et le test N15 contenait une hypothèse
fragile sur le RSS et les pages mémoires touchées. Le test N15 a été corrigé.
Les six échecs dus à l'environnement **ne constituent pas** des défauts CFDX
établis : la CI les exécute normalement après installation des dépendances.
Le test NACA supplémentaire a été découvert en élargissant la campagne, puis
stabilisé sans relâcher le contrat d'échec du lanceur. La dernière passe
Python a bien été relancée après cette correction : zéro échec sur la sélection
exécutable. Les modules non collectables faute de dépendances ne sont pas
comptés comme exécutés. Le linter/formatter Ruff n'était pas installé dans
cet environnement : sa porte CI reste à vérifier sur le runner officiel.

La compilation C++ émet aussi des avertissements de conversion de types déjà
présents dans la base, et au moins un avertissement d'optimisation GCC associé
à `MatrixFreeVcyclePreconditioner::smooth`. Leur cause n'est pas prouvée : ils
ne sont pas qualifiés ici de défauts mémoire et ne doivent pas être corrigés
à l'aveugle. La campagne complète de tests mémoire reste nécessaire.

## Mesure de performance : ne pas confondre microbenchmark et simulation

Le programme optionnel `cfdx_benchmark_anisotropic_line_apply` mesure la
préparation et le coût moyen d'une application du préconditionneur, puis le
résidu relatif vrai `||b-Az||₂/||b||₂` et un checksum. Une application de
préconditionneur **n'est pas** censée résoudre entièrement `Az=b`, d'où le
résidu non nul, identique entre les versions comparées.

Comparaison locale avec les mêmes options de compilation et les mêmes
matrices, 60 applications répétées et deux exécutions par version :

| Système | Version initiale, temps par application | Corrigée | Remarque |
|---|---:|---:|---|
| 48² = 2 304 cellules | 0,801–0,817 ms | 0,880–0,901 ms | Régression locale, donc aucune accélération universelle revendiquée |
| 64² = 4 096 cellules | 1,523–1,631 ms | 1,445–1,536 ms | Gain modéré, sensible aux variations de mesure |
| 96² = 9 216 cellules | 3,708–3,723 ms | 3,596–3,624 ms | Gain modéré, sensible aux variations de mesure |

Le checksum et le résidu relatif retournés ont été identiques entre les
deux versions pour chaque taille. Ces chiffres ne sont **pas** des gains de
temps de résolution NACA ou CFD globale. Ils ne couvrent pas non plus la
préparation répétée ou les effets cache/NUMA sur des grands problèmes.

Exécuter un benchmark local reproductible :

```bash
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_ENABLE_MPI=OFF
cmake --build build-release --target cfdx_benchmark_anisotropic_line_apply
./build-release/cfdx_benchmark_anisotropic_line_apply 64 60
./build-release/cfdx_benchmark_anisotropic_line_apply 96 60
```

## Ce que recommande l'étude et ce qui reste réellement à développer

L'étude jointe distingue utilement la qualité des opérateurs géométriques, le
conditionnement linéaire et la convergence non linéaire. En particulier :

- **Amélioration livrée :** intégrité des lignes lors des updates numériques,
  conservation de la hiérarchie validée, aliasing, allocations auxiliaires,
  instrumentation de benchmark, tests de régression et tests instrumentés.
- **Non implémenté :** semi-coarsening guidé par les directions physiques,
  interpolation à distance algébrique et/ou minimisation énergétique,
  détecteur de directions sur maillages courbes/non structurés, Line-GS
  multiplicatif et Schwarz recouvrant. La factorisation par lignes actuelle
  reste un solveur de blocs **disjoints**.
- **À qualifier avant activation globale :** nouveau dispatch adaptatif, car
  la qualification fournie indique que les lignes seules peuvent échouer à
  10³ tout en réussir à 10⁶. Une activation automatique non certifiée serait
  plus risquée que le comportement actuel.
- **Séparé du préconditionnement :** corrections non orthogonales des équations
  de quantité de mouvement, skewness, Rhie–Chow, gradients aux faces, bilans
  masse/énergie et forces. Changer l'AMG ne les corrige pas.
- **Physique et maillage :** qualificatifs NACA N64/N128/N256, convergence
  stricte, critères de conservation, Cd/Cp/Cf/Cl/Cm, GCI, couches limites
  courbes fortement étirées, SST/CHT et solveurs couplés.

Pour l'étape suivante, comparer sur matrices réelles N128/N256 un mélange
FGMRES + Line-AMG, Boomer, SA, puis mesurer les rapports de réduction des
**vrais** résidus physiques, les temps de préparation et les durées complètes.
Ne pas forcer les solveurs potentiellement non symétriques dans CG sans preuve
SPD. Pour les opérateurs de diffusion non orthogonaux, développer un oracle
manufacturé affine/quadratique et un oracle de conservation aux faces avant
toute activation par défaut.

## Reproduction locale

```bash
# Depuis la racine CFDX modifiée :
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_ENABLE_MPI=OFF
cmake --build build-release --target test_anisotropic_line_integrity test_v24_anisotropic_line_amg test_petsc_hypre_native_port_anisotropy
ctest --test-dir build-release -R 'test_anisotropic_line_integrity|test_v24_anisotropic_line_amg|test_petsc_hypre_native_port_anisotropy' --output-on-failure
PYTHONPATH=src/cfdx/python python -m pytest -q tests/python/test_n15_benchmark.py
# CI complète sur runner doté des dépendances décrites dans .github/workflows/cfdx-ci.yml.
```

Bibliographie de l'étude : Mavriplis (1998, multigrid anisotrope), Brandt et
coll. (2015, distance algébrique), Schröder (2012, interpolation anisotrope),
et les documentations officielles PETSc PCGAMG et HYPRE BoomerAMG.

**Verdict :** les correctifs locaux sécurisent des contrats identifiés et
passent les portes ciblées. Ils ne démontrent ni un meilleur solveur universel
ni la disparition de tous les échecs CI : ces conclusions nécessitent les
journaux GitHub et les campagnes physiques complètes.