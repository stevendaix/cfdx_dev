@@ -0,0 +1,102 @@
# CFDX V62 — Jacobien généralisé ancré, PCD/Schur et BFS152

**Date : 13 novembre 2026.** Base : CFDX V61 reconstituée à partir du ZIP BFS V6 de référence. Il s'agit d'un lot **de diagnostic et de recherche**, et NON d'une qualification industrielle. Les sources radiatives HO–LO/DOM V57 fusionnées en V60 sont conservées, mais n'ont pas été modifiées dans ce lot.

## Verdict vérifiable

- Le chemin Newton V62 crée **une seule perturbation déterministe de l'état physique** avant de différencier l'ensemble des colonnes du Jacobien dense `6N`. Les colonnes n'emploient donc pas des sélections indépendantes et potentiellement incompatibles de branches de limiteur. La même évaluation complète `F_phys` demeure l'oracle et le critère d'acceptation du pas ; les préconditionneurs et seuils historiques ne sont pas modifiés.
- Un test C++ à limiteur `min` avec égalité multiple, conservation de la jauge pression, `Uz=0` sur géométrie extrudée, et dérivées croisées momentum–k–omega passe. **Ce test analytique construit n'est pas le BFS physique.** Une matrice obtenue par dérivation à un état voisin n'est pas forcément un élément certifié du Jacobien généralisé de Clarke/Bouligand à l'état de base, et reste soumise aux gardes de fidélité.
- Le solveur natif a été recompilé avec le commutateur opt-in `CFDX_FPHYS_JFNK_V62_ANCHORED_JACOBIAN`. L'activation suppose les contrôles déjà requis V54, V57 et V61. Il n'est pas activé par défaut. **Le rejeu BFS152 avec ce mode a été interrompu après plus de 200 secondes, sans certificat de fidélité ni pas Newton accepté ; cet essai ne démontre donc ni acceptation ni rejet numérique d'un pas calculé.**
- Sur **les matrices de production BFS152 à la première itération** de six préconditionneurs, l'approximation sparse `S_selfp = C - D diag(A_uu)^{-1} G` a rang numérique complet (152/152). Un oracle **externe SciPy** utilisant des ILU de vitesse et de `S_selfp` atteint `||b-Ax||_2/||b||_2 = 2.2931371313e-10` en 19 itérations GMRES, avec `rtol = 1e-9`. Cet essai est un point de référence algébrique, **pas un correctif PCD du noyau de production** et **pas** un certificat de convergence physique BFS.
- Les six Schur historiques ont tous été exécutés avec des matrices `A`, des RHS `b` et des états `x0` identiques au niveau SHA-256. **Aucun ne converge physiquement le BFS152**. PCD échoue dans le solveur linéaire après 2000 itérations ; les cinq autres échouent dans les contrôles non linéaires après momentum/SST. Les références NASA/OpenFOAM, l'équivalence des fonctions de paroi et les maillages BFS712/BFS20540 n'ont pas été qualifiés.

## 1. Jacobien différentiel commun (V62)

Pour un état de base `X`, un vecteur d'échelles strictement positives et un pas `h`, V62 construit une unique configuration ancrée

`X_anchor = X + delta * diag(scales) * q`, `delta = 32h / 2^attempt`,

avec `q` déterministe (hachage entier fixe), et avec les contraintes : jauge de pression conservée, pas de `Uz` fictif lorsque tout le champ `Uz` vaut exactement zéro, `k` et `omega` perturbés vers leurs valeurs admissibles. `attempt` est borné à 12 et l'absence d'état admissible déclenche une erreur. L'évaluation du résidu complet est recalculée à `X_anchor`.

Chaque colonne de la matrice utilise ensuite **le même** `X_anchor`, avec les quatre évaluations `F(X_anchor ± he_j)` et `F(X_anchor ± h e_j/2)` pour extrapolation Richardson et contrôle des dérivées latérales. Le test de déplacement physique `||F(X_anchor)-F(X)||` reprend une garde existante (pas d'assouplissement). Une matrice physiquement trop éloignée, dont les colonnes sont non fiables ou dont la vérification Arnoldi/true residual échoue, doit rester rejetée. La recherche de pas n'évalue que le **véritable** `F_phys(X + alpha d)` (pas un résidu ancré substitué).

**Limites mathématiques :** la perturbation d'ancrage ne garantit pas que toutes les faces SOU/MUSCL ou contraintes `min/max` soient dans un secteur différentiable. La sélection n'est pas une construction formellement prouvée de toutes les matrices de Clarke ; sur une vraie fonction non lisse, la fidélité de `J_anchor v` à la dérivée directionnelle de `F_phys(X)` peut échouer. La matrice dense exige `O((6N)^2)` en mémoire et environ `4*6N` évaluations physiques, hors tests et itérations Krylov ; elle n'est donc pas envisageable telle quelle pour BFS20540, où un opérateur local assemblé ou une différentiation directionnelle structurée est requis.

## 2. PCD, Schur capturé et oracle indépendant

Référence : système couplé natif `A=[U G; D C]` de taille `4N=608`, avec `N=152` et `7080` non-zéros. Le bloc pression compte 152 inconnues. L'oracle construit `S_selfp=C-D diag(U)^-1 G` avec **le C et la jauge originaux**, refuse les diagonales vitesse non positives, la matrice non finie, les six entrées divergentes et les opérateurs singuliers. Il n'inverse pas exactement le grand bloc vitesse dans le mode ILU.

| Indicateur mesuré sur BFS152 (4N gelé) | Résultat |
|---|---:|
| Six captures `A,b,x0` | SHA-256 identiques |
| Taille `A` ; coefficients non nuls | `608` ; `7080` |
| Taille et rang numérique `S_selfp` | `152` ; `152` |
| Non-zéros `S_selfp` ; conditionnement 2-norme estimé | `1212` ; `2.01665e5` |
| GMRES bloc avec LU exact sur U et Schur approché | 15 itérations ; `3.0931258e-10` vrai résidu relatif |
| GMRES bloc avec ILU(U) + ILU(S) (`drop=1e-3, fill=5`) | 19 itérations ; `2.2931371e-10` vrai résidu relatif |
| Convergence PCD natif | NON |

La référence utilise `scipy.sparse.linalg.gmres` et une **préconditionnement gauche**, différent de l'algorithme natif FGMRES et des cycles AMG de CFDX. Le calcul du vrai résidu emploie exclusivement l'original `b-Ax` ; le résultat des itérations ou la norme de Hessenberg ne sert pas de preuve autonome. Il ne faut pas attribuer ces 19 itérations à une méthode PCD corrigée.

Un essai auxiliaire de régularisation `Fp+tau Mp` montre une amélioration possible du conditionnement de Fp mais aucune proximité suffisante avec l'inverse de Schur réel (le meilleur défaut Schur après rééchelonnement reste de l'ordre de `0.66` dans les paramètres balayés). Aucune régularisation arbitraire de Fp ou condition de pression auxiliaire n'est promue en production.

## 3. Six Schur de production, protocole commun

Le script natif `scripts/run_bfs_six_schur_v37.py` a été relancé sur le même `bfs152_v60_exact.cfdx.h5`, avec la même géométrie, 152 cellules, 634 faces, mêmes valeurs physiques et commandes `--iterations 1 --experimental-sst`, en activant la capture de la **matrice 4N originale**. Le comparateur et les empreintes SHA-256 confirment l'égalité de `A,b,x0` sur les six méthodes.

| Préconditionneur de production | Verdict observé |
|---|---|
| PCD | échec `KRYLOV_SCHUR` ; 2000 itérations, vrai résidu relatif `6.274552e-3` |
| LSC | `POST_MOMENTUM_NOT_CERTIFIED` |
| BFBt | `POST_MOMENTUM_NOT_CERTIFIED` |
| BlockLocal | `POST_MOMENTUM_NOT_CERTIFIED` |
| SIMPLE-Schur | `POST_MOMENTUM_NOT_CERTIFIED` |
| SIMPLEC-Schur | `POST_MOMENTUM_NOT_CERTIFIED` |

Les captures ne démontrent pas une convergence de l'ensemble `U,p,k,omega`. Le BFS152 est initialisé uniformément ; les fonctions de paroi SST OpenFOAM ne sont pas entièrement actives et vérifiées dans cette expérience (`SST molecular boundary diffusion: legacy`). Aucun résultat `Cf`, `Cp`, `x_r/h` n'est donc présenté comme validation NASA.

## 4. Reproduction

Prérequis : Clang/GCC C++17, CMake, Ninja, bibliothèque de construction CFDX, Python 3, NumPy et SciPy **uniquement pour l'oracle externe**.

```bash
cmake -S . -B build-v62 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON \
  -DCFDX_ENABLE_RADIATION_RESEARCH_V57=ON
cmake --build build-v62 --target \
  cfdx_production_solver test_common_active_set_anchor_v62 \
  test_pcd_schur test_pcd_pressure_inflow_boundary \
  test_frozen_physical_jacobian_v54 test_global_physical_newton_v52
ctest --test-dir build-v62 --output-on-failure \
  -R 'test_common_active_set_anchor_v62|test_pcd_v62_selfp_reference|test_pcd_schur|test_pcd_pressure_inflow_boundary|test_frozen_physical_jacobian_v54|test_global_physical_newton_v52'
# L'oracle de matrice capturée exige les six paires originales MatrixMarket :
python3 scripts/pcd_v62_selfp_fgmres_reference.py \
  --captures path/to/six_original_captures --output /tmp/cfdx-v62-selfp.json
# Les options Newton sont forensiques et NE CONVERGENT PAS le BFS par défaut :
CFDX_FPHYS_CANDIDATE_V51=1 CFDX_FPHYS_JFNK_V52=1 \
CFDX_FPHYS_JFNK_V52_FREEZE_PRESTEP=1 CFDX_FPHYS_JFNK_V53_BLOCK=1 \
CFDX_FPHYS_JFNK_V52_KRYLOV=2 CFDX_FPHYS_JFNK_V54_FROZEN=1 \
CFDX_FPHYS_JFNK_V57_TRUTH=1 CFDX_FPHYS_JFNK_V61_RICHARDSON=1 \
CFDX_FPHYS_JFNK_V62_ANCHORED_JACOBIAN=1 \
  build-v62/cfdx_production_solver \
  --mesh path/to/bfs152_v60_exact.cfdx.h5 \
  --output-dir /tmp/cfdx-bfs152-anchor --iterations 1 --experimental-sst
```

**Compilation diagnostique `-O0` pour l'unité principale : aucune qualification de performance Release.** L'essai dense complet `6N` sur BFS152 peut être coûteux ; une interruption ou un timeout signifie absence de résultat, **pas** convergence implicite.

## 5. Suite prioritaire (P0 encore ouvert)

1. Construire une linéarisation active-set **locale, analytique et par faces**, raccordée au flux convectif primal ainsi qu'aux dérivées de `rAU`, Rhie–Chow, gradient de pression, viscosité turbulente et contraintes SST. Contrôler l'opérateur complet contre les résidus physiques `F(X±hv)`, sans supposer qu'un sous-gradient ancré sera automatiquement fidèle.
2. Instrumenter les blocs `d(momentum)/dk`, `d(momentum)/domega`, `d(k)/dU`, `d(omega)/dU`, les changements des fonctions de paroi, les coefficients de face et le transport spatial. Distinguer non-différentiabilité vraie, bruit arithmétique et erreur d'assemblage.
3. Étudier le défaut de commutation PCD sous la jauge et les BC réelles ; comparer l'action `Kp^-1 Fp Mp^-1` orientée à une solution algébrique Schur vraie. Évaluer `selfp` comme **septième référence**, non comme remplacement clandestin des six variantes.
4. Rejouer les six Schur sur un `Amat` **identique, fidèle et gelé**, monitorer les résidus linéaires vrais par direction et l'efficacité du préconditionneur ; ne permettre le passage à Newton qu'après certificat de descente de `F_phys` sur un état accepté.
5. Ensuite seulement, lancer BFS712 et BFS20540 avec parois SST OpenFOAM comparables, références expérimentales Driver–Seegmiller (`x_r/h≈6.26±0.10`) et études Cf/Cp/profils, sans desserrer les tolérances. Qualifier les régressions NACA et les autres cas CFDX en parallèle.

## 6. Patches et modes d'exécution

Le ZIP source V61 livré lors du lot précédent perdait les permissions exécutables de **14 scripts** (radiation, audit BFS et tests). Pour ne pas masquer le défaut de packaging, la livraison comporte trois chemins distincts : `CFDX_V62_PATCH_CUMULATIF_BFS_V6.patch` depuis le ZIP d'origine exact ; `CFDX_V62_PATCH_INCREMENTAL_V61_V62.patch` depuis **le ZIP de sources V61** (restaure aussi ces 14 modes) ; `CFDX_V62_PATCH_INCREMENTAL_V61_GITTREE_V62.patch` pour un **arbre Git V61** dont les modes exécutables étaient déjà corrects. Les deux reconstitutions propres doivent être comparées par SHA-256 et par permissions. Ne pas confondre les deux patchs incrémentaux.

**Statut logiciel :** **16/16 CTest ciblés réussis**, en deux campagnes (11 + 5), avec compilation du binaire natif de production en mode diagnostique `-O0` ; pas de CI complète, pas de push GitHub, pas de qualification GPU/MPI ni de campagne industrielle. Toute preuve doit être reproduite après fusion sur la branche Git réellement visée.
