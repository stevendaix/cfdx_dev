# CFDX V24 — SST JFNK, préconditionnement spatial couplé et audit des dérivées réelles

**Date :** 10 novembre 2026. **État : expérimental, désactivé par défaut, non qualifié physiquement.**

## 1. Objet et état de départ

En V23, un préconditionneur de blocs locaux 2×2 capture les dérivées croisées des sources SST. Sur le BFS réduit de 712 cellules, son résidu directionnel GMRES à l'itération SST 10 était `0.74693614522`, supérieur au seuil **linéaire** de `0.25`, et son pas Newton était rejeté. Le candidat V22 à gauche/Jacobi conservait `0.24793702144`, puis un mérite SST physique réassemblé `0.8899074514`. **Ni V22 ni V23 ne convergent complètement le BFS.**

L'objectif V24 est d'améliorer l'approximation *préconditionneur* sans modifier l'opérateur physique : remplacer l'action locale pure par un balayage sur la connectivité complète des matrices FV de `k` et `omega`. Mesurer les blocs de Jacobien réels sur 357, 363, 426, et ne jamais promouvoir un pas qui échoue au test de mérite.

## 2. Code introduit

### 2.1 Préconditionneur spatial à blocs

`src/cfdx/physics/sst_spatial_block_preconditioner_v24.h` introduit `SSTSpatialBlockPreconditionerV24`. Il réutilise les matrices CSR `A_k`, `A_omega` réellement assemblées par `sst_reassembled_equations_v22`, et les **inverses locaux V23 déjà calculés**. Dans les inconnues logarithmiques `z_j=log(x_j/x_ref)` avec les facteurs de ligne gelés `s_i`, chaque coefficient de transport hors diagonale correspond à :

`P_ij = A_ij * x_j / s_i`, avec `s_i = frozen_residual_scales_i * frozen_left_diagonal_i`.

Le bloc diagonal `P_ii` conserve les dérivées croisées des sources V23. Une passe bloc Gauss–Seidel croissante puis décroissante, initialisée à zéro à **chaque application**, forme une action `P^{-1}` **linéaire, figée et déterministe** utilisable en préconditionnement droit de GMRES. La matrice complète n'est ni substituée à `J*v`, ni utilisée pour les bilans physiques ; l'opérateur Newton reste obtenu en réassemblant les équations `k` et `omega` de production à chaque essai.

La construction CSR rejette dimensions incohérentes, coefficients non finis, échelles non positives et débordements ; le solveur conserve les autres candidats si cette nouvelle proposition n'est pas valable. L'option `sst_spatial_block_preconditioner_v24_experimental` est désactivée par défaut, exige JFNK V22 et est indépendante de l'activation de la proposition V23. Elle est accessible par `numerics.raw_settings` dans le solveur de production.

**Parallélisme :** version de recherche CPU séquentielle, O(nnz) par application. Ne pas annoncer GPU / MPI / multi-GPU sans une implémentation adaptée (multi-coloring, Schwarz, relaxation locale asynchrone ou multigrille sur les matrices de transport) et mesures séparées.

### 2.2 Comparaison protégée des directions

À la première puis à chaque dixième itération interne, V24 compare facultativement Picard, V22 Jacobi, V23 blocs locaux (si demandé) et V24 spatial. Leurs éventuels pas passent les **mêmes** contrôles de positivité, actualisation de `omega` et de la production aux parois, et mérite SST physique réassemblé. Aucun seuil physique, temps d'arrêt, schéma de transport ou conditions aux limites n'est modifié. Le journal `SST_JFNK_V22` indique `spatial_probe`, `spatial_proposed`, `spatial_linear_true_relative`, `spatial_offcell`, `selected`.

### 2.3 Jacobien analytique approximatif contre Jacobien directionnel réel

`sst_local_jacobian_audit_v24` compare chaque bloc V23 au bloc réel du résidu FV reconstitué par différences finies symétriques en **log(k),log(omega)**, sur les cellules ciblées. En diagnostic `CFDX_SST_JACOBIAN_AUDIT_V24=1`, le solveur trace les cellules 357, 363 et 426 aux itérations 1 et 10 sur un maillage au plus 5 000 cellules. La mesure est **désactivée par défaut**, non intrusive, et son coût exclut toute mesure de performance dans ce mode.

## 3. Essai réel contrôlé, BFS 712 cellules

Deux HDF5 ont été comparés au niveau de leurs tableaux de topologie/géométrie : **identiques**. Les champs `case_setup_json` deviennent identiques après retrait de la seule option V24. Le cas conserve 10 itérations internes SST uniquement pour le *diagnostic* ; tolérance SST **1e-8 inchangée**, aucun succès de l'itération extérieure.

| Première tentative, itération SST 10 | V23 local | V24 spatial |
|---|---:|---:|
| Coefficients hors cellule dans P | 0 | 5 256 |
| Erreur directionnelle vraie GMRES (proposition) | 0.74693614522 | **0.16613194317** |
| Seuil directionnel **linéaire** | 0.25 | 0.25 |
| Direction préconditionnée retenue dans l'acceptation | NON | NON |
| Direction effectivement choisie | Newton–Jacobi V22 | Newton–Jacobi V22 |
| Mérite SST **physique** réassemblé | 0.8899074514 | 0.8899074514 |
| État extérieur convergé | NON | NON |

**Interprétation :** la résolution du système Newton préconditionné s'améliore, mais la globalisation rejette le candidat spatial malgré son résidu directionnel `0.166<0.25`. L'échec n'est pas caché : V22 gagne le choix de mérite, l'itération extérieure est rejetée. Une norme linéaire relative `L2` à 0.166 n'assure pas la décroissance d'un mérite physique `L∞` à échelles spatiales variables.

**Test algébrique distinct :** sur un réseau synthétique de 12 cellules couplées, la norme L2 du résidu du préconditionneur gelé passe de `4.39583` (blocs locaux) à `0.863666` (spatial). Cette démonstration du préconditionneur *approché* ne certifie aucune solution RANS.

## 4. Audit local du vrai Jacobien SST

L'écart relatif indiqué est le maximum absolu des différences des 4 coefficients, rapporté au maximum absolu des 4 coefficients du bloc vrai ; toutes les lignes utilisent les échelles **figées identiques** au JFNK V22.

| Cellule | Écart premier pas SST | Écart itération SST 10 |
|---|---:|---:|
| 357, propriétaire paroi | 11.54197 % | **11.54218 %** |
| 363, cellule intérieure | 0.37878 % | **0.19258 %** |
| 426, production turbulente sensible | 0.000378 % | **0.000472 %** |

L'erreur maximale de bloc local est donc observée en **357**, mais sa contrainte omega a un Jacobien réel diagonal quasi unitaire à ces états. L'écart apparaît notamment dans les dérivées de l'équation `k`, au voisinage de la paroi. Pour aller plus loin, il faudra décomposer la sensibilité des flux diffusions et de `G` à l'état de paroi et vérifier leur cohérence avec l'assemblage, sans modifier le modèle physique au seul motif d'une mauvaise convergence.

## 5. Tests réellement exécutés / non exécutés

- **10/10 CTest ciblés** (tests V17–V24, plus protocoles Python V23/V24) : réussis.
- **14/14 tests Python ciblés** V22/V23/V24 : réussis ; plusieurs sont des tests de protocole avec des fixtures synthétiques, pas des oracles physiques.
- Le solveur de production compile et se lie en Release **diagnostique** (translation unit principale à `-O0` pour limiter le coût de compilation). Il n'y a ni benchmark optimisé ni résultat GPU/MPI.
- Le BFS 712 cellules a produit des traces réelles aux itérations 1 et 10. **Ce rejeu à budget réduit n'atteint pas la convergence SST ; aucune convergence globale n'est revendiquée.**
- Non exécutés ici : BFS complet 800 itérations SST avec états extérieurs convergés ; grand BFS 20 540 cellules ; toutes les variantes COUPLED Schur/PCD/AMG/MGR ; validation NASA (x_r/h, Cf, Cp, profils Reynolds) ; CI intégrale.
- Les poids de `omegaWallFunction` et la manipulation matricielle OpenFOAM Foundation restent **non équivalents numériquement** à la simple ligne propriétaire entièrement remplacée dans CFDX ; ils exigent une qualification séparée. Aucun push GitHub.

## 6. Suite prioritaire justifiée

1. Faire converger la **recherche de pas non linéaire** lorsque GMRES produit déjà un résidu directionnel satisfaisant : analyser le vrai mérite le long de la direction spatiale en logarithmes et le comportement des cellules les plus sensibles. La tolérance SST ne doit pas changer.
2. Améliorer le traitement des dérivées aux parois et des coins : distinguer contraintes propriétaire, poids par face `updateCoeffs()` et `manipulateMatrix()`, recalcule de `G`, effets de `nu_t` et diffusion. Ajouter MMS/admissibilité pour les dérivées.
3. Comparer Jacobi, blocs locaux, spatial SGS et d'autres préconditionneurs parallélisables sous une **même** politique d'essais et de convergence. Qualifier ensuite 3 états SST extérieurs sur 712 cellules avec 800 itérations au plus et toutes les portes globales.
4. Seulement après succès, grand BFS et six Schur avec mêmes données physiques, suivis de la comparaison NASA ; préserver NACA, énergie, rayonnement, sans promettre que ce lot les qualifie.

## 7. Références méthodologiques

- PETSc SNES : https://petsc.org/release/manual/snes/ — distinction résidu non linéaire, Jacobien/approximation et préconditionneur.
- OpenFOAM Foundation v11 : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html — pondérations par faces et manipulation de matrice.
