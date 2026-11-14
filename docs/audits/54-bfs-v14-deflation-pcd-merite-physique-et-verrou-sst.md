# CFDX — BFS V14 : déflation bilatérale Fp, résidu vrai, contraction physique et SST

Date : 2026-11-09. Portée : solveur de production C++ de l'archive V13 fournie,
rejeux rapides d'un BFS six blocs de 712 cellules, matrices Fp capturées.
**Aucune qualification physique Driver–Seegmiller acquise et aucune tolérance assouplie.**

## 1. Incompatibilité Fp

Pour A=Fp, W base orthonormée de ker(A^T), Z base de ker(A), le second membre
compatible est b_c=(I-WW^T)b. Le plancher vrai du second membre ORIGINAL est
||WW^Tb||/||b||. Le préconditionneur PCD peut inverser l'opérateur déflaté
A_d=A+WZ^T, sans jamais changer le système physique COUPLED. Le solveur interne
est FGMRES, avec un préconditionneur ILU0 construit sur le **système ancré
auxiliaire**; le produit Arnoldi reste A_d, et non l'opérateur ancré.

L'évaluation indépendante conserve trois résultats :
- erreur originale ||A x-b||/||b|| ;
- erreur de l'équation compatible ||A x-b_c||/||b_c|| ;
- erreur de l'opérateur déflaté ||A_d x-b_c||/||b_c||.

V14 corrige une erreur d'étiquetage : la preuve de résidu auxiliaire ne peut
plus être calculée contre la matrice d'ancrage du PRÉCONDITIONNEUR quand
l'opérateur réellement résolu est A_d.

Captures Fp réelles : 712x712, 3442 coefficients, deux modes certifiés sur
la première capture; plancher ORIGINAL 0,06279468, résidu déflaté 2,29e-9
en 25 itérations pour le rejeu. Sur BFS, 632 applications Fp déflatées acceptées
par contrôle indépendant : maximum résidu vrai déflaté 9,99e-9 ; maximum
résidu original 0,212323. **Ce succès n'est pas une convergence du FGMRES
extérieur** (arrêt temporel du diagnostic).

La base V13 dépend de la matrice *exacte*. Son empreinte est vérifiée au
chargement ; toute évolution de Fp entre états externes requiert une nouvelle
certification. Aucune mise à jour automatique des bases n'est intégrée.

## 2. Contrat non linéaire indépendant

Le garde opt-in compare les mérites issus du réassemblage complet U,p,k,omega,
nu_t. Toute tentative rejetée retourne l'état antérieur et n'entre pas au
registre `GLOBAL_CONTRACTION_COMMITTED_V14`. Une nouvelle validation exige
une contraction du mérite de 0,995 au maximum ; seuils des équations inchangés.

BFS Block-local : état 1 réellement accepté, mérite 0,527537. À l'état 2,
le bloc COUPLED est résolu (résidu vrai 9,73e-13), mais le SST retourne
`committed=0` après 800 pas ; il n'existe PAS de suite de trois états validés.
`qualify_bfs_global_contraction_v14.py` refuse donc la preuve (exit 3).

## 3. Cause de stagnation SST observée

Le diagnostic `CFDX_SST_DIAGNOSTICS=1` prouve à la deuxième itération que les
changements de champs sont devenus nuls, mais que le défaut absolu des
équations réassemblées reste supérieur au seuil : omega row-scaled
2,569e-8 pour une cible absolue 1,0e-8 à partir d'environ 300 itérations.
Aucune modification du critère n'est faite pour donner un vert trompeur.
L'arrêt V14 `SST_FP64_PRECISION_STALL_V14` détecte, sur trois itérations
identiques consécutives, l'absence de progrès FP64 avec défaut physique
non fermé : il s'arrête alors **UNCONVERGED** sans épuiser artificiellement
les 800 itérations. Sur le rejeu, l'arrêt se déclenche après 248
itérations, avec défaut absolu omega de 2,5691893e-8 > 1e-8. Les réessais externes conservent leur rollback.
Étudier les compensations de flux, la normalisation dimensionnelle contrôlée,
les contraintes des parois et la cohérence des états avant tout changement.

## 4. Fonctions de paroi : contractuel mais NON intégré

Le tutoriel Driver–Seegmiller OpenFOAM utilise `omegaWallFunction`,
`nutUSpaldingWallFunction`, `kqRWallFunction`. V13/V14 possèdent un audit
indépendant du triplet omega/nu_t/G et des coins multi-face ; la valeur
SST résolue à la face `6nu/(beta1 y²)` de CFDX n'est **pas** une
implémentation équivalente. Le solveur emploie encore une diffusion de paroi
moléculaire pour le momentum, sans calcul simultané OpenFOAM de nu_t et G.

La prochaine étape doit être une implémentation *atomique et réversible* des
contraintes de paroi sur omega, nu_t et G, dans l'assemblage FV du transport
et du momentum. Ne pas remplacer un seul scalaire par sa référence en
ignorant les autres. Ensuite seulement générer le maillage `blockMesh`
OpenFOAM natif et qualifier x_r/h, Cf, U et les contraintes de Reynolds.

## 5. Preuves et reproduction

`test_pcd_bilateral_deflation_v14` rejoue une matrice 2x2 non symétrique et
une capture BFS réelle ; `test_global_nonlinear_contraction_v14` teste rejet,
rollback, borne du mérite. `qualify_pcd_bilateral_true_v14.py` revérifie
la cohérence des trois résidus sans promouvoir le COUPLED physique.
Suite ciblée : 6/6 CTest, 119 tests Python passés, 6 ignorés. Voir `preuves/` pour journaux exacts,
rejets et empreintes. Exécutable Release diagnostique : `main.cpp` compilé
avec -O0 uniquement pour accélérer le build ; aucune comparaison de vitesse
Release optimisée ni GPU/MPI réalisée. Aucune CI générale, aucune convergence
BFS longue revendiquée, aucun push GitHub.

## Sources méthodologiques
- PETSc KSP: https://petsc.org/release/manual/ksp/
- OpenFOAM Driver–Seegmiller: https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/
- OpenFOAM nutUSpaldingWallFunction: https://doc.openfoam.com/2312/tools/processing/boundary-conditions/rtm/derived/wall/nutUSpaldingWallFunction/
