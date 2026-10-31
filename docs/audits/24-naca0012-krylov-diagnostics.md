# FICHE MÉTHODE 24 — Microscope de diagnostic NACA0012 / solveurs scalaires

**Statut : INSTRUMENTATION AJOUTÉE — CAUSE RACINE EN ATTENTE DES ARTEFACTS CI**

## 1. Références sources (obligatoire : DOI ou ISBN)

- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd edition, SIAM, 2003, ISBN 978-0-89871-534-7 — méthodes de Krylov, préconditionnement et résidus.
- R. S. Varga, *Matrix Iterative Analysis*, 2nd edition, Springer, 2000, ISBN 978-3-540-66311-8 — propriétés spectrales et analyse matricielle.

## 2. Formulation continue

CG résout le système linéaire (Ax=b) sous les hypothèses de symétrie et de définie-positivité de (A). Le préconditionneur doit également induire une action compatible avec la positivité. Le scalaire initial (r_0^T M^{-1}r_0) doit être fini et strictement positif. Une valeur non positive ou non finie interdit de former un pas CG ; elle constitue un diagnostic de compatibilité, pas une raison de modifier le seuil de convergence.

Le résidu certifié est (r=b-Ax), recalculé en FP64. Le résidu récursif et les produits préconditionnés sont des diagnostics distincts et ne remplacent jamais cette certification.

## 3. Discrétisation Volume Fini

Le prédicteur momentum produit un système creux (A_U U=b_U). Pour diagnostiquer une divergence après un seul pas, il faut distinguer :
- la qualité de l'assemblage (dimensions, nnz, coefficients non finis, diagonale) ;
- les propriétés de l'opérateur (norme infinie et défaut de symétrie) ;
- l'échelle du second membre et du résidu vrai ;
- le choix du solveur Krylov et le préconditionneur réellement demandé ;
- le comportement des solveurs de repli, exécutés ou non.

Le défaut de symétrie rapporté est la norme infinie absolue de (A-A^T), calculée sur les coefficients présents. Il s'agit d'un indicateur de diagnostic, pas d'un seuil d'acceptation ajouté.

## 4. Traitement des frontières

L'instrumentation est placée au niveau du solveur scalaire, après assemblage et après les tentatives de solve. Elle observe la matrice et le RHS effectivement transmis au solveur. Elle ne change aucune condition limite, aucun flux, aucune jauge de pression, aucune ligne de matrice, ni l'ordre des fallbacks.

Les traces se limitent aux scalaires de diagnostic, dimensions, compteurs et statuts.

## 5. Pièges connus

Campagne `CFDX total validation`, job `Long validation campaigns`, commit
`ca0180472e7b57e2f616f406d96d580124754693` (31 octobre 2026).

- `test_naca0012_qualification` rapporte une valeur provisoire
  `NACA0012_RE1000_A0 level=64` : Cd=0.0901496, Cl=-9.15e-11,
  Cd pression=0.0344606, Cd visqueux=0.055689, continuité=8.84e-11,
  456 itérations. Ce message ne constitue pas un PASS final.
- Une correction de pression échoue à l'itération non linéaire 20 :
  status=1, 457 itérations, résidu=3.021189409e2, relatif=1.181456e-9.
- À l'itération 22, le solve Ux échoue :
  status=3 (**NOT_APPLICABLE**, et non DIVERGED), itérations=1,
  résidu rapporté=1.329e25. L'énumération publique vaut
  CONVERGED=0, MAX_ITER_REACHED=1, DIVERGED=2, NOT_APPLICABLE=3.
  Le fallback CG retourne NOT_APPLICABLE lorsque la courbure Krylov
  `p^T A p` est non positive ou non finie. Avant ce correctif, cette branche
  ne conservait pas d'échantillon, d'où `samples=0` et l'impossibilité
  de distinguer ces deux cas dans le log. Les retries suivants rapportent
  des relatifs 5.482 puis 2.252 avant l'abandon.
- Le test de cisaillement rapporte provisoirement Cd_visqueux GG=0.055689,
  WLS=0.0580625, 1D=0.106308 (spread relatif=0.690069), puis échoue sur le
  même problème Ux. Ce spread doit être réévalué sur des solutions convergées
  et comparables ; il ne permet pas encore d'isoler l'erreur de gradient.

Ces traces établissent une défaillance linéaire/non linéaire, mais ne permettent
pas de conclure si l'origine est la matrice assemblée, le choix Krylov, le
préconditionneur, le RHS ou l'état initial.

## 6. Ordre de précision attendu et validation exigée

Sur tout échec final d'un solve scalaire, le journal doit émettre :
- `SCALAR_LINEAR_ATTEMPT` pour le solve demandé et chaque fallback réellement
  exécuté (GMRES, Gauss-Seidel, CG), avec statut, itérations, résidus et dernier
  échantillon Krylov disponible ;
- `SCALAR_LINEAR_FAILURE` avec dimensions, nnz, norme infinie de matrice,
  défaut de symétrie absolu, bornes diagonales, norme RHS, nombre de coefficients
  et RHS non finis, solveur/préconditionneur demandés et politique de fallback ;
- `SCALAR_LINEAR_KRYLOV_SAMPLE` avec les scalaires de breakdown disponibles.

Dans CG :
- toutes les sorties numériques `DIVERGED` après le calcul du résidu initial
  recalculent le résidu vrai FP64 et conservent les scalaires Krylov ;
- le premier produit préconditionné non positif conserve également un échantillon ;
- le statut public reste `DIVERGED` ou `NOT_APPLICABLE` selon le cas ;
- aucun nouvel enum, critère, tolérance ou budget n'est introduit.

Deux tests injectent respectivement un produit préconditionné initial négatif et
un produit fini dont le produit scalaire déborde. Ils vérifient le statut,
l'itération zéro, le résidu vrai, l'échantillon diagnostique et l'absence de
modification de la solution.

Après CI sur le HEAD exact, examiner dans cet ordre :
1. les solveurs demandés et les fallbacks réellement exécutés ;
2. les dimensions, diagonales, valeurs non finies et défaut de symétrie ;
3. le premier scalaire Krylov anormal (`rMz`, `pAp`, `alpha`) ;
4. l'écart entre résidu vrai et récursif, puis l'échelle du RHS ;
5. seulement ensuite, choisir une correction de formulation, de préconditionnement
   ou de sélection de Krylov.

Release, DebugSanitizers, tests numériques et campagnes NACA0012/Ghia doivent
être relancés sur le SHA exact. Un résultat algébrique ne vaut pas qualification
physique, convergence en maillage ou validation du cisaillement pariétal.

## 7. Statut d'implémentation

- Modifié : `src/cfdx/core/linalg/cg_solver.h` — résidu vrai FP64 sur breakdown
  numérique et échantillon des produits initiaux non finis/non positifs.
- Modifié : `src/cfdx/physics/finite_volume_transport.h` — traces des tentatives,
  statistiques de matrice et scalaires Krylov uniquement sur échec.
- Modifié : `tests/unit/test_cg_solver.cpp` — tests de régression des produits
  préconditionnés initiaux négatif et non fini.
- Aucune tolérance, aucun budget d'itérations, aucun fallback et aucune porte
  physique n'a été modifié.
- Cause racine NACA : **NON DÉTERMINÉE** tant que les nouveaux artefacts ne sont
  pas disponibles et examinés.
- Décision : **PENDING**. Ne pas fusionner avant CI verte sur le HEAD exact et
  examen des diagnostics NACA ciblés.


## 8. Addendum exact-HEAD — campagne CI du 2026-10-31 (HEAD parent `298e4d01`)

Le run `41184306862` échoue uniquement dans les deux campagnes NACA longues :
`test_naca0012_qualification` et `test_naca0012_wall_shear_audit`.
La campagne Ghia profils N=20/40/80/160 passe en 3436.25 s (N=160 :
8600 itérations, résidu momentum vrai `8.31e-11`). Les autres partitions
de validation-total, dont N9, thermique/rayonnement et wall-distance, passent.

Sur N128, les traces montrent d'abord un échec de correction de pression à
l'itération non linéaire 20, puis des échecs du prédicteur Ux à l'itération 22.
Le RHS momentum et la norme de matrice croissent entre les retries alors que
les solveurs demandés/GMRES/Gauss-Seidel n'atteignent pas la tolérance.
L'échec final ne doit donc pas être attribué au seul CG : le RHS est déjà
fortement dégradé avant le fallback CG. Les statistiques globales actuelles
ne donnent toutefois ni l'indice de la ligne maximale, ni les contributions
pression/convection/diffusion/source ; la cause d'assemblage reste ouverte.

**Correction de diagnostic incluse dans le HEAD suivant :** sur une courbure
`p^T A p <= 0` ou non finie, CG conserve désormais un échantillon avec
le résidu vrai FP64, `r^T M^{-1}r`, `p^T A p` et `alpha` non applicable.
Une régression exige cet échantillon sur le cas SPD-indéfini de test.
Aucun algorithme, seuil, budget, fallback ou gate physique n'est modifié.
