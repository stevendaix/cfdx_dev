# FICHE MÉTHODE 05 — Rhie–Chow axisymétrique VMFL036

**Statut :** en qualification — cette fiche documente la dérivation et les contrôles requis ; elle ne déclare pas la méthode physiquement qualifiée avant les campagnes sur le SHA final.

## 1. Domaine

Solveur axisymétrique natif (coordonnées méridiennes (x,r)) de `src/cfdx/physics/axisymmetric_vmfl036_solver.h`. Il est distinct des opérateurs Cartesian `Mesh/FvGeometry` et de la fiche 04.

## 2. Diagonale relaxée et réponse de pression

L'assemblage de quantité de mouvement stocke la diagonale relaxée
[
a_P^r = a_P/\alpha_u.
]
En gelant les termes hors diagonale et la source explicite de relaxation, la correction de vitesse associée à une correction de pression est
[
\mathbf{u}' = -\frac{V}{a_P^r}\nabla p'.
]
Comme `apx_` et `apr_` stockent (a_P^r), la réponse directionnelle est (V/apx_) et (V/apr_). Multiplier encore le dénominateur par `alpha_u` revient à annuler la relaxation une seconde fois.

La même réponse doit être utilisée dans :
- le terme de pression interpolé de Rhie–Chow ;
- le coefficient implicite de l'équation de correction de pression ;
- la correction conservative des flux de face ;
- la correction de vitesse par (\nabla p').

## 3. Décomposition géométrique sur une face oblique

Pour un vecteur unitaire (\mathbf e) du centre propriétaire au centre voisin et la normale unitaire (\mathbf n), l'opérateur implicite choisi est la réponse tensorielle
[
D_f^{\perp} = (\mathbf D_f\mathbf e)\cdot\mathbf n,
\qquad \mathbf D_f=\operatorname{diag}(V/a_{Px}^r,V/a_{Pr}^r)_f.
]
Le coefficient de pression associé est (\rho A_f D_f^{\perp}/d_{PN}). Le reste de (\mathbf D_f\nabla p\cdot\mathbf n) est différé explicitement. Le même (D_f^{\perp}) doit apparaître dans le flux de Rhie–Chow, la matrice de pression et la correction conservative ; la projection (\mathbf n^T\mathbf D_f\mathbf n) n'est équivalente que si (\mathbf e=\mathbf n) ou dans certains cas isotropes.

À la frontière de pression fixée, (\mathbf e_b) va du centre de cellule au centre de face ; le coefficient est (\rho A_f(\mathbf D_o\mathbf e_b)\cdot\mathbf n/d_{ob}), avec (p'_b=0).

## 4. Relaxation façon Majumdar

Le terme de correction de flux additionnel n'est pas activé en production par cette PR. Le helper de test reste disponible, mais sa combinaison avec le split RC résiduel doit être qualifiée sur un balayage (\alpha_u) avant activation. Il est interdit de conclure à la nécessité ou à l'absence de cette correction à partir de la seule notation (D^{RC}=\alpha_u D^{pc}) : les diagonales utilisées par chaque opérateur doivent être définies et dérivées en premier.

## 5. Gates de validation

- Champ de pression linéaire avec valeurs de frontière cohérentes ; contrôler la décomposition implicite + différée.
- Invariance sous ajout d'une constante à (p).
- Fermeture de continuité à partir des flux conservatifs autoritaires et antisymétrie interne.
- Balayage Re=100 sur le même maillage fine 48×72 : (\alpha_u=0.3,0.5,0.65,0.7,1.0), mêmes gates et budgets ; conserver (C_d), (C_{d,p}), (C_{d,v}), résidus et nombre d'itérations dans `validation-report/vmfl036-alpha-sweep.json`.
- Campagnes VMFL036 Re=50 et Re=100 coarse/medium/fine ; aucune modification de référence, bande de traînée ou tolérance n'est autorisée pour faire passer la PR.

## 6. Références

- Patankar, *Numerical Heat Transfer and Fluid Flow*, 1980 — relaxation et méthode SIMPLE.
- Rhie & Chow, « Numerical study of the turbulent flow past an airfoil with trailing edge separation », *AIAA Journal* 21(11), 1983 — interpolation momentum-weighted sur maillage collocatif.
- Majumdar, « Role of underrelaxation in momentum interpolation for calculation of flow with nonstaggered grids », *Numerical Heat Transfer* 13, 1988 — correction de flux sous relaxation.

## 7. État de preuve

La correction de code est en cours de validation. Cette fiche ne prétend pas que l'écart de (C_d) de +9,93 % est résolu. Le bilan de quantité de mouvement, la fermeture du flux et l'indépendance à (\alpha_u) doivent être vérifiés séparément avant toute conclusion physique.
