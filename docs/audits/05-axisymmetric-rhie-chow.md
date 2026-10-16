# FICHE MÉTHODE 05 — Rhie–Chow axisymétrique VMFL036

**En-tête d'évidence**
- Scope : solveur axisymétrique natif (coordonnées méridiennes x,r), pression-vitesse, flux conservatif et traînée VMFL036.
- Statut : EN QUALIFICATION — la dérivation est documentée, mais l'indépendance à alpha_u et l'impact physique restent soumis à la campagne sur le SHA exact.
- Issues/PRs liées : #915 (preuve VMFL036), #919 (réponse pression relaxée et opérateur sur faces obliques).
- Sources : `src/cfdx/physics/axisymmetric_vmfl036_solver.h`, `tests/validation/test_vmfl036_axisymmetric.cpp`.

## 1. Références sources
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, 1980, ISBN: 978-0-89116-522-4 — relaxation et méthode SIMPLE.
- C. M. Rhie, W. L. Chow, « Numerical study of the turbulent flow past an airfoil with trailing edge separation », *AIAA Journal* 21(11), 1983.
- B. Majumdar, « Role of underrelaxation in momentum interpolation for calculation of flow with nonstaggered grids », *Numerical Heat Transfer* 13, 1988.
- Contrat de preuve local : `docs/validation/VMFL036_FINE_GAP_STUDY.md`.

## 2. Formulation continue
Pour l'écoulement incompressible axisymétrique, la conservation de masse s'écrit comme une divergence du flux de révolution sur les faces méridiennes. La traînée totale est la somme des contributions de pression et visqueuse :
[
C_D=C_{D,p}+C_{D,v}.
]
Les trois quantités doivent être enregistrées séparément : deux erreurs de signe opposé peuvent se compenser dans (C_D) et donner un total trompeusement proche de la référence.

## 3. Discrétisation Volume Fini
### Diagonale relaxée
L'assemblage de quantité de mouvement stocke la diagonale relaxée
[
a_P^r=a_P/\alpha_u.
]
En gelant les termes hors diagonale et la source explicite de relaxation, la correction de vitesse associée à une correction de pression est
[
\mathbf u'=-\frac{V}{a_P^r}\nabla p'.
]
Comme `apx_` et `apr_` stockent déjà (a_P^r), la réponse est (V/apx_) et (V/apr_). Multiplier encore le dénominateur par `alpha_u` annule la relaxation une deuxième fois.

Cette réponse doit être commune au flux Rhie–Chow, à la matrice de correction de pression, à la correction conservative du flux et à la correction de vitesse.

### Face oblique
Pour le vecteur unitaire (\mathbf e) du centre propriétaire au centre voisin, la normale unitaire (\mathbf n) et
[
\mathbf D_f=\operatorname{diag}\left((V/a_{Px}^r)_f,(V/a_{Pr}^r)_f\right),
]
le coefficient implicite choisi est
[
D_f^{\perp}=(\mathbf D_f\mathbf e)\cdot\mathbf n.
]
Le coefficient matriciel est (\rho A_f D_f^{\perp}/d_{PN}). Le reste du flux tensoriel ( (\mathbf D\nabla p)_f\cdot\mathbf S_f-D_f^{\perp}(p_N-p_P)/d_{PN}) est différé explicitement. Le même coefficient doit être employé dans la matrice, le flux de Rhie–Chow et la correction conservative. La projection (\mathbf n^T\mathbf D_f\mathbf n) n'est équivalente que dans des géométries particulières (face orthogonale ou réponse isotrope).

## 4. Traitement des frontières
- À la sortie à pression fixée, (p'_b=0). Le coefficient emploie (\mathbf e_b) du centre de cellule au centre de face : (\rho A_f(\mathbf D_o\mathbf e_b)\cdot\mathbf n/d_{ob}).
- Les mêmes réponse et distance centre-face sont réutilisées lors de la correction conservative du flux.
- Les faces d'axe ont une aire de révolution nulle ; elles ne doivent pas injecter de flux massique.
- La fermeture de continuité est calculée avec les flux conservatifs autoritaires, avec orientation propriétaire/voisin opposée sur les deux cellules.

## 5. Pièges connus
- Diagonale relaxée et non relaxée mélangées entre deux opérateurs : facteur parasite (\alpha_u) ou (1/\alpha_u).
- Projection uniquement sur la normale de face sur maillage oblique : le coefficient peut ne plus être celui de la décomposition centre-à-centre.
- Remplacer l'interpolation de (\mathbf D\nabla p) par (\operatorname{interp}(\mathbf D)\operatorname{interp}(\nabla p)) lorsque (\mathbf D) varie entre cellules : le flux reconstruit peut différer de celui qui a servi à former le RHS de pression.
- Activer une correction de Majumdar sans la dériver pour le split RC effectivement utilisé : le helper reste testable, mais cette correction additionnelle n'est pas activée dans la boucle de production par la PR #919.
- (C_D) total seul n'est pas une preuve indépendante : (C_{D,p}), (C_{D,v}), la fermeture de quantité de mouvement et la continuité doivent aussi être vérifiés.

## 6. Ordre de précision attendu et validation exigée
La correction d'opérateur ne doit pas modifier les références, gates, tolérances ou bandes de traînée pour rendre les tests verts.

La preuve requise comprend :
1. décomposition implicite + différée sur un champ linéaire avec conditions de frontière compatibles ;
2. invariance à un décalage uniforme de pression ;
3. fermeture de continuité indépendante sur les flux conservatifs ;
4. balayage Re=100 sur le même maillage fine 48×72 pour (\alpha_u=0.3,0.5,0.65,0.7,1.0), budgets inchangés ;
5. campagnes VMFL036 Re=50 et Re=100 coarse/medium/fine ;
6. étude de (C_{D,p}) et (C_{D,v}) en plus de (C_D), puis raffinement radial indépendant du raffinement tangentiel.

Le test long `test_vmfl036_alpha_sensitivity` produit `validation-report/vmfl036-alpha-sweep.json` avec le SHA, paramètres, résidus, itérations, (C_{D,p}), (C_{D,v}), (C_D) et dispersions relatives. Il échoue si un cas ne converge pas ou produit une traînée non finie. Aucun seuil d'indépendance à (\alpha_u) n'est affirmé avant examen de la première preuve.

## 7. Statut d'implémentation CFDX
La PR #919 aligne la réponse issue de la diagonale relaxée dans le flux RC, la matrice de pression, la correction du flux et la correction de vitesse. Elle remplace aussi la projection normale seule par ( (\mathbf D_f\mathbf e)\cdot\mathbf n) sur les faces obliques et apparie le terme différé au même saut de pression que la matrice.

**Statut de preuve : en attente.** Le changement est une modification du solveur physique ; il ne démontre pas encore que l'écart (C_D=+9,93\%) est résolu. La campagne alpha-sweep et les validations physiques sur le SHA final sont bloquantes. La bande actuelle de 10 % et la référence (C_D=1,0895) restent inchangées.

## 8. Journal
- 2026-10-16 — dérivation de la réponse relaxée et du split tensoriel ; ajout de la preuve alpha-sweep en PR #919.
