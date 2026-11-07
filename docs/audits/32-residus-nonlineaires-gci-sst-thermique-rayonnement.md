# FICHE MÉTHODE 32 — Résidu non linéaire indépendant, GCI NACA, thermique/SST/rayonnement

**En-tête d'évidence**
- Périmètre : SIMPLE, SIMPLEC, COUPLED, audit volume-fini en faces sur l'état accepté, famille NACA0012 64/128/256, contrôle des efforts, tests élémentaires SST/énergie/rayonnement.
- Statut : **vérification logicielle ciblée ; qualification physique NACA/SST/rayonnement non obtenue**. Le calcul de GCI de production est interdit tant que les trois solutions ne convergent pas.
- Base : ZIP CFDX transmis le 2026-11-07 ; travaux sur une copie locale, sans push ni hypothèse d'équivalence avec les PR distantes.
- Tests : `test_conservation_boundedness`, `test_steady_incompressible_solver`, `test_naca0012_quick`, `test_naca_nested_radial_mesh.py`, `test_naca_force_mesh_qualification.py`, `test_qualify_naca_grid.py`, tests SST/thermiques/radiatifs.

## 1. Références sources
- PETSc SNES, *Nonlinear Solvers*, https://petsc.org/main/manual/snes/ : distinction entre résidu non linéaire et opérateur Picard/Jacobien gelé.
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html : études de convergence spatiale et GCI.
- NASA Turbulence Modeling Resource, https://tmbwg.github.io/turbmodels/ : définition des variantes SST et vérifications de référence à établir.
- G. de Vahl Davis, *Natural convection of air in a square cavity*, IJNMF 3 (1983), DOI 10.1002/fld.1650030305.
- F. Moukalled et al., *The Finite Volume Method in Computational Fluid Dynamics*, Springer 2016, DOI 10.1007/978-3-319-16874-6.

## 2. Formulation continue
Pour la quantité de mouvement incompressible, l'opérateur physique complet est
\[
 \frac{\partial (\rho u_i)}{\partial t}
 +\nabla\cdot (\rho\mathbf{u}u_i)
 -\nabla\cdot\left[\mu\left(\nabla u_i + (\nabla\mathbf u)^T_i\right)\right]
 +\frac{\partial p}{\partial x_i}-\rho f_i=0.
\]
Le résidu algébrique d'un système de Picard $A(U^k)U^{k+1}-b(U^k)$ n'est pas la réévaluation de cet opérateur sur $(U^{k+1},p^{k+1})$.

## 3. Discrétisation Volume Fini — Contrat de convergence
Pour chaque cellule $P$, le nouvel audit reconstruit les flux de convection et de viscosité aux faces à partir de l'état accepté, soustrait les effets transpose si activés, ajoute le gradient de pression et la force volumique, puis ajoute l'accumulation physique Euler implicite ou BDF2 à pas variable. Le pseudo-transitoire de stabilisation n'est **pas** une accumulation physique.

Le certificat conserve le maximum local $\max_P|R_{i,P}|$, au lieu du seul total signé $\sum_PR_{i,P}$ où les erreurs de cellules voisines peuvent s'annuler. La normalisation emploie une échelle physique positive des flux, sources et accumulations. Les valeurs non finies, les normes négatives et les états de taille incorrecte sont rejetés. Ce certificat partage encore avec le solveur les flux de masse Rhie–Chow et certains opérateurs FV : il est indépendant **de la matrice gelée**, pas un oracle analytique universel.

`require_independent_face_residual_for_convergence` active la porte stricte au niveau du solveur, sur SIMPLE/SIMPLEC/COUPLED ; elle conserve la tolérance relative originale et laisse le comportement de production inchangé si désactivée. La campagne NACA l'active explicitement. La continuité et les résidus de chaque composante doivent également rester satisfaits.

## 4. Traitement des frontières et des maillages
Les flux aux frontières reprennent les lois de paroi, de symétrie et les conditions de vitesse/pression réellement déclarées ; les certificats de symétrie et de jauge du lot 30 restent requis. Pour NACA, le raffinement radial interpolé à partir de 24 niveaux de base se fait sur les **coordonnées** : interpoler séparément $r(\eta)$ et le poids de normale créerait un terme croisé non linéaire, faussant le rapport du premier espacement de paroi.

La famille candidate compte 3 096, 12 336 et 49 248 cellules ; les ratios mesurés de première couche sont environ 2,027 et 2,014, avec zéro cellule inversée. Ce contrôle est géométrique ; il ne prouve pas la convergence des coefficients aérodynamiques, ni l'entrée en régime asymptotique.

## 5. Pièges connus — Faux positifs désormais bloqués
- Deux résidus locaux égaux et opposés ne peuvent plus donner un Â« bilan parfaitement fermé Â».
- Un résidu de pression/AMG de $10^{-11}$ ne fait pas disparaître un défaut de momentum ou de masse.
- Une valeur `high_order_fraction = 1` est insuffisante si le schéma est encore `UPWIND` : la qualification exige un schéma TVD borné et l'opérateur complet sur le dernier état.
- Un $C_d$ monotone peut masquer l'absence de convergence de $C_{d,p}$ et $C_{d,v}$, compensant leurs erreurs. Les trois doivent présenter une signature de Richardson acceptable.
- Les grilles sans certificats d'efforts d'états convergés produisent `BLOCKED`, jamais une valeur de GCI Â« estimée Â» à partir des résidus linéaires.

## 6. Ordre de précision et validation exigée
Avec des rapports de raffinement $r_{21}$ et $r_{32}$ légèrement différents, l'ordre observé $p$ est résolu par
\[
\frac{f(h_0)-f(h_1)}{f(h_1)-f(h_2)}
 =\frac{r_{21}^p-1}{1-r_{32}^{-p}},
\]
au lieu de remplacer les deux rapports par leur moyenne géométrique. Les GCI $C_d$, $C_{d,p}$ et $C_{d,v}$ sont candidats uniquement après convergence de **toutes** les grilles et contrôle géométrique. Un triplet seul ne démontre pas rigoureusement le régime asymptotique : il faut confronter les profils $C_p/C_f$, une famille plus fine et les références expérimentales indépendantes.

Tests SST, énergie, rayonnement et S2S vérifient des contrats/équations élémentaires. La validation NASA TMR du SST, les bilans CHT sur une étude de grille et les références radiatives P1/S2S/DOM exigent encore des campagnes dédiées ; Â« test vert Â» ne signifie pas Â« modèle validé Â».

## 7. Statut d'implémentation et preuves locales
- Le certificat local et le contrôle Euler/BDF2 ont été compilés avec régressions C++ ; l'opt-in strict et les cas de zéro physique SIMPLE/SIMPLEC passent.
- Le GCI inégal, le maillage radial emboîté, les contrôles de schéma borné, la somme des efforts et les résidus indépendants sont couverts par des tests négatifs Python.
- Les diagnostics NACA de 120/18/6 itérations respectivement N64/N128/N256 restent **non convergés** ; aucune donnée de force obtenue à ces budgets n'est certifiée. Le budget était diagnostique, pas une preuve d'échec à un budget intégral.
- Convection naturelle de Vahl Davis Ra=1000 : essais isolés sur 8 et 12 cellules par côté convergés avec bilans énergétiques de l'ordre de $10^{-12}$ ; l'essai 16Ã16 a atteint une limite de 200 s sans résultat. Le binaire distingue désormais explicitement `SINGLE_GRID_DIAGNOSTIC_PASS` de `PREQUALIFICATION_PASS` sur trois grilles ; preuve sur trois grilles distinctes et convergence vers les vitesses de référence encore ouvertes.
- Les tests élémentaires SST, thermique et radiatifs passent ; la promotion scientifique des modèles n'est pas autorisée à ce stade.
