# FICHE MÃTHODE 32 â RÃ©sidu non linÃ©aire indÃ©pendant, GCI NACA, thermique/SST/rayonnement

**En-tÃªte d'Ã©vidence**
- PÃ©rimÃ¨tre : SIMPLE, SIMPLEC, COUPLED, audit volume-fini en faces sur l'Ã©tat acceptÃ©, famille NACA0012 64/128/256, contrÃ´le des efforts, tests Ã©lÃ©mentaires SST/Ã©nergie/rayonnement.
- Statut : **vÃ©rification logicielle ciblÃ©e ; qualification physique NACA/SST/rayonnement non obtenue**. Le calcul de GCI de production est interdit tant que les trois solutions ne convergent pas.
- Base : ZIP CFDX transmis le 2026-11-07 ; travaux sur une copie locale, sans push ni hypothÃ¨se d'Ã©quivalence avec les PR distantes.
- Tests : `test_conservation_boundedness`, `test_steady_incompressible_solver`, `test_naca0012_quick`, `test_naca_nested_radial_mesh.py`, `test_naca_force_mesh_qualification.py`, `test_qualify_naca_grid.py`, tests SST/thermiques/radiatifs.

## 1. RÃ©fÃ©rences sources
- PETSc SNES, *Nonlinear Solvers*, https://petsc.org/main/manual/snes/ : distinction entre rÃ©sidu non linÃ©aire et opÃ©rateur Picard/Jacobien gelÃ©.
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html : Ã©tudes de convergence spatiale et GCI.
- NASA Turbulence Modeling Resource, https://tmbwg.github.io/turbmodels/ : dÃ©finition des variantes SST et vÃ©rifications de rÃ©fÃ©rence Ã  Ã©tablir.
- G. de Vahl Davis, *Natural convection of air in a square cavity*, IJNMF 3 (1983), DOI 10.1002/fld.1650030305.
- F. Moukalled et al., *The Finite Volume Method in Computational Fluid Dynamics*, Springer 2016, DOI 10.1007/978-3-319-16874-6.

## 2. Formulation continue
Pour la quantitÃ© de mouvement incompressible, l'opÃ©rateur physique complet est
\[
 \frac{\partial (\rho u_i)}{\partial t}
 +\nabla\cdot (\rho\mathbf{u}u_i)
 -\nabla\cdot\left[\mu\left(\nabla u_i + (\nabla\mathbf u)^T_i\right)\right]
 +\frac{\partial p}{\partial x_i}-\rho f_i=0.
\]
Le rÃ©sidu algÃ©brique d'un systÃ¨me de Picard $A(U^k)U^{k+1}-b(U^k)$ n'est pas la rÃ©Ã©valuation de cet opÃ©rateur sur $(U^{k+1},p^{k+1})$.

## 3. DiscrÃ©tisation Volume Fini â Contrat de convergence
Pour chaque cellule $P$, le nouvel audit reconstruit les flux de convection et de viscositÃ© aux faces Ã  partir de l'Ã©tat acceptÃ©, soustrait les effets transpose si activÃ©s, ajoute le gradient de pression et la force volumique, puis ajoute l'accumulation physique Euler implicite ou BDF2 Ã  pas variable. Le pseudo-transitoire de stabilisation n'est **pas** une accumulation physique.

Le certificat conserve le maximum local $\max_P|R_{i,P}|$, au lieu du seul total signÃ© $\sum_PR_{i,P}$ oÃ¹ les erreurs de cellules voisines peuvent s'annuler. La normalisation emploie une Ã©chelle physique positive des flux, sources et accumulations. Les valeurs non finies, les normes nÃ©gatives et les Ã©tats de taille incorrecte sont rejetÃ©s. Ce certificat partage encore avec le solveur les flux de masse RhieâChow et certains opÃ©rateurs FV : il est indÃ©pendant **de la matrice gelÃ©e**, pas un oracle analytique universel.

`require_independent_face_residual_for_convergence` active la porte stricte au niveau du solveur, sur SIMPLE/SIMPLEC/COUPLED ; elle conserve la tolÃ©rance relative originale et laisse le comportement de production inchangÃ© si dÃ©sactivÃ©e. La campagne NACA l'active explicitement. La continuitÃ© et les rÃ©sidus de chaque composante doivent Ã©galement rester satisfaits.

## 4. Traitement des frontiÃ¨res et des maillages
Les flux aux frontiÃ¨res reprennent les lois de paroi, de symÃ©trie et les conditions de vitesse/pression rÃ©ellement dÃ©clarÃ©es ; les certificats de symÃ©trie et de jauge du lot 30 restent requis. Pour NACA, le raffinement radial interpolÃ© Ã  partir de 24 niveaux de base se fait sur les **coordonnÃ©es** : interpoler sÃ©parÃ©ment $r(\eta)$ et le poids de normale crÃ©erait un terme croisÃ© non linÃ©aire, faussant le rapport du premier espacement de paroi.

La famille candidate compte 3 096, 12 336 et 49 248 cellules ; les ratios mesurÃ©s de premiÃ¨re couche sont environ 2,027 et 2,014, avec zÃ©ro cellule inversÃ©e. Ce contrÃ´le est gÃ©omÃ©trique ; il ne prouve pas la convergence des coefficients aÃ©rodynamiques, ni l'entrÃ©e en rÃ©gime asymptotique.

## 5. PiÃ¨ges connus â Faux positifs dÃ©sormais bloquÃ©s
- Deux rÃ©sidus locaux Ã©gaux et opposÃ©s ne peuvent plus donner un Â« bilan parfaitement fermÃ© Â».
- Un rÃ©sidu de pression/AMG de $10^{-11}$ ne fait pas disparaÃ®tre un dÃ©faut de momentum ou de masse.
- Une valeur `high_order_fraction = 1` est insuffisante si le schÃ©ma est encore `UPWIND` : la qualification exige un schÃ©ma TVD bornÃ© et l'opÃ©rateur complet sur le dernier Ã©tat.
- Un $C_d$ monotone peut masquer l'absence de convergence de $C_{d,p}$ et $C_{d,v}$, compensant leurs erreurs. Les trois doivent prÃ©senter une signature de Richardson acceptable.
- Les grilles sans certificats d'efforts d'Ã©tats convergÃ©s produisent `BLOCKED`, jamais une valeur de GCI Â« estimÃ©e Â» Ã  partir des rÃ©sidus linÃ©aires.

## 6. Ordre de prÃ©cision et validation exigÃ©e
Avec des rapports de raffinement $r_{21}$ et $r_{32}$ lÃ©gÃ¨rement diffÃ©rents, l'ordre observÃ© $p$ est rÃ©solu par
\[
\frac{f(h_0)-f(h_1)}{f(h_1)-f(h_2)}
 =\frac{r_{21}^p-1}{1-r_{32}^{-p}},
\]
au lieu de remplacer les deux rapports par leur moyenne gÃ©omÃ©trique. Les GCI $C_d$, $C_{d,p}$ et $C_{d,v}$ sont candidats uniquement aprÃ¨s convergence de **toutes** les grilles et contrÃ´le gÃ©omÃ©trique. Un triplet seul ne dÃ©montre pas rigoureusement le rÃ©gime asymptotique : il faut confronter les profils $C_p/C_f$, une famille plus fine et les rÃ©fÃ©rences expÃ©rimentales indÃ©pendantes.

Tests SST, Ã©nergie, rayonnement et S2S vÃ©rifient des contrats/Ã©quations Ã©lÃ©mentaires. La validation NASA TMR du SST, les bilans CHT sur une Ã©tude de grille et les rÃ©fÃ©rences radiatives P1/S2S/DOM exigent encore des campagnes dÃ©diÃ©es ; Â« test vert Â» ne signifie pas Â« modÃ¨le validÃ© Â».

## 7. Statut d'implÃ©mentation et preuves locales
- Le certificat local et le contrÃ´le Euler/BDF2 ont Ã©tÃ© compilÃ©s avec rÃ©gressions C++ ; l'opt-in strict et les cas de zÃ©ro physique SIMPLE/SIMPLEC passent.
- Le GCI inÃ©gal, le maillage radial emboÃ®tÃ©, les contrÃ´les de schÃ©ma bornÃ©, la somme des efforts et les rÃ©sidus indÃ©pendants sont couverts par des tests nÃ©gatifs Python.
- Les diagnostics NACA de 120/18/6 itÃ©rations respectivement N64/N128/N256 restent **non convergÃ©s** ; aucune donnÃ©e de force obtenue Ã  ces budgets n'est certifiÃ©e. Le budget Ã©tait diagnostique, pas une preuve d'Ã©chec Ã  un budget intÃ©gral.
- Convection naturelle de Vahl Davis Ra=1000 : essais isolÃ©s sur 8 et 12 cellules par cÃ´tÃ© convergÃ©s avec bilans Ã©nergÃ©tiques de l'ordre de $10^{-12}$ ; l'essai 16Ã16 a atteint une limite de 200 s sans rÃ©sultat. Le binaire distingue dÃ©sormais explicitement `SINGLE_GRID_DIAGNOSTIC_PASS` de `PREQUALIFICATION_PASS` sur trois grilles ; preuve sur trois grilles distinctes et convergence vers les vitesses de rÃ©fÃ©rence encore ouvertes.
- Les tests Ã©lÃ©mentaires SST, thermique et radiatifs passent ; la promotion scientifique des modÃ¨les n'est pas autorisÃ©e Ã  ce stade.
