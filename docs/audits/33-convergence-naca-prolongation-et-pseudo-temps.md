# FICHE MÃTHODE 33 â Convergence NACA, prolongement imbriquÃ© et contrÃ´le du pseudo-temps

**En-tÃªte d'Ã©vidence**
- PÃ©rimÃ¨tre : NACA0012 laminaire, solveur SIMPLE/SIMPLEC, initialisation d'un maillage raffinÃ© Ã  partir d'un Ã©tat moins fin, pseudo-transitoire optionnel et certificats de rÃ©sidus physiques.
- Statut : **correctif implÃ©mentÃ© et tests ciblÃ©s ; N64 UPWIND convergÃ© ; qualification TVD, N128, N256 et GCI encore ouverte**. Les expÃ©riences courtes sont explicitement non qualifiantes.
- Base : archives CFDX fournies le 2026-11-07, modifications locales sans push sur les PR parallÃ¨les.
- Preuves : `test_naca_prolong_initial_state.py`, `test_naca_convergence_diagnosis.py`, `test_sst_blending`, `test_steady_incompressible_solver`, historiques NACA locaux, tableaux A/B Ã  budget identique.

## 1. RÃ©fÃ©rences sources
- Ferziger, PeriÄ, Street, *Computational Methods for Fluid Dynamics*, Springer, DOI 10.1007/978-3-319-99693-6 (convergence non linÃ©aire, transfert entre grilles et volumes finis).
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html (ne pas confondre convergence itÃ©rative et spatiale).
- W. L. Oberkampf, C. J. Roy, *Verification and Validation in Scientific Computing*, DOI 10.1017/CBO9780511760396 (sÃ©paration vÃ©rification/validation).

## 2. Formulation continue
Pour un rÃ©gime permanent incompressible, on rÃ©sout la conservation de masse et de quantitÃ© de mouvement :
\[
 \nabla\cdot\mathbf{u}=0,\qquad
 \nabla\cdot(\rho\mathbf{u}\otimes\mathbf{u}) + \nabla p
 -\nabla\cdot[\mu(\nabla\mathbf{u}+\nabla\mathbf{u}^T)] -\rho\mathbf{f}=0.
\]
La pseudo-dÃ©rivÃ©e \(\partial(\rho\mathbf u)/\partial \tau\) n'est qu'une rÃ©gularisation algorithmique. Elle doit s'annuler au point fixe et ne doit jamais Ãªtre incluse dans le rÃ©sidu physique certifiÃ©.

## 3. DiscrÃ©tisation Volume Fini
Une Ã©tape de pseudo-transitoire ajoute \(M_P (\mathbf u_P^{k+1}-\mathbf u_P^k)\) au prÃ©dicteur. Les essais N64 ont identifiÃ© un CFL adaptatif trop petit : \(M_P\) devenait trÃ¨s supÃ©rieur Ã  la diagonale physique du momentum, immobilisant le champ. Une nouvelle limite **optionnelle** est appliquÃ©e Ã  la masse artificielle **avant** son insertion :
\[
 M_P\leftarrow\min\left(M_P,\;\gamma\min_i|A_{P,i}^{\rm physique}|\right),
 \quad \gamma>0.
\]
La variante testÃ©e utilise \(\gamma=0{,}25\), sans modification des tolÃ©rances physiques ni de la solution stationnaire. Le comportement historique reste conservÃ© si la borne est absente. Le pseudo-transitoire ne doit **pas** Ãªtre activÃ© par dÃ©faut sur NACA tant qu'il n'accÃ©lÃ¨re pas la convergence complÃ¨te.

Le prolongement N64 vers N128 utilise la topologie annulaire des grilles O emboÃ®tÃ©es (2n+1 colonnes et doublement des anneaux radiaux). Quatre cellules physiquement voisines dans les mÃªmes anneaux grossiers servent de support Ã  une reconstruction affine en coordonnÃ©es \((x,y)\), dont l'intercept fournit la valeur interpolÃ©e au centre fin. Les systÃ¨mes mal conditionnÃ©s, gÃ©omÃ©tries non apparentÃ©es et valeurs non finies sont rejetÃ©s. Le prolongement fournit un **Ã©tat initial non certifiÃ©**, pas une solution sur le maillage fin.

## 4. Traitement des frontiÃ¨res
La fermeture du bord de fuite n'est pas subdivisÃ©e : \(n_{\theta,f}=2n_{\theta,c}-1\), et le dernier intervalle constitue la fermeture de la couture. Les voisins de l'interpolation sont choisis selon l'adjacence **topologique** et non simplement selon la distance euclidienne, afin d'Ã©viter de mÃ©langer les deux faces de la paroi mince. Les champs \(u_x,u_y,p\) sont interpolÃ©s dans le plan, sans modifier les conditions aux limites du problÃ¨me raffinÃ©. La troisiÃ¨me coordonnÃ©e reste celle du maillage fin et doit Ãªtre la mÃªme extrusion plane.

## 5. PiÃ¨ges connus
- Une convergence de matrice de pression \(<10^{-10}\) ne valide pas le rÃ©sidu non linÃ©aire du momentum.
- Le rÃ©sidu global du momentum peut Ãªtre petit alors que le rÃ©sidu relatif de la composante \(y\), dont la norme source est petite, reste significatif. Le suivi des **trois** composantes est obligatoire.
- Un Ã©tat initial prolongÃ©, mÃªme issu d'un maillage convergÃ©, n'autorise pas l'Ã©tude GCI avant convergence complÃ¨te et indÃ©pendante des trois maillages.
- Une initialisation par plus proches voisins peut sÃ©lectionner des cellules en vis-Ã -vis d'une paroi mince ; utiliser l'adjacence structurÃ©e lorsque la topologie est vÃ©rifiÃ©e.
- Diminuer le CFL pseudo-transitoire peut **ralentir** trÃ¨s fortement la convergence lorsque l'accumulation fictive masque les changements physiques.
- La croissance du rÃ©sidu entre deux fractions de convection diffÃ©rentes n'est **pas** la croissance sur un mÃªme opÃ©rateur. Une protection fondÃ©e sur cette comparaison rejette artificiellement le dÃ©but de la rampe TVD. En revanche, les contrÃ´les de valeurs non finies restent toujours actifs et la protection de croissance redevient active lorsque la fraction ne change plus.
- `high_order_fraction=1` est compatible avec UPWIND : la validation d'ordre deux requiert le schÃ©ma TVD/Van Leer effectivement activÃ© et convergÃ©.

## 6. Ordre de prÃ©cision et comparaison
Les comparaisons A/B de prolongement sont faites Ã  paramÃ¨tres et budget identiques sur le mÃªme maillage cible. Sur N128, 24 itÃ©rations Ã  partir de l'Ã©tat N64 itÃ©ration 300 donnent pour momentum RHS \(5{,}09\times10^{-4}\) (IDW) contre \(4{,}22\times10^{-4}\) (interpolation structurÃ©e), et dÃ©faut local en faces \(y\) \(8{,}64\times10^{-3}\) contre \(5{,}59\times10^{-3}\). Cela mesure une **amÃ©lioration du dÃ©marrage**, pas une convergence N128 prouvÃ©e. Pour le pseudo-temps, les mÃªmes 80 itÃ©rations N64 donnent environ \(5{,}59\times10^{-3}\) (pseudotemps non bornÃ©) et \(3{,}81\times10^{-4}\) (fraction 0,25), mais le SIMPLE sans pseudo est meilleur : \(1{,}70\times10^{-5}\). Le pseudo-transitoire bornÃ© reste donc expÃ©rimental et dÃ©sactivÃ© par dÃ©faut.

L'acceptation finale exige convergence de tous les rÃ©sidus physiques et des efforts sur N64/N128/N256 avec la mÃªme convection bornÃ©e d'ordre supÃ©rieur, indÃ©pendance au maillage, GCI de pression/viscositÃ© cohÃ©rent, puis comparaison avec des donnÃ©es de rÃ©fÃ©rence Ã  incertitude connue.

## 7. Statut d'implÃ©mentation et preuves locales
- N64, UPWIND, \(Re=1000\) : convergÃ© en 336 itÃ©rations acceptÃ©es sans relÃ¢chement des seuils, mais **non qualifiÃ© en ordre deux**.
- Reconstruction initiale affine sur les centres des vraies grilles N64/N128 : erreur \(<7\times10^{-14}\) sur un champ affine de test. L'export Ã©crit les fichiers de faÃ§on atomique et refuse d'Ã©craser les fichiers sources.
- Sur N128, les rÃ©sidus baissent mieux avec le prolongement structurÃ© que l'IDW ancien aprÃ¨s 24 Ã©tapes, mais restent non convergÃ©s.
- Nouveau diagnostic CSV de toutes les composantes du momentum, des corrections et de la pression. Le programme d'analyse refuse explicitement toute promotion d'un CSV en certificat du solveur.
- La borne pseudo-transitoire optionnelle est couverte par un test du point fixe physique, de la compatibilitÃ© avec l'ancien comportement et des paramÃ¨tres invalides.
- La continuation TVD depuis N64 convergÃ© en UPWIND Ã©chouait dÃ¨s la transition UPWINDâTVD (itÃ©ration 41), le contrÃ´leur comparant Ã  tort les rÃ©sidus de deux opÃ©rateurs diffÃ©rents. Le nouveau prÃ©dicat `convection_residuals_share_operator` interdit cette comparaison Ã  une transition, et la rÃ©tablit aprÃ¨s la rampe. RÃ©gression `test_convection_continuation` verte. AprÃ¨s correction, la relance N64 atteint **100 % de TVD Ã  lâitÃ©ration 120**, puis termine **non convergÃ©e Ã  lâitÃ©ration 330** : rÃ©sidu momentum RHS `7.25783e-5`, maximum par composante `5.374e-3`, face momentum y `2.06470e-3`, continuitÃ© `2.09097e-7`. La pression linÃ©aire atteint `8.43316e-12` ; ce nâest pas un certificat non linÃ©aire. Le blocage dâentrÃ©e de rampe est corrigÃ©, mais pas la stagnation physique finale.
- RÃ©glage N128 expÃ©rimental Ã  `alpha_u=0.8` contre `0.6` : aprÃ¨s 36 itÃ©rations, meilleur rÃ©sidu global de momentum (`1.51163e-4`), mais dÃ©faut de continuitÃ© et rÃ©sidu indÃ©pendant en face y plus Ã©levÃ©s quâavec `alpha_u=0.6` (`1.39525e-5` et `5.72028e-3` respectivement). Ne pas promouvoir ce rÃ©glage par dÃ©faut.
- Un prolongement N128âN256 de dimensions vÃ©rifiÃ©es existe pour **initialisation expÃ©rimentale seulement**, Ã  partir dâun Ã©tat N128 **non convergÃ©**. Aucun certificat N256 nâa Ã©tÃ© produit.
- Aucune convergence TVD complÃ¨te N64/N128/N256 et aucune valeur GCI ne sont revendiquÃ©es par cette fiche.
