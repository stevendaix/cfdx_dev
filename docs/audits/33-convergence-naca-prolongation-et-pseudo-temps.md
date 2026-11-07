# FICHE MÉTHODE 33 — Convergence NACA, prolongement imbriqué et contrôle du pseudo-temps

**En-tête d'évidence**
- Périmètre : NACA0012 laminaire, solveur SIMPLE/SIMPLEC, initialisation d'un maillage raffiné à partir d'un état moins fin, pseudo-transitoire optionnel et certificats de résidus physiques.
- Statut : **correctif implémenté et tests ciblés ; N64 UPWIND convergé ; qualification TVD, N128, N256 et GCI encore ouverte**. Les expériences courtes sont explicitement non qualifiantes.
- Base : archives CFDX fournies le 2026-11-07, modifications locales sans push sur les PR parallèles.
- Preuves : `test_naca_prolong_initial_state.py`, `test_naca_convergence_diagnosis.py`, `test_sst_blending`, `test_steady_incompressible_solver`, historiques NACA locaux, tableaux A/B à budget identique.

## 1. Références sources
- Ferziger, Perić, Street, *Computational Methods for Fluid Dynamics*, Springer, DOI 10.1007/978-3-319-99693-6 (convergence non linéaire, transfert entre grilles et volumes finis).
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html (ne pas confondre convergence itérative et spatiale).
- W. L. Oberkampf, C. J. Roy, *Verification and Validation in Scientific Computing*, DOI 10.1017/CBO9780511760396 (séparation vérification/validation).

## 2. Formulation continue
Pour un régime permanent incompressible, on résout la conservation de masse et de quantité de mouvement :
\[
 \nabla\cdot\mathbf{u}=0,\qquad
 \nabla\cdot(\rho\mathbf{u}\otimes\mathbf{u}) + \nabla p
 -\nabla\cdot[\mu(\nabla\mathbf{u}+\nabla\mathbf{u}^T)] -\rho\mathbf{f}=0.
\]
La pseudo-dérivée \(\partial(\rho\mathbf u)/\partial \tau\) n'est qu'une régularisation algorithmique. Elle doit s'annuler au point fixe et ne doit jamais être incluse dans le résidu physique certifié.

## 3. Discrétisation Volume Fini
Une étape de pseudo-transitoire ajoute \(M_P (\mathbf u_P^{k+1}-\mathbf u_P^k)\) au prédicteur. Les essais N64 ont identifié un CFL adaptatif trop petit : \(M_P\) devenait très supérieur à la diagonale physique du momentum, immobilisant le champ. Une nouvelle limite **optionnelle** est appliquée à la masse artificielle **avant** son insertion :
\[
 M_P\leftarrow\min\left(M_P,\;\gamma\min_i|A_{P,i}^{\rm physique}|\right),
 \quad \gamma>0.
\]
La variante testée utilise \(\gamma=0{,}25\), sans modification des tolérances physiques ni de la solution stationnaire. Le comportement historique reste conservé si la borne est absente. Le pseudo-transitoire ne doit **pas** être activé par défaut sur NACA tant qu'il n'accélère pas la convergence complète.

Le prolongement N64 vers N128 utilise la topologie annulaire des grilles O emboîtées (2n+1 colonnes et doublement des anneaux radiaux). Quatre cellules physiquement voisines dans les mêmes anneaux grossiers servent de support à une reconstruction affine en coordonnées \((x,y)\), dont l'intercept fournit la valeur interpolée au centre fin. Les systèmes mal conditionnés, géométries non apparentées et valeurs non finies sont rejetés. Le prolongement fournit un **état initial non certifié**, pas une solution sur le maillage fin.

## 4. Traitement des frontières
La fermeture du bord de fuite n'est pas subdivisée : \(n_{\theta,f}=2n_{\theta,c}-1\), et le dernier intervalle constitue la fermeture de la couture. Les voisins de l'interpolation sont choisis selon l'adjacence **topologique** et non simplement selon la distance euclidienne, afin d'éviter de mélanger les deux faces de la paroi mince. Les champs \(u_x,u_y,p\) sont interpolés dans le plan, sans modifier les conditions aux limites du problème raffiné. La troisième coordonnée reste celle du maillage fin et doit être la même extrusion plane.

## 5. Pièges connus
- Une convergence de matrice de pression \(<10^{-10}\) ne valide pas le résidu non linéaire du momentum.
- Le résidu global du momentum peut être petit alors que le résidu relatif de la composante \(y\), dont la norme source est petite, reste significatif. Le suivi des **trois** composantes est obligatoire.
- Un état initial prolongé, même issu d'un maillage convergé, n'autorise pas l'étude GCI avant convergence complète et indépendante des trois maillages.
- Une initialisation par plus proches voisins peut sélectionner des cellules en vis-à-vis d'une paroi mince ; utiliser l'adjacence structurée lorsque la topologie est vérifiée.
- Diminuer le CFL pseudo-transitoire peut **ralentir** très fortement la convergence lorsque l'accumulation fictive masque les changements physiques.
- La croissance du résidu entre deux fractions de convection différentes n'est **pas** la croissance sur un même opérateur. Une protection fondée sur cette comparaison rejette artificiellement le début de la rampe TVD. En revanche, les contrôles de valeurs non finies restent toujours actifs et la protection de croissance redevient active lorsque la fraction ne change plus.
- `high_order_fraction=1` est compatible avec UPWIND : la validation d'ordre deux requiert le schéma TVD/Van Leer effectivement activé et convergé.

## 6. Ordre de précision et comparaison
Les comparaisons A/B de prolongement sont faites à paramètres et budget identiques sur le même maillage cible. Sur N128, 24 itérations à partir de l'état N64 itération 300 donnent pour momentum RHS \(5{,}09\times10^{-4}\) (IDW) contre \(4{,}22\times10^{-4}\) (interpolation structurée), et défaut local en faces \(y\) \(8{,}64\times10^{-3}\) contre \(5{,}59\times10^{-3}\). Cela mesure une **amélioration du démarrage**, pas une convergence N128 prouvée. Pour le pseudo-temps, les mêmes 80 itérations N64 donnent environ \(5{,}59\times10^{-3}\) (pseudotemps non borné) et \(3{,}81\times10^{-4}\) (fraction 0,25), mais le SIMPLE sans pseudo est meilleur : \(1{,}70\times10^{-5}\). Le pseudo-transitoire borné reste donc expérimental et désactivé par défaut.

L'acceptation finale exige convergence de tous les résidus physiques et des efforts sur N64/N128/N256 avec la même convection bornée d'ordre supérieur, indépendance au maillage, GCI de pression/viscosité cohérent, puis comparaison avec des données de référence à incertitude connue.

## 7. Statut d'implémentation et preuves locales
- N64, UPWIND, \(Re=1000\) : convergé en 336 itérations acceptées sans relâchement des seuils, mais **non qualifié en ordre deux**.
- Reconstruction initiale affine sur les centres des vraies grilles N64/N128 : erreur \(<7\times10^{-14}\) sur un champ affine de test. L'export écrit les fichiers de façon atomique et refuse d'écraser les fichiers sources.
- Sur N128, les résidus baissent mieux avec le prolongement structuré que l'IDW ancien après 24 étapes, mais restent non convergés.
- Nouveau diagnostic CSV de toutes les composantes du momentum, des corrections et de la pression. Le programme d'analyse refuse explicitement toute promotion d'un CSV en certificat du solveur.
- La borne pseudo-transitoire optionnelle est couverte par un test du point fixe physique, de la compatibilité avec l'ancien comportement et des paramètres invalides.
- La continuation TVD depuis N64 convergé en UPWIND échouait dès la transition UPWIND→TVD (itération 41), le contrôleur comparant à tort les résidus de deux opérateurs différents. Le nouveau prédicat `convection_residuals_share_operator` interdit cette comparaison à une transition, et la rétablit après la rampe. Régression `test_convection_continuation` verte. Après correction, la relance N64 atteint **100 % de TVD à l’itération 120**, puis termine **non convergée à l’itération 330** : résidu momentum RHS `7.25783e-5`, maximum par composante `5.374e-3`, face momentum y `2.06470e-3`, continuité `2.09097e-7`. La pression linéaire atteint `8.43316e-12` ; ce n’est pas un certificat non linéaire. Le blocage d’entrée de rampe est corrigé, mais pas la stagnation physique finale.
- Réglage N128 expérimental à `alpha_u=0.8` contre `0.6` : après 36 itérations, meilleur résidu global de momentum (`1.51163e-4`), mais défaut de continuité et résidu indépendant en face y plus élevés qu’avec `alpha_u=0.6` (`1.39525e-5` et `5.72028e-3` respectivement). Ne pas promouvoir ce réglage par défaut.
- Un prolongement N128→N256 de dimensions vérifiées existe pour **initialisation expérimentale seulement**, à partir d’un état N128 **non convergé**. Aucun certificat N256 n’a été produit.
- Aucune convergence TVD complète N64/N128/N256 et aucune valeur GCI ne sont revendiquées par cette fiche.
