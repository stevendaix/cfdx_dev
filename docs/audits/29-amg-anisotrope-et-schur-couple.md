# FICHE MÉTHODE 29 — AMG anisotrope, interpolation énergétique et Schur couplé

**En-tête d'évidence**
- **Scope :** matrices de diffusion/pression anisotropes et système vitesse–pression collocated, préconditionneurs natifs CFDX.
- **Statut : DRAFT (2026-11-06).** À revoir et approuver avant promotion en production ; cette fiche ne certifie pas la PR #1182.
- **Registres :** `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json`.
- **Issue/PR :** #1182 (rapport de revue joint) ; la comparaison HEAD GitHub n'a pas pu être effectuée dans cette copie.
- **Base vérifiée :** archive locale `cfdx_corrige.zip` ; pas de vérification GitHub du HEAD `012840d1`.
- **Limites :** CPU sériel, preuve CFD N128/N256, SST, maillages 3D et MPI incomplète ; méthodes anisotropes opt-in.

## 1. Références sources

1. Mavriplis, D. J. (1998), *Multigrid Strategies for Viscous Flow Solvers on Anisotropic Unstructured Meshes*, JCP 145, 141–165, DOI **10.1006/jcph.1998.6036**. Couplage semi-coarsening et lissage en lignes, et non l'un seul.
2. Mavriplis (1998), *Directional Agglomeration Multigrid Techniques for High-Reynolds Number Viscous Flows*, NASA/CR-1998-206911, https://ntrs.nasa.gov/citations/19980041394 ; maillages 2D/3D non structurés et préconditionnement directionnel.
3. Brandt, Brannick, Kahl & Livshits (2015), *Algebraic Distance for Anisotropic Diffusion Problems: Multilevel Results*, ETNA 44, 472–496 ; détection non alignée par modes de test et distance algébrique.
4. Olson, Schroder & Tuminaro (2011), *A General Interpolation Strategy for Algebraic Multigrid Using Energy Minimization*, SIAM JSC 33, 966–991, DOI **10.1137/100803031** ; contrainte de support, conservation des modes et minimisation colonne par colonne. **La minimisation CFDX à pas borné n'est pas une réimplémentation équivalente de l'algorithme de l'article.**
5. Braess & Sarazin (1997), *An Efficient Smoother for the Stokes Problem*, Applied Numerical Mathematics 23, 3–19, DOI **10.1016/S0168-9274(96)00059-1**.
6. Vanka (1986), *Block-Implicit Multigrid Solution of Navier-Stokes Equations in Primitive Variables*, JCP 65, 138–158, DOI **10.1016/0021-9991(86)90008-2**. L'article original emploie des variables décalées ; la variante CFDX collocated exige des preuves propres.
7. Kay, Loghin & Wathen (2002), *A Preconditioner for the Steady-State Navier–Stokes Equations*, SIAM JSC 24, 237–256 ; PCD.
8. Elman et al. (2006), *Block Preconditioners Based on Approximate Commutators*, SIAM JSC, DOI **10.1137/040608817** ; LSC.
9. Elman et al. (2008), *A Taxonomy and Comparison of Parallel Block Multi-Level Preconditioners for the Incompressible Navier–Stokes Equations*, JCP 227, 1790–1808 ; MGR / bloc Schur.
10. PETSc PCFIELDSPLIT, https://petsc.org/main/manualpages/PC/PCFIELDSPLIT/ ; HYPRE BoomerAMG, https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html . Références de conception uniquement : CFDX n'emploie pas ces bibliothèques en exécution.

## 2. Formulation continue

Pour une diffusion éventuellement tensorielle, `-div(K grad(phi)) = f` avec `K` symétrique définie positive dans le domaine elliptique qualifié. La qualité spectrale de sa discrétisation dépend à la fois de `K`, de l'étirement et de la non-orthogonalité des cellules. Aucune correction d'AMG ne peut garantir la précision d'un flux géométriquement faux.

Pour l'écoulement incompressible, `rho (u · grad) u - div(2 mu_eff D(u)) + grad(p) = f` et `div(u)=0` (CL, turbulence et flux Rhie–Chow indépendamment qualifiés).

## 3. Discrétisation Volume Fini

- Matrice scalaire CSR `A x=b`, pression avec jauge ou projection de l'espace nul constant selon CL. La validation impose `||b-Ax||_2 / ||b||_2` **sur la matrice originale** ; la norme inf et les bilans de flux sont mesurés séparément.
- Pour un gradient de face, le flux physique s'écrit `-Gamma S_f · grad(phi)_f` : partie orthogonale implicitement assemblée et correction non orthogonale/skewness, explicite ou linéarisée selon le contrat de l'opérateur. Vérifier les flux face par face et leur annulation réciproque aux faces internes.
- AMG directionnel natif : construction de lignes à partir d'arêtes CSR fortes, substitution tridiagonale, Line-Gauss–Seidel multiplicatif, agrégation prioritaire le long des couplages faibles. Repli explicite vers l'agrégation classique si les lignes ne couvrent pas la structure. Aucun raccourci sur les résidus.
- Interpolation à énergie contrôlée : support à distance deux limité, conservation du mode constant `P 1=1`, pas d'optimisation borné selon une fonctionnelle d'énergie sur le prolongateur. La conservation de `P 1` n'implique pas à elle seule un bon spectre ni une accélération ; mesurer `A_c = P^T A P` et la complexité.
- Système couplé `[[M,G],[D,C]] [u,p]^T=[r_u,r_p]^T`. Schur de référence `S=C-D M^{-1} G` ; approximations par blocs (PCD, LSC, BFBt, Vanka ou Braess–Sarazin) n'ont pas la même action que `S` ou `S^{-1}`. Les oracles doivent distinguer les deux.
- Les variantes non symétriques restent **FGMRES opt-in** ; CG est réservé aux systèmes SPD et préconditionneurs compatibles avec sa preuve.

## 4. Traitement des frontières

- Pression Neumann pure : fixer explicitement la jauge/projection et vérifier la compatibilité du second membre, jamais réduire artificiellement le résidu en projetant sans mesurer le flux conservatif.
- Dirichlet, Neumann, Robin : assemblage et reconstruction des flux doivent utiliser les mêmes valeurs ou flux imposés ; vérifier les faces obliques et skewées indépendamment.
- WALL/EMPTY/périodique : distinguer les dimensions intrinsèques et les contraintes de cisaillement ; une méthode WLS fondée sur un stencil réellement tridimensionnel doit refuser une dimension intrinsèque erronée.
- Les blocs pression–vitesse doivent préserver les conditions de bord et le noyau du système selle de cheval ; un préconditionneur qui converge sur une matrice synthétique ne prouve pas leur traitement complet.

## 5. Pièges connus

- Confondre normes physique et mise à l'échelle (le résidu initial de 400 peut diminuer après une résolution plus précise ; un test doit recalculer `Ax-b` sur le **résultat**).
- Confondre interpolation améliorée et solveur plus rapide : le setup peut devenir dominant sur 1 024 inconnues.
- Faire confiance aux sorties `Status=CONVERGED` si le résidu original ou les bilans de conservation n'ont pas été recalculés.
- Confondre source SST et opérateur WLS : les deux chemins du code peuvent être affine-exacts ; exiger une différence strictement positive entre deux solutions affines est un oracle invalide.
- Recyclage numérique de lignes après rotation de la direction dominante ; hiérarchie transactionnelle à préserver en cas de rejet.
- Confondre l'action d'un préconditionneur Schur avec l'inverse du Schur exact ; les formulations LSC/BFBt nécessitent leurs hypothèses de stabilité propres.

## 6. Ordre de précision attendu et validation exigée

Les préconditionneurs n'ont **aucun ordre spatial** à revendiquer. La discrétisation FV doit être qualifiée indépendamment par MMS affine et quadratique sur trois maillages et faces skew/non orthogonales, puis par budgets locaux. Pour chaque famille de matrices (alignée, oblique 15/30/45°, directions tournantes, ratio 1/10/100/1 000/10 000), publier : `setup`, `apply`, `solve`, temps total, itérations, nombre de niveaux, complexité `NNZ(A_l)`, vraie erreur de solution lorsqu'elle existe, résidu vrai et causes de fallback.

Pour COUPLED, vérifier numériquement `S v` **et** l'action du solveur `S^{-1}b` sur des blocs de référence, ainsi que les bilans physiques pression-vitesse. Tests finaux NACA N64/N128/N256, Ghia, canal turbulent SST avec références indépendantes, convergence résiduelle, conservation et GCI. MPI/GPU à qualifier distinctement. Les tolérances existantes ne sont pas assouplies par cette fiche.

## 7. Statut d'implémentation CFDX

**Implémentation présente dans l'archive :**
- `src/cfdx/core/linalg/algebraic_line_preconditioner.h`, `anisotropic_line_amg.h`, `amg_preconditioner.h`, `hypre_amg.cpp` : lignes, AMG et interpolation énergie (natif, opt-in).
- `src/cfdx/core/linalg/coupled_amg_schur.h`, `vanka_patch_preconditioner.h`, `pcd_schur.h`, `lsc_bfbt_schur.h` : préconditionnement couplé (certification incomplète sur matrices CFD réelles).
- `src/cfdx/physics/sst_solver.h`, `scalar_transport_gradient.h`, `finite_volume_transport.h` : sources SST, WLS et flux de transport.

**Preuves locales de la présente correction CI :** tests CTest ciblés SST/skewness et harness, plus compilation/exécution du test transport avec Clang (voir rapport du lot). Ce ne sont ni les tests complets de #1182 sur GitHub, ni une campagne de validation CFD convergée.

**Bloquants :** replay du commit exact `012840d1`, CI Release + DebugSanitizers + Python, campagne N9, matrices réelles NACA/COUPLED, reproductibilité de la vitesse MPI et qualification physique SST. L'acceptation de cette fiche nécessite une revue humaine ; n'activer aucune nouvelle méthode automatiquement.

## 8. Journal

- 2026-11-06 : création sur archive locale lors du traitement des échecs CI signalés dans le document joint ; **DRAFT**.