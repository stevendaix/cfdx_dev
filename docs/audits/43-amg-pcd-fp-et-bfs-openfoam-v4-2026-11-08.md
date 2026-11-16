# CFDX — AMG PCD Fp : correction du grossissement et conservation du BFS OpenFOAM V4

**Date :** 8 novembre 2026. **Périmètre :** solveur natif AMG utilise comme preconditionneur FGMRES du seul `Fp` nonsymetrique de PCD, et controle de non-regression du maillage TH-03 OpenFOAM. Aucun seuil physique n'est modifie.

## Defaut AMG mesure dans les captures de production

L'AMG natif a pu construire 25 niveaux sur NACA N64 (`3096, 2092, 1594, 1368, 1225, ... ,171`). Ce grossissement de plus en plus faible produit une succession de corrections grossieres instables sur un Fp tres dissymetrique, avec un quasi-mode constant presque nul (pressions Neumann avec jauge d'une cellule).

Sur la premiere application N64, `||Fp z||/||r|| ~= 3.47e14`, alors que l'AMG tronque aux premiers niveaux garde un rapport proche de 1. Ce n'est pas un probleme resolu en ajoutant des iterations FGMRES.

**Correctif :** `MatrixFreeVcyclePreconditioner` possede une limite optionnelle de fraction de grossissement ; la valeur historique est `1.0`, donc ne modifie aucun autre appelant. Pour le seul Fp-COUPLED/PCD dont l'inverse est demande par FGMRES, `NativeAMGPreconditioner::configure_max_coarsening_fraction(0.90)` arrete la hierarchie avant le premier niveau dont le ratio atteint 90%. Le dernier niveau restant emploie l'ancien relaxateur clairseme (24 passes), sans fabrication de pseudo-inverse. L'operateur fin, les conditions limites, Kp/CG et les seuils sont inchanges.

Un garde optionnel `configure_bounded_nonsymmetric_action(1e6)` verifie `||A*M^-1*r||/||r||` sur la matrice originelle lors de chaque application et ne remplace la direction par Jacobi amorti que si l'amplification est catastrophique. Ce mecanisme non lineaire est exclusivement utilise avec **FGMRES**, pas CG. Le repli reste verifie et peut echouer explicitement. Une protection plus agressive a ete essayee puis rejetee parce qu'elle deteriorait les solveurs qui convergeaient. Aucun changement global des autres utilisateurs AMG.

## Rejeux reels de captures CSR + RHS + x0 (tolerance 1e-8)

| Capture | AMG precedent: residu vrai relatif | AMG Fp corrige: residu vrai relatif | ILU0: residu vrai relatif | Interpretation |
|---|---:|---:|---:|---|
| NACA N64 ordinaire, 3096 inconnues | 0.9951 | **4.6466e-9** | 1.9246e-9 | AMG corrige converge, 42 iterations |
| NACA N128 ordinaire, 12336 inconnues | ~1 | **5.9133e-9** | 6.1567e-9 | AMG corrige converge |
| NACA N256 ordinaire, 49248 inconnues | 7.8147e-5 | **1.1939e-9** | 8.6341e-9 | AMG corrige converge, 23 iterations |
| N9 Couette petit Fp | 4.5369e-10 | **4.5369e-10** | 7.9882e-9 | resultat anterieur conserve |
| NACA N64 difficile | ~0.99993 | **2.6314e-2** | 1.2468e-3 | PAS converge |
| NACA N128 difficile | ~1 | **1.0219e-2** | 1.1679e-3 | PAS converge |

Ces tests concernent les **resolutions internes Fp** sur les matrices capturees, et non N64/N128/N256 COUPLED converges ou les efforts aerodynamiques. La capture difficile N256 n'a pas termine son replay dans le budget local ; aucune reussite n'est revendiquee. Sur le solveur de production NACA recompilé, l'instrumentation confirme que des seconds membres Fp restent defavorables et atteignent le plafond 100 iterations. L'AMG corrige n'est pas encore une qualification physique de COUPLED.

## BFS — travaux OpenFOAM preexistants a conserver

Les fichiers `scripts/generate_th03_openfoam_blockmesh.py`, `scripts/audit_th03_openfoam_mesh.py` et `scripts/initialize_th03_hybrid.py` sont **deja** des travaux importants dans CFDX ; ils ne sont ni supprimes ni remplaces. Le generateur V4 utilise les **six blocs** de `backwardFacingStep2D`, avec **20540 cellules**, bornes `x/h=-130..50`, `h=0.0127 m`, patches `inlet`, `outlet`, `lowerWallStartup`, `upperWallStartup`, `upperWall`, `lowerWall`, `front`, `back` et extrusion `EMPTY`.

`tests/python/test_th03_openfoam_existing_mesh_regression.py` verrouille ces valeurs, le type des patches et la distinction de Reynolds. La recreation complete HDF5 du maillage V4 a ete effectuee dans ce lot ; l'audit retourne `PASS_GEOMETRY_ONLY`, 0 face interne ou frontiere inversee, et conserve `NOT_QUALIFIED`.

**Attention :** la recreation native de la repartition `blockMesh` n'est pas identique sommet par sommet au fichier `polyMesh` OpenFOAM ; l'importeur `--openfoam-case` existe pour cet usage. `blockMesh` n'etait pas installe dans l'environnement de cette passe. Le tuto OpenFOAM fixe `U=44.2 m/s`, `h=0.0127 m`, `nu=1.56e-5 m²/s`, donc `Re_h=35983.3333`; le contrat **experimental** TH-03 du depot fixe `Re_h=37500`, et exige un profil d'entree turbulent mesure. Les campagnes a 192/402 cellules avec entree uniforme ne sont **pas** celles de ce maillage V4 et leurs `x_r/h` ne le qualifient pas.

La campagne OpenFOAM V4 a ete lancee historiquement mais son solveur pression a echoue sur le maillage full, avant qualification SST. Reprendre cette campagne a geometrique et CL identiques, verifier d'abord le correcteur de pression et les vrais residus avant de comparer le reattachement `6.26 +/- 0.10`, `Cf`, `Cp`, et les stations `x/h=-4,1,4,6,10`. Ne jamais publier une comparaison d'efforts a partir d'un test de maillage.

## Regression / limitations

- Compilation Release `test_naca0012_qualification` : **OK**.
- CTest C++ specialises : **8/8 PASS** (AMG, PCD, Coupled Schur, selection NACA, turbulence).
- Python BFS/OpenFOAM + PCD replay + matrice COUPLED : **24/24 PASS**.
- Generation HDF5 du maillage 20540 : **OK** ; audit `PASS_GEOMETRY_ONLY`.
- CI complete Release/DebugSanitizers : non executee.
- Convergence complete COUPLED/PCD NACA N64/N128/N256 : non demontree.
- Qualification SST/BFS Driver–Seegmiller et comparaison CFDX/OpenFOAM sur `polyMesh` byte-exact : non demontrees.

## Reproduction

```bash
cmake -S . -B build-amg -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF
cmake --build build-amg --target test_pcd_fp_replay test_naca0012_qualification -j 2
ctest --test-dir build-amg -R 'test_pcd_fp_replay_capture|test_pcd_schur|test_amg_preconditioner_qualification' --output-on-failure
# Avec une vraie capture Fp en entree, compare amg ET amg_robust :
build-amg/test_pcd_fp_replay /chemin/vers/naca_n64_rhs1.cfp 100 1e-8
PYTHONPATH=src/cfdx/python:. python -m pytest -q tests/python/test_th03_openfoam_blockmesh.py tests/python/test_th03_openfoam_existing_mesh_regression.py
PYTHONPATH=src/cfdx/python:. python scripts/generate_th03_openfoam_blockmesh.py --resolution full --output build-amg/bfs_v4.h5 --max-iterations 2
PYTHONPATH=src/cfdx/python:. python scripts/audit_th03_openfoam_mesh.py build-amg/bfs_v4.h5 --output build-amg/bfs_mesh_audit.json
```

Sources externes utilisees comme references de contrat :
- OpenFOAM, `backward-facing-step`: https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/
- NASA/TMR, `2DBFS`: https://tmbwg.github.io/turbmodels/backstep_val.html
- hypre, BoomerAMG: https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html
- PETSc, solveurs et AMG: https://petsc.org/main/manual/ksp/
