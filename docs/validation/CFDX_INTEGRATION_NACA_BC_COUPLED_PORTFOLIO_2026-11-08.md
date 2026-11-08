# CFDX - Integration NACA / BC / COUPLED / qualification (2026-11-08)

## Base et contrat de livraison

Source de reference: `cfdx_dev-master (2)(1).zip`, SHA-256
`c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c`.
Les autres archives et les rapports d'audit concernent des bases parfois differentes.
Le patch cumulatif correspondant a ce document vise **exclusivement** la source ci-dessus.
Ne pas forcer les hunks sur un HEAD GitHub evolue; rebaser et reexecuter les gates.

## Modifications integrees

1. NACA: correction de `surface(n, trailing_edge_regularization)`, campagnes
   quasi-isotropes, certification `certify(..., raw_faces)`, controle de
   l'integrale independante des efforts, concordance entre schema final et
   metadonnees, continuation UPWIND -> ordre superieur sans retour en arriere.
2. Contrat de campagne: le lanceur GCI choisit explicitement SIMPLE et le
   schema final MUSCL par defaut. Van Leer reste selectable explicitement;
   ni variable COUPLED issue du shell ni sorties CSV anciennes ne sont
   heritees. L'indicateur GCI *provisoire* N64/N128 exige lui aussi un schema
   identique; N64/N128/N256 et deux familles n'acceptent pas les schemas mixtes.
3. BC type Fluent: presets natifs, verification preflight, validation C++
   fail-closed, CLI natif `Case.run` et GUI. Les operateurs non implementes
   (slip vrai, mass-flow, pressure-inlet, periodic/interface, backflow
   universel) restent refuses. La GUI Qt complete et le pont GUI -> schema
   natif aller-retour restent des chantiers ouverts.
4. Qualification multi-domaines: registre de 25 cas, six candidats non
   inscrits, vue des contrats / campagnes / preuves physiques separees,
   workflows CI d'activation par profil. Les PASS historiques ne sont pas
   promus automatiquement.
5. COUPLED: le pilote NACA accepte un mode experimental avec selection Schur
   obligatoire (`block_local`, `pcd`, `lsc`, `bfbt`, `simple`, `simplec`).
   Les metadonnees d'historique mentionnent correctement `COUPLED`; le
   solveur verifie le Schur reellement utilise. Il s'agit du solveur
   monolithique existant et PAS d'une integration Newton-Krylov.
6. Tests Python: correction de la priorite d'import entre `scripts/cfdx.py`
   et le package Python natif, revelee par l'execution des suites ensemble.

## Resultats effectivement reproduits dans cette integration

- Python: 285 tests PASS, 3 SKIP, 7 sous-tests PASS (NACA, GCI, BC et suivi).
- CTest ciblage BC/NACA/PCD: 9/9 PASS; build NACA Release: reussi.
- Politique CI structurelle: PASS, 23 workflows et 8 partitions sous garde.
- Registre CTest avec `CFDX_ENABLE_LONG_VALIDATION=ON`: 393 cas inscrits,
  9 resultats JUnit importes, 0 manque d'ancrage; 0/25 preuves physiques.
  393 tests inscrits **ne signifie pas** 393 tests executes.
- NACA N64 Re 1000, alpha 0, SIMPLE + MUSCL + p-Neumann + traction faces WLS,
  alpha_U=0.6 alpha_p=0.25, tolerance lineaire 1e-12:
  `converged=1` apres 305 iterations, Cd=0.11195569962907033,
  Cdp=0.023750712460525777, Cdv=0.088204987168544455,
  momentum RHS relatif=7.7018106932817084e-11,
  farfield `kinematic_match=1`, ecart max flux/vitesse nul.
  Le certificat dit `ONE_GRID_ONLY`: pas de GCI prouve.
- NACA COUPLED + block_local, meme N64/MUSCL/Neumann, 3 iterations:
  selection Schur `exact_selection=1`, calcul non converge a 3 pas (attendu).
- NACA COUPLED + PCD, meme N64/MUSCL/Neumann: premiere resolution
  monolithique non convergee (status=2; 79 iterations internes,
  residu relatif env. 2.585405e-3). Le chemin est expose en diagnostique,
  PAS qualifie en physique. Un essai PCD Dirichlet s'est arrete sur timeout.
- Le build complet `cfdx_production_solver` Release a ete interrompu par
  timeout; ne pas lui attribuer de PASS. DebugSanitizers, CI GitHub,
  campagne N128/N256, GCI et GUI Qt: non executes dans ce lot.

## Pour reproduire

Appliquer le patch cumulatif sur la base exacte, puis:

```bash
cmake -S . -B build-integrated -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_LONG_VALIDATION=ON -DCFDX_ENABLE_MPI=OFF
cmake --build build-integrated --target test_naca0012_qualification \
  test_naca_coupled_selection test_production_boundary_contract -j 2
ctest --test-dir build-integrated --output-on-failure \
  -R '^(test_naca_coupled_selection|test_production_boundary_contract)$'
PYTHONPATH=src/cfdx/python python3 -m pytest -q \
  tests/python/test_naca*.py tests/python/test_qualification_portfolio.py \
  tests/python/test_qualification_evidence_gate.py \
  tests/python/test_fluent_boundary_presets.py \
  tests/python/test_native_solver_preflight.py
python3 scripts/ci/ci_policy.py
```

NACA SIMPLE/MUSCL qualifie sur **un maillage** (variables a definir depuis la
racine, chemins de maillage et d'executable a adapter):

```bash
CFDX_NACA_ALGORITHM=SIMPLE CFDX_NACA_LINEAR_TOLERANCE=1e-12 \
CFDX_NACA_CONVECTION_CELL_LIMITED=1 CFDX_NACA_FARFIELD_PRESSURE_NEUMANN=1 \
CFDX_NACA_REQUIRE_KINEMATIC_FACE_FLUX=1 CFDX_NACA_WALL_TRACTION_FACE=1 \
CFDX_NACA_WALL_GRADIENT_WLS=1 \
build-integrated/test_naca0012_qualification --single "$MESH_N64" 0.6 0.25 1000
```

Etude COUPLED **diagnostique** (ne transforme jamais une non-convergence en PASS):

```bash
CFDX_NACA_ALGORITHM=COUPLED CFDX_NACA_COUPLED_SCHUR=block_local \
CFDX_NACA_CONVECTION_CELL_LIMITED=1 CFDX_NACA_FARFIELD_PRESSURE_NEUMANN=1 \
build-integrated/test_naca0012_qualification --single "$MESH_N64" 0.6 0.25 3
```

Rapport de qualification transversale (un inventaire doit avoir les cibles
longues *inscrites*, meme si elles n'ont pas ete executees):

```bash
ctest --test-dir build-integrated --show-only=json-v1 > manifest.json
python3 scripts/qualification_portfolio.py --ctest-manifest manifest.json \
 --evidence-dir qualification_evidence \
 --output-dir out/qualification-portfolio --strict-links
```

## Limites pour promotion

Ne pas fusionner uniquement pour des tests synthetiques ou une certification
N64. Sur le HEAD exact de la PR: construire la production Release et
DebugSanitizers, executer tests Python/PyFluent et CTest complets, verifier
les vraies BC sur le CLI HDF5 natif, puis produire N64/N128/N256 converges
sur un meme operateur et une meme famille de maillages. Certifier Cp, Cf,
Cd, Cl, Cm et GCI seulement apres les oracles physiques independants.
Les campagnes SST, conduction/CHT, rayonnement, LES et VOF progressent en
parallele par leurs propres portes de validation.