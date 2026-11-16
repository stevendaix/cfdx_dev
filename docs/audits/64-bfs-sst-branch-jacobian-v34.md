# CFDX V34 — dérivées SST par branche et mesures du Jacobien intercellulaire

Date : 11 novembre 2026. Base : V33 exacte (archive sources).

**Statut : expérimental, non convergé physiquement, défaut OFF.**

## Hypothèse auditée

Le changement de branche du limiteur de production SST en cellule 341 aux
itérations intérieures 3/4 (V33) pouvait perturber la diagonale du
préconditionneur k–omega V23/V24. V34 compare donc les blocs 2x2 locaux aux
véritables dérivées du résidu FV réassemblé, en conservant la référence
stationnaire et les contraintes de paroi de V33.

## Ce qui a été raccordé au solveur de production

- `sst_branch_preconditioner_v34_experimental` : option HDF5, désactivée par
  défaut, dépend des options expérimentales V22 JFNK et V24 spatial.
- Sélection automatique des cellules aux plus grands résidus FV normalisés
  dans chacune des équations `k`, `omega`, ou sélection de deux cellules
  diagnostiques `CFDX_SST_V34_CELL_IDS=341,329`.
- Dérivation des vrais quatre coefficients locaux dans les variables
  `log(k), log(omega)` et avec les mêmes échelles figées que le Newton.
- Vérification des branches SST de V33 : refus explicite aux planchers ou
  près du changement de `max`/`min` ; comparaison des dérivées unilatérales
  avant toute insertion dans le préconditionneur.
- Si l'écart est inférieur à 1 %, **conservation du bloc V23**. Si le bloc
  réassemblé est singulier/quasi singulier : **conservation du bloc V23**.
- La direction candidate V24 continue à être comparée à Newton-Jacobi et
  Picard ; **le vrai résidu physique réassemblé décide**. Le Jv Newton,
  les flux Rhie–Chow, les BC et la limite `1e-8` sont inchangés.
- `CFDX_SST_V34_CROSS_AUDIT=341:425` : oracle optionnel du véritable bloc
  intercellulaire 2x2, confronté au transport gelé V24. Cet oracle **ne
  remplace aucun coefficient** et ne prétend pas qualifier OpenFOAM.

Les calculs de différence finie sont plafonnés à 4 096 cellules, deux
hotspots et cinq assemblages FV par cellule à chaque probe (itérations 1,10,
puis cadence V28 éventuelle). Le grand BFS 20 540 cellules exige une stratégie
scalable par coloration, échantillonnage top-K et buffers préalloués.

## Résultats sur le BFS de 712 cellules

Les deux HDF5 ont exactement neuf jeux de données identiques ; seul le
paramètre V34 diffère. Le même binaire V34 est employé. Première tentative
extérieure et limite intérieure de **dix itérations**, pas qualification longue.

| Mesure | V34 OFF | V34 ON, cellules 341/329 |
|---|---:|---:|
| Mérite physique L_inf à itération 10 | 0.889874349 | 0.889874349 |
| Newton–Jacobi retenu | oui | oui |
| Blocs corrigés à itération 10 | — | 0 |
| Écart maximal local au vrai Jacobien, itération 1 | — | 1.09737e-5 |
| Écart maximal local au vrai Jacobien, itération 10 | — | 2.78159e-3 |
| Différence entre les deux mérites | — | **0** |
| Première itération extérieure engagée | non | non |

Sélection automatique supplémentaire : cellules 210/284 à itération 10,
écart maximal 2.54082e-5, bloc inchangé ; à itération 1, une non-régularité
unilatérale est refusée (sélection 347/55).

### Vérité du couplage 341 / 425

Le produit de dérivation réassemblé 4x est comparé aux coefficients de
convection/diffusion gelés réellement assemblés dans V24.

| Bloc intercellulaire, ligne cellule 341 colonne 425 | Iter. 1 | Iter. 10 |
|---|---:|---:|
| Écart relatif maximal, échelle du système Newton | 6.91926e-5 | 8.05593e-4 |
| dRk(341)/dlog(omega(425)) | 0 | -2.99379e-5 |
| dRomega(341)/dlog(k(425)) | -6.61142e-8 | -4.24279e-5 |

Le fait que ces **deux blocs particuliers** soient assez fidèles n'autorise
pas à conclure que toutes les dérivées hors diagonale du Jacobien sont
qualifiées. Une contraction insuffisante peut venir d'une combinaison de
plusieurs cellules, du couplage avec U/p, de la contrainte de paroi, des
candidats et/ou de la globalisation elle-même.

## Limites de qualification

- Aucune convergence SST complète (seuil physique 1e-8) ni état extérieur
  U,p,k,omega accepté ; première tentative rejetée avec rollback.
- Un coefficient préconditionneur fidèle **ne prouve pas** la précision
  du Jv global ; ce dernier doit être comparé à des perturbations
  simultanées de U,p,k,omega, flux Rhie–Chow et BC au même état.
- La moyenne des coins et le traitement G/omega/matrice OpenFOAM Foundation
  11 ne sont pas qualifiés faute d'une exécution OpenFOAM native.
- Les six variantes COUPLED du BFS 20 540 et les campagnes NASA/NACA,
  thermique, GPU, MPI n'ont pas été exécutées dans ce lot.

## Reproduction

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build --target cfdx_production_solver \
  test_sst_branch_preconditioner_v34
ctest --test-dir build -R 'v34' --output-on-failure
CFDX_SST_DIAGNOSTICS=1 CFDX_SST_V34_CELL_IDS=341,329 \
  CFDX_SST_V34_CROSS_AUDIT=341:425 \
  build/cfdx_production_solver --mesh bfs_v34_on.cfdx.h5 \
  --output-dir bfs_v34_out --iterations 1 --experimental-sst
```

Le paramètre HDF5 `sst_branch_preconditioner_v34_experimental=true`
nécessite `sst_spatial_block_preconditioner_v24_experimental=true` et
`sst_jfnk_v22_experimental=true`.

## Références

- PETSc SNES, séparation vrai Jacobien/approximation de préconditionnement,
  JFNK et test des dérivées : https://petsc.org/release/manual/snes/
- NASA TMR, SST-2003 et production limitée :
  https://tmbwg.github.io/turbmodels/sst.html
- OpenFOAM Foundation v11, omegaWallFunction updateCoeffs/manipulateMatrix :
  https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html
