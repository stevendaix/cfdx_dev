# TH-03 — V4 : maillage OpenFOAM six blocs, import natif, calculs négatifs

Date : 2026-11-03. **STATUT PHYSIQUE : NON QUALIFIÉ.**

## Référence primaire (v2606)

`https://gitlab.com/openfoam/core/openfoam/-/raw/OpenFOAM-v2606/tutorials/incompressible/simpleFoam/backwardFacingStep2D/system/blockMeshDict`

Le `blockMeshDict` comporte six blocs hexagonaux, `scale=0.0127`, les longueurs
`x/h = -130 .. 50`, `y/h = 0 .. 9`, et huit patches : `inlet`, `outlet`,
`lowerWallStartup`, `upperWallStartup`, `upperWall`, `lowerWall`, `front`, `back`.
Les résolutions natives des six blocs sont `(2,65,1)`, `(88,65,1)`,
`(97,48,1)`, `(33,48,1)`, `(97,65,1)`, `(33,65,1)` : **20 540 cellules**.
`front/back` sont `empty`; les startup walls sont `symmetryPlane`.

**Écart de contrat important :** la fiche expérimentale TH-03 du dépôt gèle
`Re_h=37 500` et `nu=1.4969066667e-5 m²/s`, avec un profil de couche limite
mesuré à `x/h=-4`. Le tutoriel OpenFOAM fixe `Re_h≈35 983` avec `nu=1.56e-5
m²/s` et sa propre turbulence/entrée. Il s'agit d'une **reproduction du cas
tutoriel OpenFOAM**, pas d'une reproduction expérimentale TH-03 déjà
qualifiée. `x_r/h=6.26` est seulement une valeur de comparaison documentaire
jusqu'à harmonisation complète du contrat.

Le générateur `scripts/generate_th03_openfoam_blockmesh.py` reproduit le
**domaine, les blocs, les nombres de cellules, les noms/types des patches**.
Les gradings `simpleGrading`/`edgeGrading` sont reconstruits par interpolation
bilinéaire et sous-segments géométriques indépendants ; cette reconstruction
n'est **pas la sortie exacte** de `blockMesh`. Pour conserver les faces et
sommets exacts, exécuter l'outil OpenFOAM `blockMesh` sur la machine cible, puis
utiliser `--openfoam-case /chemin/cas` pour importer les fichiers ASCII
`constant/polyMesh/{points,faces,owner,neighbour,boundary}`.
Un fixture `polyMesh` ASCII de synthèse a été importé de bout en bout, mais
aucun binaire `blockMesh` OpenFOAM n'était disponible dans cet environnement.

## Génération

```bash
PYTHONPATH=src/cfdx/python python scripts/generate_th03_openfoam_blockmesh.py \
  --resolution full --max-iterations 300 --output build/th03_v4/full.h5

PYTHONPATH=src/cfdx/python python scripts/audit_th03_openfoam_mesh.py \
  build/th03_v4/full.h5 --output build/th03_v4/mesh_audit.json
```

Import de l'exact `polyMesh` OpenFOAM, **si les fichiers existent** :

```bash
PYTHONPATH=src/cfdx/python python scripts/generate_th03_openfoam_blockmesh.py \
  --openfoam-case /chemin/vers/backwardFacingStep2D \
  --output build/th03_v4/imported.h5 --max-iterations 300
```

Calcul (SST encore expérimental, **un échec n'est pas une validation**) :

```bash
build/cfdx_production_solver --mesh build/th03_v4/full.h5 \
  --output-dir build/th03_v4/out --iterations 300 \
  --experimental-sst --benchmark-export
```

## Géométrie et indépendance des tests

| Maillage | Cellules | Faces | Faces inversées | Non-orthogonalité max |
|---|---:|---:|---:|---:|
| smoke, 6 blocs | 712 | 2907 | 0 | 8,403° |
| full, 6 blocs | 20540 | 82493 | 0 | 9,413° |

La validation géométrique contrôle les indices CSR, normales, aires des
cellules extrudées, couverture des patches, domaine et orthogonalité. Elle ne
certifie pas `y+`, la SST, le schéma de convection, ni la précision de `Cf`.
Le test `tests/python/test_th03_openfoam_blockmesh.py` comprend des
contre-exemples (face inversée, patch erroné, `polyMesh` binaire non supporté).

## Correctifs CFDX

1. **Flux de symétrie** : les deux calculs autoritatifs de débit massique
   (`make_mass_flux`, `make_rhie_chow_mass_flux`) imposent un débit normal nul
   sur `PatchType::SYMMETRY` comme sur `EMPTY`; test C++ dédié.
2. **Typage OpenFOAM** : `patch` générique n'est plus abusivement une paroi,
   les libellés `symmetryPlane`, `pressureOutlet`, etc. sont comparés après
   mise en minuscules. L'import TH-03 exige les huit noms/types attendus.
3. **Conditionnement testable** : deux réglages de cas optionnels sont
   exposés : `numerics.raw_settings.pressure_gmres_restart` et
   `numerics.raw_settings.pressure_linear_scaling`. L'acceptation du solveur
   continue de dépendre du vrai résidu de l'opérateur original ; aucune
   tolérance de convergence ni gate physique n'est desserrée.

## Campagne physique exécutée — résultats négatifs conservés

Propriétés reproduites : `h = 0.0127 m`, `U_inlet = 44.2 m/s`, `nu =
1.56e-5 m²/s`, `Re_h ≈ 35983`, `k_inlet = 1.09e-3 m²/s²`,
`omega_inlet = 181728 s⁻¹`. Les modèles, les profils et les conditions
pariétales ne sont pas encore numériquement équivalents au tutoriel OpenFOAM.

| Essai | Contrôle linéaire (maximum 1000) | Dernier résidu pression relatif | Verdict |
|---|---|---:|---|
| smoke 712 cellules, U=44.2, SST | FGMRES, restart=160, tol=1e-11 | 1,836e-11 | ÉCHEC |
| coarse 1280 cellules, U=44.2, SST | FGMRES, restart=160, tol=1e-11 | 2,142e-11 | ÉCHEC |
| full 20540 cellules, U=44.2, SST | FGMRES, restart=160, tol=1e-11 | 2,377e-10 | ÉCHEC |
| fixture ASCII polyMesh importé, 712 cellules | FGMRES, restart=160, tol=1e-11 | 2,092e-11 | ÉCHEC |
| smoke 712 cellules, laminaire diagnostique U=1 m/s | 3 itérations de transport achevées | non pertinent | NON CONVERGÉ, diagnostic de pipeline uniquement |

Les trois premières lignes sont des cas **à Re_h de référence** ; aucun des
calculs ne franchit le correcteur de pression, et **aucun `x_r/h` n'est
présenté comme valide**. La tentative de mise à l'échelle diagonale conserve
les vrais résidus physiques, mais n'a pas résolu le blocage. Le fait que le
même blocage existe en laminaire `U=44.2 m/s` confirme que ce défaut ne dépend
pas exclusivement de la fermeture SST. En revanche, un cas `U=1 m/s` n'est
pas le benchmark et ne permet aucune comparaison à Driver–Seegmiller.

Le solveur de production est recompilé ; 45 tests Python ciblés et cinq tests
C++ d'observables ont été exécutés avec succès. La compilation ne certifie
aucune des prédictions physiques.

## Verrous restants, par ordre de priorité

1. Réparer/certifier le solveur de correction de pression sur systèmes très
   étirés : valeurs propres/conditionnement, résidu **original** à chaque
   itération, réutilisation AMG applicable aux matrices SPD, stratégie de
   pression indépendante du solveur momentum, et preuve de conservation.
2. Reproduire **exactement** les conditions turbulentes et les fonctions de
   paroi OpenFOAM, ou adopter et qualifier explicitement une formulation SST
   réellement résolue jusqu'au mur avec un maillage `y+` approprié.
3. Implémenter / qualifier la convection bornée d'ordre deux effectivement
   utilisée à la fin de la campagne, non un simple schéma amont ordre un.
4. Obtenir des calculs convergés sur **trois maillages** et caractériser
   l'incertitude spatiale (GCI), avec les bilans locaux et globaux.
5. Extraire `Cf`, `Cp`, profils et réattachement par diagnostics indépendants ;
   comparer aux données Driver–Seegmiller dans les mêmes conditions de Reynolds
   et avec les incertitudes expérimentales associées.

**Un `PASS_GEOMETRY_ONLY` ne devient jamais un `PASS_PHYSICAL`.**
