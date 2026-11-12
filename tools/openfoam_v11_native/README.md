# Native OpenFOAM Foundation 11 wall-function matrix oracle — experimental

The `cfdxOmegaWallSnapshotV31.H` is **instrumentation for OpenFOAM itself**;
this CFDX environment does not contain a native OpenFOAM v11 installation.
It must be compiled inside the exact OpenFOAM Foundation 11 solver sources
(e.g. a local research copy with a `fvScalarMatrix omegaEqn` in scope), after
`omegaEqn.boundaryManipulate(omega.boundaryFieldRef())` and before the solve.
Call `cfdxOmegaWallSnapshotV31(mesh, omega, G, nut, omegaEqn, "native.json")`.
Do NOT calculate coefficients from CFDX formulas or export after an alternate
matrix stage. This hook is **not tested against OpenFOAM's compiler here**;
OpenFOAM code integration and matrix-algebra validation are outstanding.

The original instrumentation exports owner cells, wall faces, corner-face counts,
`omega`, `G`, `nut`, face coverage and the Foundation 11 corrected weights,
and the source/diagonal/full LDU entries for all touched rows. The source
commit and actual native log are recorded separately by
`package_native_snapshot_v31.py`. Use identical canonical mesh/state snapshots
on the two codes; a common hash is **not** proof of matching physical setup
unless the canonical artifacts actually encode identical ordered DOFs/fields.

Run the existing `scripts/qualify_openfoam_wall_v30.py native.json cfdx.json
--native-run-log openfoam.log`. The comparator rejects absent fields, rows,
wrong weights, mismatched provenance and numeric errors. It only certifies
complete input *on the declared scope*; independently verifying completeness
of mesh topology and stage is part of the native integration test. No native
OpenFOAM execution or successful parity claim is made by this deliverable.

Reference (Foundation v11, NOT the OpenCFD fork):
https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html

The header uses original CFDX-owned instrumentation code and does not contain
copied OpenFOAM implementation code.

## V32 : contrôle de disponibilité de l'environnement natif

Dans le shell **OpenFOAM Foundation 11 réellement chargé**, exécuter :

```bash
python3 tools/openfoam_v11_native/check_native_environment_v32.py
# Nécessite WM_PROJECT=OpenFOAM, WM_PROJECT_VERSION=11,
# WM_PROJECT_DIR, wmake et foamRun.
```

Le statut `NATIVE_RUNTIME_AVAILABLE_NOT_YET_QUALIFIED` est uniquement une
précondition. Il faut ensuite **compiler** l'instrumentation `.H` V31 dans un
solveur Foundation 11, exécuter le cas avec les champs / maillage / ordre des
cellules canonisés, exporter la matrice `omegaEqn` après manipulation des
conditions limites, puis comparer avec une capture CFDX indépendante. Toute
capture synthétique ou issue du code CFDX ne constitue **pas** un oracle natif.

Si le précontrôle renvoie `NATIVE_RUNTIME_UNAVAILABLE`, la qualification
OpenFOAM est impossible dans cet environnement : ne pas générer de référence
artificielle. Cette procédure n'installe pas OpenFOAM et ne se substitue pas à
une exécution de `foamRun` ou à une compilation par `wmake`.

# Native OpenFOAM Foundation 11 wall-function matrix oracle — experimental

The `cfdxOmegaWallSnapshotV31.H` is **instrumentation for OpenFOAM itself**;
this CFDX environment does not contain a native OpenFOAM v11 installation.
It must be compiled inside the exact OpenFOAM Foundation 11 solver sources
(e.g. a local research copy with a `fvScalarMatrix omegaEqn` in scope), after
`omegaEqn.boundaryManipulate(omega.boundaryFieldRef())` and before the solve.
Call `cfdxOmegaWallSnapshotV31(mesh, omega, G, nut, omegaEqn, "native.json")`.
Do NOT calculate coefficients from CFDX formulas or export after an alternate
matrix stage. This hook is **not tested against OpenFOAM's compiler here**;
OpenFOAM code integration and matrix-algebra validation are outstanding.

The original instrumentation exports owner cells, wall faces, corner-face counts,
`omega`, `G`, `nut`, face coverage and the Foundation 11 corrected weights,
and the source/diagonal/full LDU entries for all touched rows. The source
commit and actual native log are recorded separately by
`package_native_snapshot_v31.py`. Use identical canonical mesh/state snapshots
on the two codes; a common hash is **not** proof of matching physical setup
unless the canonical artifacts actually encode identical ordered DOFs/fields.

Run the existing `scripts/qualify_openfoam_wall_v30.py native.json cfdx.json
--native-run-log openfoam.log`. The comparator rejects absent fields, rows,
wrong weights, mismatched provenance and numeric errors. It only certifies
complete input *on the declared scope*; independently verifying completeness
of mesh topology and stage is part of the native integration test. No native
OpenFOAM execution or successful parity claim is made by this deliverable.

Reference (Foundation v11, NOT the OpenCFD fork):
https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html

The header uses original CFDX-owned instrumentation code and does not contain
copied OpenFOAM implementation code.

## V32 : contrôle de disponibilité de l'environnement natif

Dans le shell **OpenFOAM Foundation 11 réellement chargé**, exécuter :

```bash
python3 tools/openfoam_v11_native/check_native_environment_v32.py
# Nécessite WM_PROJECT=OpenFOAM, WM_PROJECT_VERSION=11,
# WM_PROJECT_DIR, wmake et foamRun.
```

Le statut `NATIVE_RUNTIME_AVAILABLE_NOT_YET_QUALIFIED` est uniquement une
précondition. Il faut ensuite **compiler** l'instrumentation `.H` V31 dans un
solveur Foundation 11, exécuter le cas avec les champs / maillage / ordre des
cellules canonisés, exporter la matrice `omegaEqn` après manipulation des
conditions limites, puis comparer avec une capture CFDX indépendante. Toute
capture synthétique ou issue du code CFDX ne constitue **pas** un oracle natif.

Si le précontrôle renvoie `NATIVE_RUNTIME_UNAVAILABLE`, la qualification
OpenFOAM est impossible dans cet environnement : ne pas générer de référence
artificielle. Cette procédure n'installe pas OpenFOAM et ne se substitue pas à
une exécution de `foamRun` ou à une compilation par `wmake`.
