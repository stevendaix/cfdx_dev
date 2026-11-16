# CFDX V30 — pseudo-temps couplé SST et oracle OpenFOAM indépendant

Le pseudo-temps couplé V30 est une **option expérimentale** qui ne modifie pas les équations stationnaires physiques ni les seuils. Il construit la direction Newton du système assemblé `F(k,omega)` à l'aide du décalage `F(x) + D_tau(x-x_old)` où, pour chaque ligne, `D_tau=max(|A_ii|, 0.5 sum_j |A_ij|, 0.5 sum_faces |phi_f|)/CFL`. La ligne contrainte de `omegaWallFunction` ne reçoit pas de pseudo-masse. Ce décalage est une approximation algébrique de rayon spectral FV, pas une équivalence de pseudo-temps avec SU2 ou Fluent.

`sst_coupled_ptc_v30_experimental=true` nécessite JFNK V22. La direction V30 est soumise à la même sélection finale sur le vrai résidu SST non décalé et réassemblé, avec positivity et rollback. À chaque pas SST vraiment accepté (y compris via Newton–Jacobi ou Picard), le CFL est adapté ; il est réduit en l'absence de pas admissible. Le contrôleur n'est publié vers l'itération extérieure qu'après convergence du transport SST. Les transactions du solveur U/p/flux/SST sauvegardent et restaurent également ce contrôleur. L'option V29 et les autres préconditionneurs restent disponibles sans altération.

## Oracle indépendant OpenFOAM — contrat bloquant

`scripts/qualify_openfoam_wall_v30.py` compare deux exports JSON au schéma `CFDX_OPENFOAM_OMEGA_V30` : l'un calculé **par une vraie exécution OpenFOAM Foundation v11**, l'autre par CFDX sur le même maillage et **le même état figé**. Ils doivent déclarer les mêmes empreintes SHA-256 de géométrie et d'état. Les exports contiennent les mêmes IDs de cellules, faces et lignes matricielles. Les valeurs comparées comprennent `omega`, `G`, `nut`, le nombre de faces de coin, la couverture géométrique/pondération de chaque face, et tous les coefficients CSR de chaque ligne `omega` (diagonale, RHS, colonnes). La validation rejette les IDs manquants, les poids invalides, les coefficients absents, les états différents et les erreurs numériques hors tolérances. Elle refuse la certification si l'export de référence n'est pas déclaré natif avec un exécutable, une référence de cas et une révision source vérifiables.

Cet outil est un **comparateur indépendant**, pas un exécutable OpenFOAM ni une reproduction intégrale de `omegaWallFunction`. La propriété `native_export` est une déclaration de provenance insuffisante en soi. Le CLI exige aussi `--native-run-log` et vérifie l'empreinte SHA-256 du fichier présent contre `provenance.run_log_sha256`; sans cette preuve, il refuse `PASS_NATIVE`. Cette vérification d'intégrité ne garantit toutefois pas à elle seule l'authenticité scientifique du journal. En validation sérieuse, l'équipe doit attacher le journal de build, la révision, les fichiers du cas et leurs empreintes, ainsi que la trace de l'utilitaire d'export natif. **Sans ces pièces et les deux exports complets, l'équivalence reste NON QUALIFIÉE**. Les seules fixtures synthétiques vérifient le comportement à échec explicite du comparateur.

Exemple de ligne de commande :

```bash
python scripts/qualify_openfoam_wall_v30.py references/openfoam_v11_native.json references/cfdx_same_state.json --native-run-log references/openfoam_run.log --output preuves/oracle_wall_v30.json
```

Clés requises pour chaque JSON :

```json
{
 "schema": "CFDX_OPENFOAM_OMEGA_V30",
 "scope": {"cells": [357], "faces": [100], "matrix_rows": [357]},
 "mesh_sha256": "<64 caractères hexadécimaux>",
 "state_sha256": "<64 caractères hexadécimaux>",
 "provenance": {"engine": "OpenFOAM Foundation", "version": "11", "native_export": true,
                 "executable": "simpleFoam", "case": "BFS_canonical", "source_commit": "rev",
                 "run_log_sha256": "<SHA-256 du journal OpenFOAM natif>"},
 "cells": [{"cell": 357, "omega": 1.0, "G": 1.0, "nut": 0.0, "corner_faces": 2}],
 "faces": [{"face": 100, "owner": 357, "coverage": 1.0, "weight": 1.0,
            "omega": 1.0, "G": 1.0, "nut": 0.0}],
 "matrix_rows": [{"row": 357, "diagonal": 1.0, "rhs": 1.0,
                  "entries": [[357, 1.0], [363, 0.0]]}]
}
```

Pour l'export CFDX, remplacer `provenance.engine` par `CFDX` et renseigner les champs à l'assemblage de production identique ; **ne jamais fabriquer le JSON de référence à partir des données CFDX**.

Source de la méthodologie : [OpenFOAM Foundation v11 omegaWallFunction](https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html), [OpenFOAM fvMatrix](https://cpp.openfoam.org/v11/fvMatrix_8C_source.html), [PETSc TSPSEUDO](https://petsc.org/release/manualpages/TS/TSPSEUDO/), [SU2 RANS](https://su2code.github.io/docs_v7/Setup-Guide-RANS-External-Aero/). Ni tolérance physique ni politique de validation modifiée.

Le schéma impose `scope.cells`, `scope.faces` et `scope.matrix_rows` (listes d'IDs triés, sans doublon). Les trois listes doivent correspondre exactement aux données fournies; une face de paroi dont la cellule propriétaire ou la ligne matricielle manquerait est rejetée. Cela empêche de qualifier deux exports tronqués de manière identique, mais la preuve que le *scope* couvre toutes les cellules touchées dépend encore de la qualité de l'export OpenFOAM natif.

Le préconditionnement V30 a été renforcé pour équilibrer le Jacobien LOG du **système décalé** avec `(A_ii+D_tau)*x_i/(norme_ligne_i*échelle_physique_i)` (mêmes bornes numériques `0.02` et `50` qu'en V22). Ce changement ne touche ni aux équations stationnaires, ni aux critères de convergence; le vrai résidu restera la seule preuve de qualification.
