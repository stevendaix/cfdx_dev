# CFDX Baseline

## État de référence

La branche `master` contient le socle M0 fusionné précédemment. La PR #3 (`audit/mesh-importers`) ajoute la validation approfondie des imports de maillages et reste la référence de travail jusqu'à son merge.

## Builds de référence

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure

cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=DebugSanitizers
cmake --build build-sanitize -j4
ctest --test-dir build-sanitize --output-on-failure
```

## Validation M0

Le statut `DONE` signifie que l'implémentation existe. Le statut de validation finale doit être confirmé par CI ; les anciens chiffres historiques ne sont pas utilisés comme preuve de l'état courant.

### M0.10 Import/export

La PR #3 ajoute :

- lecteur natif OpenFOAM `constant/polyMesh` ;
- dispatcher universel ;
- bridge meshio → CFDX-HDF5 ;
- Gmsh via meshio ;
- support des topologies courantes d'ordre supérieur par réduction aux sommets ;
- support des polyèdres lorsque meshio les expose ;
- mapping des groupes physiques vers les patches ;
- validation topologique ;
- tests OpenFOAM et meshio ;
- matrice géométrique cube / tetra / pyramid / wedge en Gmsh ASCII, Gmsh binaire et VTU ASCII.

## Règle de validation

Aucun résultat de test historique dans cette documentation ne doit être présenté comme une validation du commit courant. La source de vérité pour le statut CI est GitHub Actions sur le commit de la branche concernée.

## Prochaine étape

Après validation CI de la PR #3 et son merge : démarrer M1 — Navier-Stokes incompressible.
