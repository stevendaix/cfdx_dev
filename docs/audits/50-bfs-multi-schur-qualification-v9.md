# CFDX — Lot V9 : campagne transverse BFS / six Schur, contrat de preuve

**Date : 9 novembre 2026.** Base : archive CFDX `cfdx_dev-master (2)(1).zip`, SHA256
`c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c` ;
V8 appliquée sans conflits avec le patch cumulatif identifié par
`c345faafb2dadff97d79b644c4dc0cf766b17636dcc5618cf205ad88d6594bc0`.
Ne pas assimiler ces sources à un HEAD GitHub non vérifié.

## 1. Portée exacte et état de la qualification

Ce lot ne change **aucune tolérance de convergence**, n'active pas de solution
forcée, et ne valide pas physiquement les méthodes par déclaration.
La campagne BFS existante n'exécutait que `test_bfs_quick`; elle ne faisait pas
fonctionner les six variantes Schur du solveur de production en conditions SST.
Ce lot ajoute le lancement effectif de ces configurations et le contrôle
indépendant de leurs preuves persistées. **Un lancement possible n'est pas une
convergence mesurée**.

| Preuve | Statut dans cette livraison |
|---|---|
| Vérification génération de la même géométrie pour douze variantes | Tests Python PASS |
| Sélection et contrat d'export C++ | CTest ciblés PASS |
| Exécution vraie de `cfdx_production_solver` sur BFS | **NON EFFECTUÉE** — compilation complète du gros binaire non terminée dans cette session |
| Fermeture du résidu indépendant U,p,k,omega | **NON DÉMONTRÉE** pour les six Schur |
| Résidu vrai du bloc linéaire, conservation et SST | **NON QUALIFIÉS** sur la nouvelle matrice de campagnes |
| Réattachement, Cf, Cp, vitesse, contraintes de Reynolds | **NON ÉVALUÉS** |
| Étude NACA, N9, maillage CPU/GPU | **NON RÉEXÉCUTÉE** dans ce lot |

## 2. Matrice des configurations

`SIMPLE`, `SIMPLEC`, ainsi que :

| Schur COUPLED | Préconditionneur initialement demandé |
|---|---|
| PCD | Auto, PCD, bloc Schur/AMG |
| Block-local | Bloc Schur/AMG, MGR, Braess–Sarazin |
| LSC | Bloc Schur/AMG |
| BFBt | Bloc Schur/AMG |
| SIMPLE-Schur | Bloc Schur/AMG |
| SIMPLEC-Schur | Bloc Schur/AMG |

12 variantes distinctes, **10 configurations COUPLED**. L'identité de la
méthode **effectivement résolue** est exigée via `linear_plan` (dont le solveur
Krylov et le préconditionneur), et ne peut être remplacée par la sélection
`schur` demandée dans les entrées. Les configurations non compatibles sont
rejetées à la génération. Le registre CFDX classe encore
`preconditioner.pcd` « Planned », malgré sa présence dans le dispatch : il
faut le requalifier sur un calcul réel avant de changer ce statut documentaire.

## 3. Trois niveaux indépendants de preuve

1. **Exécution** : fichier `execution.json` strict, code de retour et champ
   `process_exit_code` cohérents avec le processus courant.
2. **Convergence non linéaire** : `convergence.json` effectivement convergé,
   dernier état accepté, `linear_plan` réellement résolu, certificat indépendant
   `global_nonlinear_v8` avec U,p,k,omega, fermeture viscosité turbulente,
   positivité, conservation et opérateur d'ordre deux sélectionné. L'absence de
   champ, les valeurs null/non finies et un Schur ou préconditionneur différent
   produisent un **FAIL**.
3. **Validation physique** : nécessite de véritables profils d'entrée et
   conditions de paroi NASA/Driver–Seegmiller, au moins trois maillages,
   `Cf(x/h)`, `Cp(x/h)`, profils de `U`, contraintes de Reynolds et réattachement.
   **Jamais inférée de 1 ou de 2.** Aucun de ces oracles n'est ajouté
   artificiellement par ce lot.

La métrique `full_order_eligible` seule ne distingue pas une solution convergée
sur une convection UPWIND d'une solution sur un opérateur convectif réellement
d'ordre deux. Le générateur de diagnostic choisit désormais explicitement
`numerics.convection.second_order_upwind` ; le contrôleur recalcule
l'identité du schéma à partir du HDF5 et refuse UPWIND.

## 4. Cas de diagnostic contre cas de référence

Le maillage rapide `generate_th03_bfs_case.py` est un cas six-blocs à inlet
uniforme (intensité turbulente 5 %), **non équivalent au profil de couche limite
expérimental**. `--bfs-reference-case /chemin/vers/cas_BFS.h5` permet de
préparer les mêmes douze variantes sur un vrai cas importé, en copiant
l'intégralité des datasets et en ne changeant que les sélections numériques,
le budget d'itérations et l'activation des deux diagnostics V7/V8. Les données
initiales, conditions aux limites, dimensions et champs du cas fourni restent
inchangés. Cette option ne garantit toutefois ni la fidélité du HDF5 importé
à l'expérience, ni sa convergence.

La référence Driver–Seegmiller publiée par OpenFOAM et TMR NASA donne
`Re_h≈36000`, `U≈44.2 m/s`, `h=0.0127 m` et `x_r/h=6.26 ± 0.10`.
Les jeux NASA TMR fournissent aussi `Cf`, `Cp`, les profils de vitesse et les
contraintes turbulentes. Références :

- https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/
- https://tmbwg.github.io/turbmodels/backstep_val.html

## 5. Reproduction

Installer l'environnement Python CFDX (`h5py`, `numpy`, `pydantic`, `pytest`).

```bash
PYTHONPATH=.:src/cfdx/python pytest -q \
  tests/python/test_bfs_schur_v9_campaign.py \
  tests/python/test_th03_bfs_case_generator.py \
  tests/python/test_coupled_matrix_campaign.py

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_ENABLE_GPU=OFF
cmake --build build --target cfdx_production_solver -j 4

# Matrice de diagnostic six Schur (pas une validation Driver-Seegmiller)
PYTHONPATH=.:src/cfdx/python python scripts/run_cfdx_coupled_campaign.py \
  --build build --out runs/bfs_debug --only bfs --iterations 1500

# Pour le cas importé : réutilise physiquement le même maillage et les mêmes BC
PYTHONPATH=.:src/cfdx/python python scripts/run_cfdx_coupled_campaign.py \
  --build build --out runs/bfs_reference --only bfs \
  --bfs-reference-case data/reference_bfs_driver_seegmiller.h5 \
  --iterations 1500
```

La commande de qualification termine avec un code non nul si un calcul a
échoué ou si le certificat V8 complet ne correspond pas aux résultats
persistés. `--bfs-reference-case` ne fabrique pas de référence NASA et exige
un véritable HDF5 CFDX déjà construit. Pour l'étude de maillage, répéter la
campagne avec au moins trois HDF5 de référence de résolutions différentes,
**sans assouplir les seuils** et calculer un GCI justifié.

## 6. Échecs historiques à conserver

Les observations du diagnostic initial (fournies par l'utilisateur, pas
reproduites par ce lot) étaient : SIMPLE pression `4.24e-11 > 1e-11`, COUPLED
+ PCD automatique `1.06e-2` relatif après 52 itérations, et une même panne
sur PCD en bloc Schur/AMG. La version V8 rapportait également un SST non
terminé et aucun `x_r/h` certifié. Ces observations ne justifient aucun
relâchement de tolérance ni substitution silencieuse de préconditionneur.

**Limites de cette livraison :** pas de binaire de production recompilé et lié
à partir de cette copie V9, pas de nouvelles simulations BFS réelles ni de
preuve d'accord expérimental. Toute affirmation contraire doit être rejetée.


## Bibliographie

- Driver, D. M.; Seegmiller, H. L. (1985). “Features of a reattaching turbulent shear layer in divergent channel flow.” *AIAA Journal*, 23(2), 163–171. DOI: [10.2514/3.8890](https://doi.org/10.2514/3.8890). NASA NTRS record: https://ntrs.nasa.gov/citations/19850039701.
