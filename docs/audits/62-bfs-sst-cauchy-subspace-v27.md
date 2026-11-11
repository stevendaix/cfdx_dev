# CFDX V27 — Cauchy dans un sous-espace SST avec mérite FV projeté

Date : 2026-11-10. Base : CFDX V26 exacte. **Mode expérimental désactivé par défaut**.

## Problème observé

Sur le BFS diagnostique de 712 cellules, la direction Newton spatiale V26 satisfait le test linéaire GMRES (résidu vrai directionnel environ 0,1661, seuil *linéaire* 0,25), mais son L2 FV dynamique augmente de 0,178007285 à 0,178007476. Rechercher de nouveaux alpha le long de cette direction ne garantit donc pas une descente de l'objectif lisse, et le recul de `omegaWallFunction` peut faire varier le véritable opérateur.

## Nouvelle proposition mathématique

On conserve le **même** résidu FV physique reassemblé de production, ses fonctions de paroi, ses limites de positivité et son contrôle terminal Linf. La nouvelle fonction `sst_cauchy_subspace_v27` explore un sous-espace de dimension au plus quatre : direction logarithmique du Newton spatial, direction logarithmique Picard si admissible, et deux vecteurs canoniques au niveau des lignes FV les plus mauvaises de k et d'omega. Aucun indice BFS n'est codé en dur.

Pour l'objectif dynamique `Phi(z)=0.5*M2(exp(z))^2`, la projection des cellules proprietaires et la mise a jour de `nu_t/G/omega` sont incluses dans chaque evaluation de `M2`. On sonde, en difference centree, les derivees `g_j=dPhi(z+h*b_j)/dh` avec `h=1e-3`, puis une direction de Cauchy sous-espace `d = -sum_j g_j b_j` normalisee. Un autre candidat utilise le gradient mesure du Linf actif. Chaque direction est re-verifiee par une difference directionnelle centree du **merite physique projete** avant toute recherche de region de confiance.

Les rayons testés sont `0.5 / 2^j, j=0..7`. La norme infinie logarithmique du deplacement *apres projection de paroi* doit rester dans le rayon. Le rapport `(Phi_avant-Phi_essai)/(-rayon*derivee_Phi)` doit dépasser `1e-3`, ET **les deux mesures** `M2` et `Linf` doivent diminuer strictement. Une proposition acceptée reste soumise à la sélection de production existante, avec réassemblage et `sst_nonlinear_step_reduces_merit`. Le retour Picard et le candidat Newton-Jacobi V22 restent disponibles. Ceci est une adaptation matriceless de Cauchy sous-espace; **ce n'est ni SNESNEWTONTR ni une reproduction de la théorie de région de confiance de PETSc**.

## Expérience sur le BFS de 712 cellules

Même HDF5 V26 excepté le drapeau `sst_subspace_trust_v27_experimental=true`; budget diagnostique 10 itérations SST, tolerance physique nominale `1e-8` **inchangée**.

À l'itération interne 10 du **premier essai extérieur** :

| Indicateur | V26 | V27 candidat Cauchy |
|---|---:|---:|
| Linf avant | 0.92596107454 | 0.92596107454 |
| Linf candidat | 0.92595972137 | 0.92178503383 |
| L2 avant | 0.17800728499 | 0.17800728499 |
| L2 candidat | 0.17800747609 (**hausse**) | 0.17516771348 (**baisse**) |
| Rayon | — | 0.5 |
| Ratio réel/prédit | — | 0.7774181139 |
| Nombre de bases | — | 4 |
| Réassemblages additionnels V27 | — | 29 |
| Cellules actives selon le résidu physique | — | k:210, omega:284 |

**Limite décisive** : le Newton-Jacobi V22 déjà sauvegardé obtient un Linf final encore meilleur de **0.88990745136** et reste sélectionné. La sortie n'est PAS convergée à 10 itérations, et V27 n'apporte pas encore de meilleur état final. Le calcul extérieur est rejeté, avec rollback. Aucune preuve à 800 itérations ni sur 20 540 cellules.

## Validation logiciel et limites

- Code de production compilé et lié (GCC 14.2.0, diagnostic `-O0`, sans mesure de performance Release), drapeau lu du HDF5 et chemin ancien conservé quand il est désactivé.
- Contrôles unitaires C++ V27 (8 cas) : objectif réassemblé, double porte, baseline périmée, paroi projetée hors trust, certificats invalides, état appelant inchangé, indices/floors, dépendance opt-in.
- Protocoles Python V27 : refuse certificats falsifiés et faux état « convergé ».
- Rejeu A/B : aucun seuil changé ; le seul facteur de configuration est le drapeau V27.
- Le coût est **important** : jusqu'à environ 29 réassemblages supplémentaires par itération échantillonnée, non optimisé CPU/GPU.
- Le modèle de prédiction de la réduction est seulement de premier ordre et les directions de Cauchy ne sont mesurées que dans un sous-espace restreint ; pas de jacobien complet, pas de qualification formelle d'une stratégie trust-region PETSc.
- La fonction propriétaire `omegaWallFunction` et ses pondérations `updateCoeffs` / `manipulateMatrix` ne sont pas encore équivalentes à OpenFOAM.
- BFS 20 540, six variantes Schur, NACA, SST à 800 itérations, NASA Driver–Seegmiller, CI globale : **non qualifiés dans ce lot**.

## Reproduction

```sh
cmake -S . -B build-v27 -G Ninja -DCMAKE_BUILD_TYPE=Release \
   -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_TESTS=ON \
   -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v27 --target cfdx_production_solver test_sst_trust_subspace_v27
ctest --test-dir build-v27 -R 'test_sst_trust_subspace_v27|test_bfs_sst_v27_protocol' --output-on-failure
# Ajouter au numerics.raw_settings du HDF5 identique a V26 :
# "sst_subspace_trust_v27_experimental":"true"
CFDX_SST_DIAGNOSTICS=1 build-v27/cfdx_production_solver \
  --mesh bfs_v27_on_10.cfdx.h5 --output-dir out_v27 \
  --iterations 1 --experimental-sst
python3 scripts/diagnose_bfs_sst_v27.py v27_on.log
```

### Références méthodologiques

- PETSc `SNESNEWTONTR` : https://petsc.gitlab.io/petsc/release/manualpages/SNES/SNESNEWTONTR/
- PETSc `SNESNEWTONLS` : https://petsc.gitlab.io/petsc/release/manualpages/SNES/SNESNEWTONLS/
- OpenFOAM Foundation v11 `omegaWallFunction` : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html
