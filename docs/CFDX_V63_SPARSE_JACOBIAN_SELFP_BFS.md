# CFDX V63 — Jacobien physique creux et Schur selfp natif (recherche)

Date : 13 novembre 2026. Base : V62 reunifiee BFS/SST + rayonnement HO–LO/DOM V57. Pas de changement des seuils, pas de promotion automatique.

## Objectif et statut

Deux fonctionnalites C++ de production, accessibles **uniquement sur opt-in**, sont ajoutees :

1. `SparseLocalPhysicalJacobianV63` : construction creuse de l'action du Jacobien du **residu physique reevalue** `F_phys(U,p,k,omega)` par coloration du voisinage du maillage reel. Toutes les equations momentum, pression/Rhie–Chow, SST et viscosite participent a `F_phys` ; les effets locaux manquants sont recherches par une verification indépendante hors stencil. Ce n'est **pas** un Jacobien analytique par flux ni une qualification Newton sur BFS.
2. `NativeSelfpILUSchurV63` : septieme preconditionneur de recherche construit a partir de la matrice COUPLED `4N` reelle : `S_tilde=C-D diag(Auu)^(-1) G`, application block-LDU avec `ILU(0)` sur vitesse et pression. Les six variantes historiques PCD/LSC/BFBt/BlockLocal/SIMPLE-Schur/SIMPLEC-Schur restent inchangées et sont rejouées separement.

Le solveur refuse toujours les colonnes non lisses et les erreurs de fidélité ; aucun seuil existant n'a ete assoupli. Les deux modes sont OFF par defaut, ne constituent pas une validation NASA et ne doivent pas etre commercialises comme solveurs certifies.

## Resultats observes

- Test manufacturé (12 cellules, 72 inconnues) : 30 couleurs, 248 coefficients non nuls, ecart directionnel independant environ `7.95e-9`, derivées croisées momentum/k/omega presentes.
- BFS152, candidat initial, 912 inconnues : 97 couleurs ; sans ancrage 194 re-evaluations Fphys, 16 396 coefficients, 221 colonnes rejetees, `max_sided=1.0396266`; avec ancrage 195 re-evaluations, 32 505 coefficients, 341 colonnes rejetees, `max_sided=0.929989`. **Newton non accepté dans les deux cas**. Les essais d'ancrage ne rendent pas le point non lisse differentiable.
- BFS152, matrice COUPLED 4N gelee (608x608, 7080 coefficients), selfp/ILU0 C++ + FGMRES : **28 iterations, residu VRAI** `5.28750496399499e-10` sur une cible de diagnostic `1e-9`. Les six captures A/b/x0 sont identiques SHA-256. Ce resultat est strictement **algebrique et hors convergence physique**.
- BFS152, selfp en solveur de production sur une iteration extérieure : construction native, `schur_nnz=1212`, `triple_terms=10026`, **non convergence** momentum apres SST `1.89867`; BlockLocal de référence sur ce meme cas reste aussi non converge.
- Les six Schur de production sur une iteration : PCD bloque dans Krylov, les cinq autres echouent au controle physique apres SST ; aucun certificat BFS. L'instrumentation V63 OFF conserve le comportement de reference.
- 10/10 tests CTest ciblés; deux nouveaux tests C++ sous ASAN/UBSAN. Pas de CI complete; pas de BFS712/BFS20540; pas de qualification RANS/SST de parois actives dans OpenFOAM natif.

## Activation optionnelle

### Selfp dans la production

Le cas doit utiliser `schur.block_local` et `coupled_block_schur` (infrastructure 4N) ; le remplacement est explicite dans le journal `COUPLED_SELFP_V63`. Les six méthodes historiques ne sont pas remplacées par défaut :

```sh
CFDX_COUPLED_SELFP_ILU_V63=1 ./build/cfdx_production_solver \
  --mesh BFS152_block_local.cfdx.h5 --iterations 1 \
  --experimental-sst --output-dir /tmp/cfdx-v63-selfp
```

### Jacobien physique creux expérimental

En plus des conditions et options pré-existantes de V51/V52/V53/V57 :

```sh
CFDX_FPHYS_CANDIDATE_V51=1 \
CFDX_FPHYS_JFNK_V52=1 \
CFDX_FPHYS_JFNK_V52_FREEZE_PRESTEP=1 \
CFDX_FPHYS_JFNK_V53_BLOCK=1 \
CFDX_FPHYS_JFNK_V57_TRUTH=1 \
CFDX_FPHYS_JFNK_V63_SPARSE_LOCAL=1 \
./build/cfdx_production_solver --mesh BFS152_block_local.cfdx.h5 \
  --iterations 1 --experimental-sst --output-dir /tmp/cfdx-v63-sparse
```

L'option `CFDX_FPHYS_JFNK_V62_ANCHORED_JACOBIAN=1` active l'ancrage experimental commun. Sa fidelite echoue aussi sur BFS152. Les gardes font refuser les directions, sans activer Newton par defaut.

### Essais et audit de la matrice capturee

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build --target cfdx_production_solver \
  test_native_selfp_schur_v63 test_sparse_local_physical_jacobian_v63 \
  cfdx_native_selfp_replay_v63
ctest --test-dir build -R 'test_(native_selfp_schur_v63|sparse_local_physical_jacobian_v63|native_selfp_v63_contract)' --output-on-failure
python3 scripts/qualify_native_selfp_v63.py \
  --replay ./build/cfdx_native_selfp_replay_v63 \
  --captures /path/to/real-same-input-six-captures \
  --output /tmp/sel-fp-v63.json
```

Les captures sont issues de `scripts/run_bfs_six_schur_v37.py` avec `--fingerprint-operator-v40 --capture-original-operator`. Ce script V63 neutralise explicitement les opt-ins selfp et Jacobien dans ses **six** executions de comparaison. Une campagne de six dont certains tests échouent retourne un code de sortie non nul : ceci est attendu pour un constat de non-qualification et ne doit pas être converti en succes.

## Prochain verrous

- Construire les contributions **analytiques par face** pour la reconstruction limitee SOU/TVD/MUSCL, `rAU`, les flux Rhie–Chow, viscosité/tau et branches min/max SST ; choisir un membre du Jacobien généralisé et le figer pour toute la base Arnoldi.
- Relier ce membre a `F_phys` indépendant, qualifier fidélité/additivité et pas Newton effectivement **accepté**, puis seulement exploiter l'acceleration creuse.
- Comparer 7 préconditionneurs sur mêmes matrices et budgets lineaires sans confondre le 4N algebrique et le 6N Newton ; tester le comportement AMG anisotrope, MPI/CPU multicoeur puis GPU.
- Qualifier parois SST, BFS712/BFS20 540, puis campagnes multi-maillages Cf/Cp/reattachement avec oracles NASA/OpenFOAM indépendants. Aucune de ces preuves n'est fournie par V63.