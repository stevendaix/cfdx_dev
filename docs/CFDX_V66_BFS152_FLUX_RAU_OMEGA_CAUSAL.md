# CFDX V66 — BFS152, flux porteur, rAU, viscosité et cellule SST omega

**Date :** 13 novembre 2026. **Base cumulative :** ZIP BFS V6 exact fourni par l'utilisateur, évolutions V57–V65 conservées. **Statut :** audit causal C++ *en lecture seule*, uniquement si `CFDX_FPHYS_V66_OMEGA_FLUX=1` et si le mérite V65 est déjà activé. Aucune variante Schur, AMG ou Newton n'a été ajoutée ; les critères et modèles physiques sont inchangés.

## 1. Expérience, contrat de preuve

Le pilote est le cas HDF5 V65 `TH03_debug_only`, 152 cellules et 634 faces, sans qualification expérimentale Driver–Seegmiller ni équivalence des fonctions de paroi SST OpenFOAM. Les anciens solveurs COUPLED/BlockLocal et SST gardent leur discretisation, conditions limites, relaxation, limites de positivité et budgets.

La V66 mesure simultanément, **à chaque essai déjà calculé des itérations extérieures 22 et 23** :

- `phi_carrier` conservé par la transaction COUPLED et `phi_phys(X)` réévalué à partir des champs de cet état, y compris l'écart de divergence orientée par cellule ;
- les diagonales inverses de momentum `rAU` issues du couplage COUPLED et celles issues du point fixe Rhie–Chow du résidu physique ;
- la viscosité turbulente stockée après SST et celle réévaluée à partir du même candidat `U,k,omega` ;
- le classement des trois cellules qui maximisent respectivement l'évolution brute de R_omega, le défaut brut de R_omega du candidat et celui de l'état précédent ;
- aux pas bêta=1 et aux derniers pas de backtracking, un budget **volumes finis réassemblé** à la cellule la plus défectueuse : convection, diffusion, production, destruction, diffusion croisée, correction de paroi, somme et comparaison avec le résidu signé indépendant `F_phys`.

L'oracle refuse un maillage/flux non conforme, des coefficients non finis, un bilan SST non fermé, des indices de cellule non valides et un pseudo-certificat physique. L'option V66 exige V65 ; V65 exige lui-même V51, V47, V13/V14, l'absence de V48 et un maillage de 152 cellules au plus. **Le diagnostic ne modifie pas les états physiques acceptés.**

## 2. Résultat du rejeu réel (C++ diagnostique)

| Essai | V65 : décision, référence | Mesures V66 |
|---|---|---|
| itération 22, bêta=1 | Rejet : bloc omega, croissance 71,8× | Cellule omega 0 ; flux porteur/état-seul 8,662688e-4 |
| itération 22, bêta=1/256 | Acceptation V65 | Diminution du défaut de flux ; vérification omega disponible dans les traces |
| itération 23, bêta=1 | Rejet : bloc omega, croissance 126,95× | Cellule omega 0 ; flux porteur/état-seul 8,678287e-4 |
| itération 23, bêta=1/1024 | Rejet : rapport de mérite 0,998296 > 0,995 | Cellule omega 12 ; flux porteur/état-seul 6,545691e-7 |

La cellule 0 du cas réel a pour centre `(0,0211667, 0,0031750, 0,0063500) m`. Sur le bêta=1 de l'itération 23, son résidu omega passe d'environ `1,37e-15` à `-0,364717` ; sa décomposition signée (variation de résidu) vaut :

- convection : `-0,947909` ; diffusion : `+2,262611` ;
- production : `-8,417458` ; destruction : `+6,737691` ; diffusion croisée : environ `+0,00035` ;
- somme : `-0,364717`, fermeture numérique au niveau des arrondis.

La cellule 12 apparaît au petit pas de l'itération 23, avec un résidu passant de `+0,00287282` à `+0,00131514` et une variation de `-0,00155768`. Les termes variationnels sont : convection `-0,00403106`, diffusion `-0,00268027`, production `+0,01075597`, destruction `-0,00546820` et diffusion croisée `-0,00013412`.

**Interprétation bornée :** la cellule 0 (proche de la paroi sur le pilote) répond fortement aux sources de production/destruction et à la diffusion lors du grand pas. Le petit pas n'échoue plus sur le même foyer de omega ; il échoue sur la contraction du mérite global. Aucune de ces observations ne prouve une erreur de signe de source, un bug de condition limite ou une équivalence physique NASA.


### 2.1 Coefficients et interface géométrique critique

Sur le **même état candidat**, `rAU` du flux COUPLED diffère de `rAU`
reconstruit par l'évaluateur état-seul ; le premier est lié à la diagonale
issue de l'assemblage Picard ayant fabriqué la direction. Ce décalage
**n'établit pas un bug** tant que le contrat de l'opérateur n'a pas été
confronté, à état figé, aux assemblages réels.

| Comparaison au même état | Itération 22, beta=1 | Itération 23, beta=1 |
|---|---:|---:|
| Écart relatif maximal `rAU` porteur / état-seul | 0,4927673 | 0,4934951 |
| Cellule, composante du maximum | 0, `Ux` | 0, `Ux` |
| Écart `rAU` du dernier état accepté / état-seul | 0,0062199 | 0,0035656 |
| Écart `nu_t` stocké / `nu_t` état-seul au candidat | 0 | 0 |
| Variation maximale de `nu_t` entre ancien et candidat | 0,6363837 | 0,6344245 |
| Écart relatif flux porteur / état-seul | 0,00086627 | 0,00086783 |

La face interne **3** relie exactement les cellules **0 (owner)** et
**12 (neighbour)** dans la topologie BFS152 fournie, et concentre le
plus grand écart absolu des flux dans cette fenêtre. La cellule 0
correspond aussi au maximum `rAU`, tandis que la cellule 12 devient le
foyer du résidu `omega` au pas très réduit. Il s'agit d'une **zone de
correlation** à étudier et non d'une relation causale démontrée.

La différence de `nu_t` entre ancien et candidat est une évolution de
la solution ; son identité avec la version état-seul à état identique
exclut seulement, dans cette fenêtre, une divergence directe de stockage
de la viscosité.

### 2.2 Reproductibilité et tests

Deux rejouages instrumentés ont atteint l'échec intentionnel V65 à
l'itération 23. Le rejeu final fournit **149 traces d'essais V65**,
**22 transactions acceptées** et **18 sondes ciblées V66** (itérations
22 et 23) dont **6 décompositions FV signées de `omega`**. L'auditeur
V66 a vérifié **textuellement** les 149 essais V65 et les deux jeux
de 22 certificats de commit contre le journal V65 original. Le test
Python contient 13 cas, y compris la fermeture négative due au bruit
à l'arrondi de l'ordre de 1e-15, qui doit être acceptée, et la
fermeture négative significative, qui doit être rejetée.

Les **8/8 CTest ciblés** passent. Ces preuves montrent une
**non-régression du diagnostic** et la cohérence de la décomposition
sur cet essai ; elles ne prouvent ni une correction de la stagnation
non linéaire ni une convergence physique. Le code de sortie 2 de la
campagne est maintenu comme échec physique.

## 3. Distinction opérateur / préconditionneur / convergence

- Le C++ de l'oracle reconstitue le flux Rhie–Chow avec les `rAU` de chaque état, et calcule **à part** l'écart au flux porteur COUPLED.
- Une différence entre `rAU` porteur et `rAU` de `F_phys(X)` peut être attendue lors du réassemblage des coefficients. **Elle n'est pas automatiquement un bug.** Il faut relier sa localisation à un défaut de conservation ou à une erreur de transaction avant de changer l'opérateur.
- Les bilans SST signés sont des **décompositions d'un résidu FV**, pas une preuve d'une solution physique convergée ni d'une dérivée Newton.
- La comparaison V65 reste soumise aux gardes historiques : croissance maximale d'un bloc `1,20` et mérite relatif `<=0,995`. Le mode V66 n'intervient ni sur ces tests ni sur le choix de bêta.
- `selfp`, PCD, LSC, BFBt, BlockLocal, SIMPLE-Schur et SIMPLEC-Schur sont conservés sans changement. Aucun nouveau Schur, Jacobien ni Newton.

## 4. Qualification logiciel et portée de la livraison

Les nouveaux tests C++ et Python comprennent des falsifications de flux, d'indices de face/cellule, de coefficients non finis, de bilans de SST et de faux certificats physiques. L'auditeur peut demander une **parité textuelle stricte des traces V65** avec la campagne V65 archivée. Un audit de traces ne vaut pas certification de convergence.

Reproduction (configuration de diagnostic, non-performance) :

```bash
cmake -S . -B build-v66 -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v66 --target cfdx_production_solver \
    test_bfs_flux_sst_causal_v66 test_comparable_physical_merit_v65 -j2
ctest --test-dir build-v66 -R 'test_bfs_flux_sst_causal_v66|test_bfs_causal_v66_evidence|test_comparable_physical_merit_v65' --output-on-failure
CFDX_FPHYS_CANDIDATE_V51=1 CFDX_FPHYS_V65_SAME_OPERATOR=1 \
CFDX_FPHYS_V66_OMEGA_FLUX=1 CFDX_BFS_STAGE_TRACE=1 \
    ./build-v66/cfdx_production_solver \
    --mesh path/to/bfs152_v47_v65_budget40.cfdx.h5 \
    --experimental-sst --iterations 23 --output-dir build-v66/bfs-v66
python3 scripts/audit_bfs_causal_v66.py \
    --log path/to/bfs-v66.log \
    --reference-log path/to/bfs152_v65_globalisation_40_requested_failed_at23.log \
    --out build-v66/bfs-v66-audit.json
```

Le binaire V65/V66 doit échouer à l'itération 23 sur le pilote étudié ; **un code de sortie d'échec n'est pas transformé en PASS physique**. La CLI d'origine a une limite par défaut de deux itérations : passer explicitement `--iterations 23`, sans quoi la campagne ne reproduit pas le verrou !

**Limitations :** aucune convergence de BFS152 certifiée, aucun BFS712/BFS20 540, aucune validation NASA/OpenFOAM avec parois SST actives, aucun GPU/MPI, aucune CI complète, aucun push GitHub.