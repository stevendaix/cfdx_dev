# CFDX V69 — boucle de transport SST avec flux Rhie–Chow réévalué (expérience uniquement)

## Hypothèse falsifiable

V68 a établi sur BFS152, à état six champs **rigoureusement identique**, que
l'écart entre le flux porteur COUPLED et le flux Rhie–Chow de `F_phys` explique
pratiquement tout le défaut brut de l'équation `omega` à l'itération 23.
Ce fait ne valide ni un nouveau flux porteur de production ni un relâchement
des gardes de globalisation V13/V14/V65.

V69 utilise **les opérateurs déjà présents** (`advance_sst_production_iterate_transactionally`,
`evaluate_physical_candidate_v51`, `sst_coupled_residual_v22`, fermeture Rhie–Chow V49).
L'expérience est strictement **en lecture seule**, sous contrôle explicite
`CFDX_FPHYS_V69_SST_RC_PICARD=1` (requiert V51, V65 et V68).

## Protocole de comparaison

Sur les essais de l'itération 23 `retry=0,8,10` :

1. Copier le candidat `X=(U,p,k,omega)` post-SST et prendre `phi_RC(X)` de V51.
2. Sans changer `U,p`, résoudre le transport **réel** de `k` et `omega` avec
   **ce seul flux**. Les conditions aux limites et tolérances de transport sont
   celles de l'exécution de production. La copie du SST dérivé est réactualisée
   à partir du maillage et de `U` par la routine standard.
3. Reconstruire `nu_t`, la matrice de momentum, `rAU` puis `phi_RC` depuis
   l'état après transport, via le même `F_phys(X)` que V51.
4. Refaire une seconde fois (2 itérations SST/flux). Consigner les normes
   **brutes signées** des blocs `k,omega`, la variation de `nu_t`, la variation
   du flux, le défaut de continuité et le mérite physique complet.
5. Rejouer **bit-à-bit** `F_phys` pour chacun des deux états. En cas de flux RC
   non fermé, transport non convergé, sortie non finie, U/p modifiés, état
   turbulent hors bornes ou opérateur dépendant de l'historique, la sonde
   échoue fermée. Aucun état contrefactuel n'est engagé.

Attention : une fermeture du point fixe `rAU`/Rhie–Chow **n'est pas** un
certificat de divergence de masse nulle. V69 suit donc la divergence séparément
et ne remplace **jamais** le flux conservatif COUPLED.

## Contrôles conservés

- Aucun changement de Schur, PCD, AMG, Newton, CFL, propriété physique,
  équation SST ou contrôle de convergence.
- Les essais V65 restent soumis aux mêmes seuils historiques (croissance
  maximale par bloc 1,20 ; mérite relatif <= 0,995).
- La fonctionnalité n'est pas activée par défaut, et le modèle de paroi V17
  actif ne fait pas partie de ce premier pilote : une telle demande émet une
  sonde non qualifiée sans changer le calcul.
- Les états de production restent gouvernés par la transaction COUPLED
  et les décisions de V65, indépendamment du succès des sondes privées.

## Tests et limites de preuve

- `test_bfs_sst_flux_picard_v69` : sondes lisses, replay déterministe,
  préservation U/p et de l'état original, refus des états invalides.
- `test_bfs_sst_rc_picard_v69_evidence` : falsifications du journal et du
  certificat, vérification stricte des champs `committed=0/qualified=0`.
- `scripts/audit_bfs_v69_sst_rc_picard.py` compare les 149 décisions V65 et
  trois traces V68 à un rejeu de référence et exige six traces V69 pour
  certifier que le **protocole** a été entièrement exécuté.

La campagne n'est PAS une validation SST/Driver-Seegmiller. Elle utilise le
pilote TH03 BFS152 avec traitement de paroi historique ; BF712, BF20540,
parois OpenFOAM actives, comparaison NASA et stabilité d'une véritable
itération extérieure restent à qualifier séparément.

## Reproduction

```bash
cmake -S . -B build-v69 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v69 --target cfdx_production_solver \
  test_bfs_sst_flux_picard_v69 -j2
ctest --test-dir build-v69 --output-on-failure -R 'test_bfs_sst_flux_picard_v69|test_bfs_sst_rc_picard_v69_evidence'
CFDX_FPHYS_CANDIDATE_V51=1 CFDX_FPHYS_V65_SAME_OPERATOR=1 \
CFDX_FPHYS_V68_FLUX_SWAP=1 CFDX_FPHYS_V69_SST_RC_PICARD=1 \
CFDX_BFS_STAGE_TRACE=1 build-v69/cfdx_production_solver \
  --mesh path/to/bfs152.cfdx.h5 --experimental-sst \
  --iterations 23 --output-dir bfs-v69 2>bfs-v69.log
python3 scripts/audit_bfs_v69_sst_rc_picard.py \
  --log bfs-v69.log --baseline bfs-v68.log --output audit-v69.json
```

Le code de sortie non nul du solveur à l'itération 23 est attendu pour le
pilote historique V65 et ne doit **jamais** être compté comme une convergence.