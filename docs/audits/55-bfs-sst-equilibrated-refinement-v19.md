# CFDX BFS V19 — vrai résidu SST, pression compensée et coins OpenFOAM

Date : 2026-11-10. Source : V18 exacte issue du ZIP CFDX original `cfdx_dev-master (2)(1).zip`.

## Résultats vérifiés / limites

- **V19 fonctionne pour le premier transport SST sur le BFS de diagnostic à 712 cellules.** V18 était rejeté à 226 itérations : résidu brut omega `6.438040145440e-8`, défaut solution normalisé par ligne `1.585159911771e-6`. V19, mêmes données physiques et mêmes seuils, satisfait les critères à 217 itérations : résidu omega `3.221884981031e-10`, défaut solution `7.932853472280e-9`, premier état **committed** (`U=0.5251832153953`, continuité `0.001771416233814`, état global NON convergé).
- **La deuxième itération n'est PAS fermée** : l'action COUPLED extérieure atteint son critère vrai, mais SST revient à la limite de 800 itérations sans convergence (`k_raw~0.00780424`, `omega_raw~396.4`), rollback. Une seule itération externe acceptée ; la preuve de contraction sur trois états échoue, code retour 3.
- **V18 OpenFOAM owner-cells reste expérimentale**. Le mode V19 `v19_omega_corner_equal_faces_experimental` remplace la moyenne surfacique aux coins par la moyenne arithmétique des faces incidentes (conforme à `createAveragingWeights()` OpenFOAM Foundation pour cellules multi-faces). Au premier état, la différence maximale sur omega dans la cellule de paroi 208 est `1.659113485302e6 s^-1`, alors que le champ de vitesse est identique. Cette correction ne reproduit pas encore tous les choix `updateCoeffs` / `manipulateMatrix` OpenFOAM, le blending ou ses actualisations temporelles. Le mode est **opt-in**.
- **20 540 cellules : la pression peut être refermée à la première correction sans changer 1e-11**, avec le chemin CFDX déjà existant `v18_sparse_flux_pressure_jacobian=true`, `v19_full_pressure_flux_step=true`, `pressure_high_low_production=true` : résidu CSR standard `3.08558533495e-10`, puis vrai résidu des flux Rhie–Chow après raffinage compensé `7.71613300214e-18`, erreur de flux physique `1.06780065179e-13`, continuité normalisée `2.39525079963e-16`. Exécution arrêtée dans le SST après le premier solveur pression : pas d'état externe certifié. **Ceci ne valide pas les six Schur sur ce maillage.**
- **Justesse physique distincte** : l'état V19 accepté du BFS réduit possède encore `k_max=534.3988710451 m²/s²` en cellule 426, rapport estimé `P_k/D_k=3.6736`, différence locale `Cp(p*)-Cp(p_stat)=0.3647201`. Cet état n'a aucune qualification Driver–Seegmiller ; pas de `x_r/h`, `Cf`, `Cp` ou profils NASA déclarés bons.

## Nouvelles méthodes et contrat d'acceptation

### Correction équilibrée de l'équation omega sur l'opérateur original

La très forte valeur de l'omega de paroi biaise l'arrêt relatif du solveur linéaire, alors que l'équation dans une cellule intérieure à petite échelle conserve un défaut absolu non acceptable. Sur l'équation gelée originalement assemblée `A omega = b` :

1. Calculer chaque `r_i=b_i-Σ_j A_ij omega_j` en `long double`, `d_i=Σ_j |A_ij| > 0`.
2. Résoudre en correction `D^-1 A delta = D^-1 r` sans modifier `A` ni `b` en production.
3. Rejeter toute valeur non finie, tout candidat non positif ou toute correction qui n'améliore pas d'au moins 10 % le **défaut original** `sst_equation_solution_defect(original, omega)`.
4. Garder le seuil original du transport `1e-8`, l'exigence de mise à jour des champs et l'équation non linéaire réassemblée ; ne pas substituer le résidu auxiliaire au certificat global.
5. Option **désactivée par défaut** : `v19_sst_equilibrated_defect_correction_experimental=enabled` dans `numerics.raw_settings`.

L'action est limitée aux itérations où `max_abs(omega_solution - omega_old) <= 1e-5`, et ne constitue pas une régularisation physique ni un plafonnement turbulent.

### Coins de paroi

La moyenne OpenFOAM Foundation `omegaWallFunction::createAveragingWeights()` est liée au nombre de faces adjacentes, pas nécessairement à leurs aires. V19 ajoute la moyenne par face au triplet omega/G/nu_t de l'état commun, en préservant par défaut la moyenne surfacique V18. Activer seulement avec les conditions `kqRWallFunction`, `omegaWallFunction`, paroi immobile, l'opérateur V17 et la contrainte V18 déjà actifs.

## Tests vérifiables

- CTest : `test_sst_wall_operator_v17`, `test_sst_wall_owner_v18`, `test_sst_equilibrated_v19`.
- Python : `tests/python/test_bfs_v19_sst_precision.py`, `test_th03_bfs_case_generator.py`, `test_bfs_v10_diagnostics.py`.
- Rejeu : `cfdx_production_solver --mesh matched_b.cfdx.h5 --iterations 3 --experimental-sst`, avec `CFDX_SST_DIAGNOSTICS=1 CFDX_BFS_STAGE_TRACE=1`.
- Certificat externe : `python scripts/qualify_bfs_committed_v13.py --log matched_b3.log --output committed.json` ; contraction : `python scripts/qualify_bfs_global_contraction_v14.py --log matched_b3.log --output contraction.json` (**code 3 attendu**).
- Le ZIP d'origine et les copies source V18 servent de bases exactes pour les patchs Git cumulatif et incrémental.

## Bloquants après V19

P0 — Faire converger **la deuxième itération SST**, où la mise à jour des champs oscille malgré un COUPLED linéaire convergé. Vérifier correction omega près de la cellule 77 et près des propriétaires de paroi, imposer un pas physiquement acceptable avec résidu non linéaire totalement réévalué, puis obtenir >=3 états externes acceptés et contractants.

P0 — Répéter sur le grand BFS avec maillage officiel, sans confondre correction SIMPLE Rhie–Chow et préconditionnement des six Schur de COUPLED ; vrai résidu, coût, budget et conservation identiques entre variantes.

P1 — Reproduire toutes les règles des fonctions de paroi OpenFOAM, les profil d'entrée et son Reynolds de référence ; qualifier `p_static=p*-(2/3)rho*k`, `Cf`, `Cp`, `x_r/h`, vitesses, tensions de Reynolds et GCI 3 grilles. **L'état physique reste NON_QUALIFIED.**

Références : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html ; https://tmbwg.github.io/turbmodels/backstep_val.html
