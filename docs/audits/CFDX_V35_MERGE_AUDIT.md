# CFDX V35 — audit de fusion sur archive BFS V6 fournie

Date : 2026-11-11. Préimage : `cfdx_dev-physics-bfs-v6-hilo-rhie-chow-production.zip`, SHA-256 dans `CFDX_V35_MERGE_PROVENANCE.json`.

## Travail réellement intégré
- Conserve les fichiers et interfaces spécifiques de la base BFS V6 (high/low Rhie-Chow, UI, comparateurs et tests GUI).
- Ajoute 67 fichiers présents seulement dans les sources V34 (diagnostics, oracles et tests SST V23-V34).
- Remplace six modules de calcul/couplage SST et JFNK par les versions V34.
- Inscrit 14 tests SST complémentaires dans CMake et supprime les cibles CTest répétées dans la base fournie.
- Supprime un include et des options CLI dupliqués dans `apps/production_solver/main.cpp`.

## Qualification effectuée
- Configuration CMake Release : PASS (MPI désactivé ; HDF5 non imposé par les options de lancement).
- Compilation de trois exécutables unitaires : PASS.
- Tests `test_sst_jfnk_v22`, `test_sst_spatial_block_preconditioner_v24`, `test_sst_branch_preconditioner_v34` : 3/3 PASS.
- **Compilation `cfdx_production_solver` : FAIL**. Le ZIP de base contient des redéclarations dans `pressure_velocity_algorithms.h` et des incohérences structurelles dans `steady_incompressible_solver.h`. La fusion seule ne les répare pas. Le log est livré.
- BFS 712/20540, six Schur, CI complète, MPI/GPU, non exécutés. Aucune convergence BFS démontrée.

## Décision
Patch de **fusion expérimentale non qualifiée pour la production**. Ne pas activer les fonctionnalités V34 par défaut ni promouvoir le code en production avant correction des erreurs du chemin BFS de base et tests end-to-end.

## Reproductibilité
Sur les octets de l'archive fournie, appliquer `git apply --check` puis `git apply` au patch cumulatif ; ne pas appliquer un deuxième patch en plus.
