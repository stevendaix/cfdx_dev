# Fiche — Activation du build backend CUDA (PR #930)

## 1. Contexte

Le backend CUDA de CFDX (chemin `src/cfdx/runtime/gpu/`) n'était jamais compilé
localement : `CFDX_ENABLE_GPU` était désactivé dans tous les builds existants
(`build-n9pr836` etc.). Les tests GPU (`test_cuda_runtime`,
`test_n13_cuda_equivalence`, `test_phase8_cuda_poisson`) n'étaient donc jamais
exécutés, et le code CUDA contenait des bugs que cette compilation révèle.

## 2. Formulation continue

Le code CUDA est une réécritite d'opérateurs numériques (gradient de Gauss,
divergence, solveur CG de Poisson) déjà existants en C++ hôte. Aucune
formulation mathématique n'est modifiée : les noyaux CUDA implémentent la même
discrétisation volume fini que les routines CPU de référence, évaluées dans
`test_n13_cuda_equivalence` avec une tolérance relative L2/Linf <= 1e-12.

La chaîne de dépendances compromise est :
`cuda_poisson.cu` -> `cuda_poisson.h` -> `cg_solver.h` -> `preconditioner.h`
(structured bindings C++17) -> `krylov_reductions.h` -> `mpi_utils.h`
(`if constexpr` / `is_same_v`). Le toolkit CUDA 11.5 avec libstdc++ GCC 11
ne supporte pas ces fonctionnalités C++17 en mode device.

## 3. Discrétisation Volume Fini

Aucune discrétisation n'est modifiée. Les noyaux CUDA (CSR SPMV, axpy, Jacobi,
dot-product réduit) opèrent sur les mêmes structures CSR que le solveur CPU.

## 4. Traitement des frontières

Aucun traitement de frontière n'est modifié. Les tests MPI (`test_mpi_partition`,
`test_phase8_mpi_poisson`, `test_n13_mpi_equivalence`, `test_n8_mpi_boundary`,
`test_n13_restart_equivalence`) vérifient que l'échange de halo et la
reconstruction globale sont inchangés.

## 5. Pièges connus

- `CMAKE_CUDA_STANDARD=14` est imposé par la toolchain locale (CUDA 11.5 +
  GCC 11). Les builds avec CUDA >= 12 et/ou libstdc++ récent devraient pouvoir
  revenir à C++17 si le code le justifie.
- Les tests de performance GPU (N15) ne sont pas couverts par cette PR.
- `CudaDeviceBuffer` a besoin d'opérateurs de déplacement explicites car ses
  copies sont supprimées ; `std::vector<CudaDeviceBuffer>::resize()` échouait
  sans cela.

## 6. Références bibliographiques

- NVIDIA, CUDA Toolkit Documentation (cuSPARSE CSR format, launched kernels),
  accessed 2026-10-17.
- ISO/IEC/IEEE 12207:2017 — Software life cycle processes (verification,
  validation, and testing), referenced by CFDX V&V governance
  (`docs/validation/CFDX_VV_GOVERNANCE.md`).

## 7. Statut d'implémentation

Statut: **implémenté et vérifié**.

- Build GPU (`build-gpu`, CUDA 11.5 + RTX 4070) : compile et 11/11 tests
  GPU+MPI passent.
- 54 tests unitaires noyau/algèbre linéaire passent.
- Build non-GPU (`build-n9pr836`) : compile toujours, `test_n13_cuda_equivalence`
  n'est plus enregistré.
- Pre-commit, audit gate, CFDX CI, CFDX GUI CI : tous passent.

## Statut

- Implémenté : oui
- Testé : oui (11 tests GPU+MPI, 54 tests unitaires)
- Vérifié : oui (équivalence CPU/CUDA L2/Linf <= 1e-12)
- Validé : non (hors portée de cette PR)
- Qualifié : non (hors portée de cette PR)