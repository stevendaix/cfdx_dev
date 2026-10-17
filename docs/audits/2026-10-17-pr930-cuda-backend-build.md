# Fiche — CUDA backend build enablement (PR #930)

## Contexte

Le backend CUDA de CFDX (chemin `src/cfdx/runtime/gpu/`) n'était jamais compilé
localement : `CFDX_ENABLE_GPU` était désactivé dans tous les builds existants
(`build-n9pr836` etc.). Les tests GPU (`test_cuda_runtime`,
`test_n13_cuda_equivalence`, `test_phase8_cuda_poisson`) n'étaient donc jamais
exécutés, et le code CUDA contenait des bugs que cette compilation révèle.

## Problème racine

Lors de l'activation de `CFDX_ENABLE_GPU=ON` avec CUDA 11.5 et GCC 11 (libstdc++
11.4), `nvcc` échouait à compiler `cuda_poisson.cu` :

- `std_function.h:435` — `_Functor&& __f` non expansé (`parameter packs not
  expanded with '...'`). Cela vient de `std::function` dans
  `cfdx/core/linalg/linear_operator.h`.
- `mpi_utils.h:154` — `if constexpr (std::is_same_v<T, ...>)` indisponible.

Ces erreurs surviennent parce que `cuda_poisson.cu` -> `cuda_poisson.h` ->
`cg_solver.h` -> `preconditioner.h` (structured bindings C++17) ->
`krylov_reductions.h` -> `mpi_utils.h` (`if constexpr`/`is_same_v`). Le toolkit
CUDA 11.5 + libstdc++ GCC 11 ne supporte pas ces fonctionnalités C++17 en mode
device.

## Correction

1. **`CMakeLists.txt`** : `CMAKE_CUDA_STANDARD` 17 -> 14. Aucune fonctionnalité
   C++17 n'est utilisée sur le chemin device CUDA.

2. **`src/cfdx/core/linalg/solver_status.h`** (nouveau) : header léger,
   compatible C++14, contenant `SolverStatus`, `SolverResult`,
   `SolverDiagnosticSample` et `to_string`. `cuda_poisson.cu` inclut ce header
   au lieu de `cg_solver.h`, coupant la chaîne C++17.

3. **`src/cfdx/core/linalg/cg_solver.h`** : inclut `solver_status.h`,
   déduplique les types.

4. **`src/cfdx/runtime/gpu/cuda_device_buffer.h`** : ajoute les opérateurs de
   déplacement. Les copies étant supprimées, `std::vector<CudaDeviceBuffer>::resize()`
   échouait dans `CudaFieldMirror::allocate()`.

5. **`src/cfdx/runtime/gpu/gpu_execution.h`** : qualifie
   `CudaPinnedBufferPool` (namespace `cfdx::runtime::ooc`).

6. **`tests/validation/test_n13_cuda_equivalence.cpp`** : ajoute
   `using namespace cfdx::runtime::gpu`. Le test est désormais conditionné à
   `CFDX_ENABLE_GPU` dans `CMakeLists.txt` (plus de `SKIP_RETURN_CODE`).

## Vérification

- Build GPU (`build-gpu`, CUDA 11.5 + RTX 4070) : compile et 11/11 tests GPU+MPI
  passent.
- 54 tests unitaires noyau/algèbre linéaire passent.
- Build non-GPU (`build-n9pr836`) : compile toujours, `test_n13_cuda_equivalence`
  n'est plus enregistré.
- Pre-commit, audit gate (fiche présente), CFDX CI, CFDX GUI CI : tous passent.

## Limites

- `CMAKE_CUDA_STANDARD=14` est imposé par la toolchain locale (CUDA 11.5 +
  GCC 11). Les builds avec CUDA >= 12 et/ou libstdc++ récent devraient pouvoir
  revenir à C++17 si le code le justifie.
- Les tests de performance GPU (N15) ne sont pas couverts par cette PR.