# NACA N128/N256 — AMG C/F complexity and CSR Galerkin qualification (V4, 2026-11-05)

## Scope and version

This qualification is built **on top of V3 commit `5f1943e`** and its six captured production pressure-correction matrices. The new implementation is self-contained C++17, with no PETSc/HYPRE runtime dependency and no third-party code copied.

- N128: 12,336 pressure unknowns, 61,166 CSR entries.
- N256: 49,248 pressure unknowns, 245,214 CSR entries.
- Three successive correction-pressure matrices per mesh, with the **original** RHS and initial guess, exactly as assembled by the NACA production driver.
- Inner Krylov true-residual threshold is unchanged at `1e-10`; no row is accepted on the basis of AMG's internal residual alone.

The scope is *linear qualification and short production-trajectory regression*, not fully converged NACA aerodynamic validation.

## What changed

### 1. Limit direct C/F interpolation fill, with constant preservation

The classical C/F interpolation may produce many coarse neighbors per fine row. On the existing Boomer-style path, an **opt-in** limit `NativeAMGPreconditioner::configure_cf_interpolation_limit(4)` retains the four largest absolute interpolation entries. Deterministic ties break on coarse-column index. The retained signed weights are renormalized, maintaining `P * 1 = 1` for pressure near-nullspace compatibility. A near-zero retained sum is a **setup failure**, not a silent fallback. The first aggregate diagnostic is reconciled with the kept row.

The default cap is **0 (unlimited)**. Keeping it opt-in avoids unqualified changes for other CFDX physics/meshes. It is independent of the already supported smoothed-aggregation prolongator filter.

### 2. Symbolic/numeric CSR Galerkin product

Previously every coarse operator and every numeric refresh accumulated `P^T A P` into `std::vector<std::map<column,double>>`, then staged COO and sorted again. The new `GalerkinCSRProduct` performs:

1. **Symbolic** CSR graph construction for `A*P` and `P^T(A*P)` using sparse per-row markers, with sorted unique columns.
2. **Numeric** sparse accumulation using cached graph positions, followed by direct sorted-CSR assembly (`SparseMatrix::from_sorted_csr`).
3. **Reuse** of the symbolic structures on subsequent matrix-value-only changes, with the previous transactional AMG update lifecycle unchanged.
4. **Fail-closed** graph checks. If the structural pattern of a coarse matrix drifts (e.g., a coefficient crosses the pre-existing `1e-30` pruning threshold), the numeric refresh aborts. No old hierarchy state is mutated; callers may perform a full setup.

No non-Galerkin coarse dropping or tolerance loosening has been introduced. The code still uses `std::map` for *some interpolation weight construction*; only the expensive Galerkin accumulation was replaced.

## Measured NACA results

Paired V3 vs V4 replay, **four repetitions**, rotated run order, same six pressure captures, same residual gate:

| Mesh/method | V3 3 solves | V4 CSR (no cap) | V4 CSR + cap 4 | Operator complexity V3 → cap 4 |
|---|---:|---:|---:|---:|
| N128 FGMRES + Boomer Fast | 315.77 ms | 151.38 ms | 139.25 ms | 4.561 → 3.711 |
| N256 FGMRES + Boomer Fast | 1497.14 ms | 696.30 ms | 658.64 ms | 4.466 → 3.624 |
| N128 FGMRES + SA Fast | 191.26 ms | 195.66 ms | not applicable | 2.327 → 2.316 (CSR only) |
| N256 FGMRES + SA Fast | 1038.72 ms | 1031.31 ms | not applicable | 2.357 → 2.348 (CSR only) |

Time includes AMG setup, two numeric updates and all three FGMRES resolutions, **excluding snapshot disk loading**. For Boomer the V4 CSR uncapped version is 2.09× (N128) and 2.15× (N256) faster than V3; with cap 4 these ratios are 2.27× and 2.27× in this local campaign. Boomer iteration totals across three captures change from 29 to 32 on N128 and from 32 to 35 on N256 with cap 4. The SA path is approximately neutral; N128 shows a small **regression** (~2.3%) and is not reported as an improvement.

All **120 replayed linear solves** passed the independently computed true-residual gate. These ratios are local benchmark results, **not** universal acceleration factors and **not** complete CFD simulation speedups.

A single comparable N256 Boomer replay gave peak RSS:
- V3: 171,424 kB;
- V4 CSR uncapped: 152,380 kB;
- V4 CSR + cap 4: 143,152 kB.

This is a single-process observation; memory performance remains subject to repeatability/allocator variance.

## Physical-path regression on production NACA driver

The production NACA `--capture-pressure` pilot ran for three outer iterations on each mesh. Compared to the V3 reference history in `evidence/next/physical_final/`, the maximum **absolute** difference across the first three iterations in each of continuity, momentum RHS residual, mass imbalance, local mass Linf, and flux-velocity mismatch was:

- N128: at most **1.28e-13** (largest of those five measures).
- N256: at most **8.91e-14**.

Pressure Krylov iteration counts on the unmodified default pilot were identical to the V3 reference (N128: 15/14/15; N256: 17/17/18). Physical output was not compared for a cap-4-enabled production driver; the cap has been qualified on its original captured pressure systems and anisotropic manufactured matrices only.

## Tests, counterexamples and limits

- 17 targeted Release CTest groups passed in the final campaign (including the new benchmark fail-closed Python regression). All eight dedicated Galerkin/CF-cap C++ scenarios also passed in Release.
- The eight dedicated C++ tests passed with AddressSanitizer + UndefinedBehaviorSanitizer, with leak detection enabled.
- Three Python checks of the benchmark parser (accept valid result, reject explicit failure, reject falsely reported success) passed.
- The Galerkin tests include equivalence against the **legacy map-based formula**, a value-only update with retained pattern, two forms of structural mismatch (fine and coarse), transaction non-mutation, and `P * 1` consistency.
- The interpolation cap was also exercised on ratios `1e3` and `1e6` in both directions, a 45-degree-dominant discrete diffusion operator, and FGMRES with a genuine projected Neumann pressure nullspace.

**Not proven:** full NACA N128/N256 steady convergence, pressure/viscous force accuracy, three-mesh GCI, CPU thread scaling, MPI or GPU; runtime validation against the official HYPRE/PETSc distributions. Those are separate qualification gates.

## Reproduction

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_galerkin_csr_symbolic_numeric cfdx_replay_pressure_sequence -j 4
ctest --test-dir build -R '^(test_galerkin_csr_symbolic_numeric|test_amg_preconditioner_qualification|test_neumann_anisotropic_pressure_qualification)$' --output-on-failure

# Use local V3 and V4 binaries; captures are not synthetic stencils.
python3 tools/benchmarks/compare_naca_amg_csr.py \
  --baseline /path/to/v3/cfdx_replay_pressure_sequence \
  --candidate build/cfdx_replay_pressure_sequence \
  --captures /path/to/N128_N256_captures \
  --repetitions 4 \
  --raw /tmp/naca_csr_raw.csv --summary /tmp/naca_csr_summary.csv

# Select the experimental cap on an identical NACA replay:
build/cfdx_replay_pressure_sequence fgmres_boomer_fast_cap4 reuse 1 /path/to/N128_portable*.manifest
```

## Algorithm references

- HYPRE BoomerAMG documentation: interpolation truncation and maximum row entries, PMIS/HMIS coarsening, strong connections: https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html
- PETSc `MatProductSymbolic` / `MatProductNumeric` design: https://petsc.org/main/manualpages/Mat/MatProductSymbolic/
- PETSc `MatPtAP` reuse capability: https://petsc.org/main/manualpages/Mat/MatPtAP/

These are references for algorithmic ideas; CFDX implementations are not vendored upstream PETSc or HYPRE source code.