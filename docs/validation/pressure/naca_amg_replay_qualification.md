# NACA N128/N256 — FGMRES/AMG memory refresh and replay qualification (2026-11-05)

## Scope / source of evidence

The cases are the actual **production correction-pressure systems** captured
at the NACA0012 linear solve call, not manufactured 5-point stencils:

| Case | pressure unknowns | matrix CSR entries | tested correctors |
|---|---:|---:|---:|
| N128 | 12,336 | 61,166 | first corrector of three successive outer iterations |
| N256 | 49,248 | 245,214 | first corrector of three successive outer iterations |

The capture manifest provides the exact matrix, right-hand side, initial
pressure iterate, reference-pressure metadata and nonlinear iteration.
Capture with the existing `test_naca0012_qualification --capture-pressure`
entrypoint; this is **not** a substitute for fully converged NACA flow.

The replay always recomputes the FP64 true residual `||b-A*x||/||b||` on the
original matrix, with a required threshold of `1e-10`. A failed setup or
uncertified solve must return nonzero, even if CSV output was emitted.

## Changes

1. `MatrixFreeVcyclePreconditioner::update_values` stages **numeric-only**
   Galerkin operators, diagonals and algebraic line factorizations. It leaves
   the existing coarse-grid interpolation and aggregation untouched until a
   successful commit. Rejected rotated anisotropy leaves the old V-cycle
   bitwise unchanged.
2. `AlgebraicLinePreconditioner::apply` performs its Thomas forward and
   backward substitutions in the correction vector, eliminating a heap
   allocation at each line-smoother call. Regression covers separate and
   aliased RHS/output on both anisotropy directions.
3. The new AMG hierarchy diagnostics expose the number of levels, operator
   complexity `sum_l nnz(A_l)/nnz(A_0)` and grid complexity
   `sum_l N_l/N_0`. These are emitted in the sequence replay CSV.
4. Both replay executables fail closed on explicit linear qualification
   failure. The paired comparison script handles *different baseline and
   candidate CSV columns* while preserving all measurements.

## Run in a local build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cfdx_replay_pressure_sequence cfdx_replay_naca_pressure test_reusable_fgmres test_petsc_hypre_native_port_anisotropy
ctest --test-dir build -R '^(test_reusable_fgmres|test_petsc_hypre_native_port_anisotropy|test_naca_pressure_replay_fail_closed)$' --output-on-failure

# Capture examples (mesh files and output paths are environment specific):
CFDX_NACA_PRESSURE_SOLVER=fgmres_sa build/test_naca0012_qualification \
  --capture-pressure meshes/naca0012_laminar_n128.h5 /tmp/N128_portable 3

# Replay exactly the same three N128 physical systems:
build/cfdx_replay_pressure_sequence fgmres_sa_fast reuse 1 /tmp/N128_portable*.manifest

# Two executable versions, paired in alternating order:
python3 tools/benchmarks/compare_naca_amg_refresh.py \
  --baseline /path/to/old/cfdx_replay_pressure_sequence \
  --candidate build/cfdx_replay_pressure_sequence \
  --captures /path/to/six/capture/manifests \
  --raw /tmp/paired_raw.csv --summary /tmp/paired_summary.csv \
  --repetitions 3
```

## Measured results / limits

The local targeted CTest subset passed **16/16** (Release). The new
numeric-reuse path passed the targeted ASan/UBSan unit test. Do not infer
that the full source tree, full sanitizers suite, or an MPI build passed.

The paired campaign had **72/72 certified linear solves**, on three physical
snapshots per mesh, with two AMG variants and three alternating repetitions.
Every true relative residual was <= `1e-10`. Median replay time results from
this run were:

| Mesh and preconditioner | Baseline, three solves | Candidate, three solves | Difference |
|---|---:|---:|---:|
| N128 FGMRES + smoothed-aggregation Fast | 219.00 ms | 207.11 ms | -5.4% |
| N128 FGMRES + Boomer-style Fast | 336.99 ms | 335.85 ms | -0.3% |
| N256 FGMRES + smoothed-aggregation Fast | 1004.17 ms | 954.02 ms | -5.0% |
| N256 FGMRES + Boomer-style Fast | 1553.80 ms | 1492.99 ms | -3.9% |

**Do not interpret these as universal speedups.** Repeated campaigns showed
noticeable timing variability and sometimes slower candidates. On one
single-shot N256 SA replay, peak RSS was 153,816 kB baseline and 133,260 kB
candidate (13.4% less); wall time was 1.28 s versus 1.27 s.

The physical NACA pilot compared the first **three outer iterations** on
both meshes after the final code changes. All five logged metrics had zero
printed-precision difference versus the recorded reference histories:
continuity, momentum RHS residual, mass imbalance, local mass Linf and
flux-velocity mismatch. This is **not** full NACA convergence or validation
of Cd, Cl, Cp, Cf, or mesh GCI.

The AMG hierarchy diagnostics for the *first* real captured pressure matrix:

| Mesh | AMG method | levels | operator complexity | grid complexity |
|---|---|---:|---:|---:|
| N128 | smoothed aggregation | 7 | 2.326 | 1.452 |
| N128 | Boomer-style | 11 | 4.561 | 1.948 |
| N256 | smoothed aggregation | 8 | 2.356 | 1.451 |
| N256 | Boomer-style | 11 | 4.466 | 1.911 |

The Boomer-style variant cuts the Krylov iterations but pays roughly twice
the coarse-operator complexity. This is the next performance bottleneck.

## Follow-up implementation priorities (NOT yet qualified)

1. *P0, coarse-grid complexity:* alternative HMIS/PMIS-style strength
   coarsening, distance-two interpolation, bounded interpolation stencil and
   controlled Galerkin coarse sparsification, preserving near-null modes;
   compare true residuals and whole-CFD outcomes.
2. *P0, symbolic/numeric Galerkin split:* the current `galerkin_coarse`
   builds each coarse row in `std::map`. Evaluate a CSR two-phase SpGEMM-like
   accumulator with symbolic reuse across pressure matrix refreshes, and
   check numerical equivalence and setup-time/memory improvements.
3. *P1, adaptive Krylov recycling:* compare restart lengths and recycling
   methods (GCRO-DR/FGMRES-compatible) on successive **real** pressure
   systems, with fail-closed stagnation checks.
4. *P1, smoothers:* compare l1-Jacobi, Chebyshev and line blocks, including
   coarse-grid residual modes, symmetric CG requirements and stretched cells.
5. *P0 physical qualification:* run NACA N128 and N256 to full nonlinear
   convergence; check physical fluxes, conservation, Cp/Cf/Cd/Cl/Cm, grid
   sensitivity and time to converged solution. Add actual SST energy/radiation
   coupled campaigns once their full paths are stable.
6. *P1 parallel scalability:* benchmark OpenMP/MPI reductions and memory
   locality. No MPI speedup is currently claimed.

External references: PETSc `PCGAMGSetReuseInterpolation`, PETSc `KSPFGMRES`,
HYPRE BoomerAMG users' guide. CFDX's native methods are independently
implemented and **do not include** PETSc or HYPRE as installed dependencies.

- https://petsc.org/release/manualpages/PC/PCGAMGSetReuseInterpolation/
- https://petsc.org/release/manualpages/KSP/KSPFGMRES/
- https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html