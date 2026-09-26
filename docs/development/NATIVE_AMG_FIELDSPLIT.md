# Native AMG and field splitting

## Decision

CFDX does not vendor HYPRE or PETSc source code. The native path implements the
small set of documented algebraic methods needed by the current solvers and
keeps the existing `SparseMatrix`, `Vector` and `Preconditioner` contracts.

The names deliberately say `Native`: these classes are not drop-in HYPRE or
PETSc backends and must not be presented as such.

## Boomer-style AMG recipe

`NativeBoomerAMGPreconditioner` reuses the CFDX Galerkin V-cycle. Its setup
pipeline follows the same high-level order as BoomerAMG:

- a signed strength-of-connection graph with threshold 0.25 and the default
  0.9 row-sum guard;
- a deterministic serial C/F maximal-independent-set split;
- normalized direct interpolation from C points to F points;
- transpose interpolation for restriction;
- Galerkin coarse operators;
- equal weighted-Jacobi sweeps before and after coarse correction, preserving
  the symmetric application required by CG;
- a pivoted direct solve on the smallest level; pure-Neumann systems use a
  symmetric zero-mean saddle-point augmentation instead of pinning a row.

The policy changes sweep count and maximum hierarchy depth. The serial C/F
split captures the independent-set completion principle, but it is not the
distributed Ruge + PMIS sequence used by HMIS. Likewise, normalized direct
interpolation is intentionally smaller than extended+i interpolation.

Reference: [HYPRE BoomerAMG documentation](https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html).

Source paths audited against HYPRE commit
`5b582610a113ede9162a5744d33873b860498a11`:

- `src/parcsr_ls/par_amg.c`: defaults and solver lifecycle;
- `src/parcsr_ls/par_amg_setup.c`: strength, coarsening, interpolation and
  `RAP` setup order;
- `src/parcsr_ls/par_strength.c`: signed strength and row-sum criteria;
- `src/parcsr_ls/par_coarsen.c`: HMIS as Ruge coarsening plus PMIS completion;
- `src/parcsr_ls/par_interp.c`: direct interpolation normalization;
- `src/parcsr_ls/par_cycle.c`: smooth, residual, restrict, coarse solve,
  prolong and post-smooth order.

## Smoothed-aggregation AMG recipe

`NativeSmoothedAggregationAMGPreconditioner` is the second selectable native
AMG family. It follows the transferable setup concepts of PETSc GAMG without
using PETSc source or runtime objects:

- deterministic uncoupled aggregates built from the signed strength graph;
- a piecewise-constant tentative prolongator `P0`;
- one damped-Jacobi interpolation step
  `P = (I - (2/3) D^-1 A) P0`;
- row normalization that preserves the constant candidate mode;
- transpose restriction, Galerkin coarse operators, symmetric weighted-Jacobi
  smoothing and the same pivoted coarse solve as the Boomer-style path.

The model is available as `PreconditionerModel::SmoothedAggregationAMG` and is
qualified with CG for pressure-Poisson and diffusion operators. Explicit
selection is propagated through `IncompressibleSolverControls`, so every
segregated pressure-velocity algorithm can choose it. The hierarchy supports
the same unchanged-operator reuse and same-pattern numeric refresh as
`NativeAMG`.

Current near-null-space support is limited to the built-in constant mode. For a
pure-Neumann pressure/diffusion operator, the dispatcher passes that mode into
either native AMG family. Every residual, smoother correction and coarse right-
hand side is projected; prolongation reproduces constants; and the smallest
singular operator is solved with a symmetric zero-mean constraint. Automatic
selection therefore uses projected CG/AMG for systems of at least 24 equations
and retains projected CG/Jacobi for smaller systems. Arbitrary user vectors,
MPI aggregation and GPU kernels remain separate work.

Reference: [PETSc GAMG manual page](https://petsc.org/main/manualpages/PC/PCGAMG/).

## FieldSplit-style recipe

`NativeFieldSplitPreconditioner` accepts two explicit, disjoint index sets that
must cover the matrix. It supports:

- additive block Jacobi;
- diagonal Schur with the conventional `-S` sign;
- lower and upper Schur triangular factors;
- full Schur factorization.

For a matrix `[A B; C E]`, the implementation forms
`S = E - C A^-1 B`. The current field inverses and Schur complement are dense,
so this is a correctness baseline for small coupled systems, not a scalable
replacement for PETSc. Larger systems should replace the inner dense solves
through the existing CFDX preconditioner interface.

References: [PETSc field splitting manual](https://petsc.org/main/manual/ksp/) and
[PETSc Schur factorization API](https://petsc.org/main/manualpages/PC/PCFieldSplitSetSchurFactType/).

The application order was also checked against PETSc commit
`74bbc03299ffb37f3e3b88d943fcbbe2cdac1b1b`, file
`src/ksp/pc/impls/fieldsplit/fieldsplit.c`: diagonal applies
`diag(A^-1, -S^-1)`, lower and upper retain the corresponding triangular
factor, and full performs the two first-field solves required by `LDU`.

## Solver integration

GMRES and BiCGStab already consume `Preconditioner`. CG now has an overload that
accepts the same interface; the original overload still selects its historical
Jacobi path. All convergence tests use an independently recomputed true
residual.

`ReusableCgContext` and `ReusableGmresContext` split the preconditioner lifecycle
into full symbolic setup, same-pattern numeric refresh and unchanged-operator
reuse. The segregated incompressible solver owns one CG context for its whole
run, so SIMPLE, SIMPLEC, PISO, PIMPLE and fractional-step pressure corrections
do not rebuild the AMG hierarchy on every solve. A changed CSR pattern triggers
a full setup; changed values with the same pattern refresh level coefficients;
an unchanged matrix reuses the prepared hierarchy directly. The counters in
`IncompressibleSolveResult::pressure_linear_context` make this behavior
observable in tests and performance measurements.

The coupled pressure-velocity path does not use this scalar context. Its Schur
preconditioner needs a separate reusable block lifecycle before reuse can be
enabled without weakening its matrix consistency.

## Limits

- no HYPRE or PETSc runtime, ABI, option database or object lifecycle;
- no claim of numerical identity with either project;
- no distributed AMG hierarchy in this dependency-free implementation;
- no arbitrary near-null-space vectors (the constant pressure mode is supported);
- no GPU AMG setup or apply kernels yet;
- no performance claim without a separate reproducible benchmark.
