# Native AMG and field splitting

## Decision

CFDX does not vendor HYPRE or PETSc source code. The native path implements the
small set of documented algebraic methods needed by the current solvers and
keeps the existing `SparseMatrix`, `Vector` and `Preconditioner` contracts.

The names deliberately say `Native`: these classes are not drop-in HYPRE or
PETSc backends and must not be presented as such.

## Boomer-style AMG recipe

`NativeBoomerAMGPreconditioner` reuses the CFDX Galerkin V-cycle and adds a
strength-of-connection threshold. Its hierarchy uses:

- strong-neighbour pair aggregation with a default threshold of 0.25;
- piecewise-constant prolongation and its transpose for restriction;
- Galerkin coarse operators;
- equal weighted-Jacobi sweeps before and after coarse correction, preserving
  the symmetric application required by CG;
- a pivoted direct solve on the smallest level.

The policy changes sweep count and maximum hierarchy depth. It does not claim
to reproduce HMIS coarsening or extended+i interpolation. Those algorithms are
substantially more involved than the dependency-free vertical implemented here.

Reference: [HYPRE BoomerAMG documentation](https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html).

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

## Solver integration

GMRES and BiCGStab already consume `Preconditioner`. CG now has an overload that
accepts the same interface; the original overload still selects its historical
Jacobi path. All convergence tests use an independently recomputed true
residual.

## Limits

- no HYPRE or PETSc runtime, ABI, option database or object lifecycle;
- no claim of numerical identity with either project;
- no distributed AMG hierarchy in this dependency-free implementation;
- no performance claim without a separate reproducible benchmark.
