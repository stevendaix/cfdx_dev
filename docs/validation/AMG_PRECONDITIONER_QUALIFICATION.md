# AMG and preconditioner qualification

This document defines the quantitative qualification slice for the native linear solver/preconditioner stack.

## Scope

The current CFDX production-native methods are:

- Jacobi;
- Gauss-Seidel;
- ILU(0);
- native Boomer-style AMG;
- native smoothed-aggregation AMG;
- native FieldSplit / block-Schur for coupled systems.

The native AMG variants are dependency-free CFDX implementations. They are inspired by documented HYPRE/PETSc concepts but are not HYPRE or PETSc implementations.

## Required evidence

Every benchmark records:

- matrix/problem family;
- number of unknowns;
- Krylov method;
- preconditioner;
- solver status;
- Krylov iterations;
- reported residual;
- independently recomputed true residual;
- wall time;
- AMG hierarchy size where applicable;
- hierarchy build count and numerical-update count where applicable.

The performance fields are diagnostic measurements. A PASS is based on solver correctness and independent residual checks, not on an arbitrary timing target.

## Elliptic systems

Poisson/diffusion systems use CG. The qualification ladder compares Jacobi with native AMG on increasing one-dimensional Poisson systems.

The algebraic tolerance must be sufficiently tighter than the discretisation error when the benchmark is subsequently coupled to a finite-volume refinement campaign. Solver convergence must therefore not be used to mask spatial discretisation error.

## Transport systems

Nonsymmetric transport-like systems use BiCGStab. ILU(0) is compared with Jacobi. The benchmark is a robustness/performance diagnostic; it does not claim that ILU(0) is universally optimal.

## Hierarchy reuse

For an unchanged CSR sparsity pattern:

1. setup(A) builds the symbolic hierarchy;
2. update_values(A_new) refreshes numerical coefficients;
3. the hierarchy build count must remain unchanged;
4. the numerical-update count must increase;
5. a true residual check must still pass.

If the CSR pattern changes, the coupled Schur numeric-update operation fails explicitly and requires a new setup(). The caller owns the symbolic rebuild decision; the preconditioner must not silently rebuild or substitute another method.

## Selection policy

The problem-aware dispatcher currently selects:

- pressure/diffusion: CG + Jacobi for tiny systems, native AMG otherwise;
- momentum/scalar transport: BiCGStab + Jacobi for tiny systems, ILU(0) otherwise;
- coupled pressure-velocity: FGMRES + native block-Schur.

Explicitly requested unavailable methods remain unavailable. In particular, ILUT, FSAI, RAS, LSC and MGR are catalogue entries for future work and must not be silently mapped to another preconditioner.

## Current limits

The native AMG implementation is serial CPU infrastructure at this stage. MPI/GPU AMG, HMIS/extended interpolation and external HYPRE/PETSc backends are separate work and are not inferred from the current native implementation.

## Validation relationship

This campaign complements, rather than replaces:

- Poiseuille spatial-order validation;
- Couette pressure-velocity validation;
- Ghia cavity validation;
- MMS/operator verification;
- conservation and boundedness diagnostics.

A faster linear solve is not evidence of a more accurate CFD discretisation.


## Block Schur design reference

The coupled 4N preconditioner follows the same numerical decomposition used by
established block-preconditioning frameworks, without depending on their
implementations.

- **PETSc PCFIELDSPLIT:** CFDX now exposes the four Schur factorization forms
  diagonal, lower, upper, and full. The default is full, corresponding to the
  approximate block-LDU inverse using the velocity block and Schur solve.
  PETSc documents these variants as the standard Schur block-factorization
  choices.
- **PETSc Schur approximations:** CFDX supports the full cell-local velocity
  block in the Schur approximation C-D M_b^{-1} G, and a block-diagonal
  approximation C-D diag(M_b)^{-1} G. The latter is deliberately explicit
  rather than an implicit fallback.
- **hypre BoomerAMG/MGR:** the pressure Schur system is treated as a separate
  AMG problem, while the coupled variables retain explicit block information.
  This follows the MGR principle that strongly coupled PDE systems benefit
  from physics-aware block reduction rather than treating every unknown as an
  unrelated scalar AMG degree of freedom.
- **Future MGR path:** native MGR-style C/F reduction, MPI distributed block
  metadata, and GPU block AMG remain separate qualification work. They are not
  claimed to be implemented by this PR.

These references are design guidance only; CFDX contains no PETSc or hypre
runtime dependency in this preconditioner.


## Coupled 4N factorization contract

For A = [M G; D C] and S~ = C - D M^{-1} G, the 4N preconditioner
implements the standard block factorization actions:

- diagonal: M^{-1} r_u and S~^{-1} r_p independently;
- lower: M^{-1} r_u followed by S~^{-1}(r_p-D M^{-1}r_u);
- upper: S~^{-1} r_p followed by M^{-1}(r_u-G S~^{-1}r_p);
- full: lower Schur elimination followed by the same upper velocity correction.

The velocity inverse used by the apply phase is the exact cell-local 3x3
block inverse. The optional diagonal velocity approximation is used only when
forming the Schur approximation; this distinction is intentional and is
covered by a non-diagonal block qualification test.

A Schur CSR graph change is a symbolic invalidation event. update_values()
therefore rejects it instead of silently rebuilding the hierarchy. An
unchanged graph may use the native AMG numeric-refresh path.
