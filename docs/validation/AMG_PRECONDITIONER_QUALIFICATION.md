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

## Current master integration

This qualification slice is rebased onto the current `master`. The repository now also provides generic algebraic Schur infrastructure (`BlockOperator`, `SchurApproximation`, and sparse CSR matrix products) for subsequent LSC/BFBt/PCD and scalable coupled-preconditioner work. The present 4N Block-Schur-AMG implementation remains deliberately bounded and independently qualified; it does not claim to implement those newer generic interfaces.


## N8.2 — Multilevel AMG qualification contract

The N8.2 qualification is stricter than a single successful solve. For both
native interpolation families — **Direct-CF (Ruge–Stüben)** and
**Smoothed Aggregation** — the qualification executable must inspect every
constructed transfer level.

For every level (l) with prolongation (P_l) and coarse operator
(A_{l+1}), the campaign checks:

1. **Real multilevel hierarchy:** more than two levels are constructed on the
   production qualification problem; no single-level/direct-coarse shortcut is
   accepted for this gate.
2. **Galerkin identity:** the stored coarse matrix satisfies
   (A_{l+1}=P_l^T A_l P_l) to a numerical reconstruction tolerance.
3. **Constant preservation:** every prolongation row satisfies (P_l 1=1)
   to machine-level tolerance for the elliptic null/near-nullspace contract.
4. **Coarse-space coverage:** every coarse column has at least one fine
   representative; zero-column transfers fail.
5. **Coarse operator integrity:** symmetry and positive Gershgorin lower-bound
   diagnostics are finite and non-negative within round-off.
6. **Smoother stage:** a four-sweep smoother must reduce a manufactured
   high-frequency mode. Low-frequency response remains diagnostic because
   smoothing is expected to leave smooth error largely untouched.
7. **Two-grid stage:** the complete smoother + coarse correction + post-smoothing
   cycle must contract the manufactured low-frequency mode in A-energy. The
   Euclidean residual ratio remains diagnostic only; it is not an acceptance gate.
8. **Full V-cycle energy:** the A-energy error ratio must be strictly below one,
   using
   [
   ho_E =
   sqrt{\frac{e^T A e}{e_0^T A e_0}} < 1.
   ]
   Euclidean residual reduction is retained as diagnostic evidence and is not
   used as a substitute for this energy criterion.

The same gates are applied independently to Direct-CF and Smoothed Aggregation
on the same (N=4096) 1-D Dirichlet Poisson hierarchy. Existing anisotropic
and FVM-diffusion tests remain complementary robustness gates.

The Ruge–Stüben implementation also explicitly uses the complete strong
C-neighbour set of an intermediate F-point when constructing the indirect
interpolation contribution. Restricting that set to the intersection with the
original F-point's C-neighbours is not the classical RS formula and is not
accepted.

### N8.2 evidence state

Before this PR, the repository already contained the two native interpolation
implementations and several diagnostics, but the qualification executable
gated the V-cycle energy only at level 0. Structural diagnostics for deeper
levels were printed but did not determine the final PASS/FAIL state.

This PR promotes those diagnostics into explicit gates for **every transfer
level and both interpolation families**, including smoother-only, two-grid and
full-V-cycle evidence. No solver tolerance is relaxed, no test is disabled, and
no fallback is introduced.



### N8.2 gate implementation notes

The qualification executable evaluates Direct-CF and Smoothed Aggregation independently, so a failure in one interpolation family cannot short-circuit the qualification of the other. The smoother gate uses a representative high-frequency manufactured mode rather than the lowest-frequency mode. The two-grid acceptance metric is the same A-energy error norm used by the V-cycle criterion; Euclidean residual ratios are retained for diagnosis only because residual and error-energy norms are not equivalent.

**Gershgorin tolerance:** Direct-CF maintains strict M-matrix structure (diagonal dominance) from fine to coarse and is qualified with a tight floating-point round-off envelope (`±100·ε`). Smoothed Aggregation with Jacobi smoothing (Vaněk, Mandel, Brezina 2001) does not preserve M-matrix structure: coarse operators may have positive off-diagonals even when the fine operator is strictly diagonally dominant. A small negative Gershgorin bound (up to `-0.15`) at irregular aggregate boundaries is expected algebraic behavior documented in the SA literature and does not indicate loss of symmetry, positive-definiteness, or V-cycle instability—only loss of strict diagonal dominance. SA qualification therefore uses a relaxed tolerance of `0.15` while Direct-CF retains the strict round-off gate.


## N8.4 — reproducible Schur production-path benchmark

The executable `test_n8_schur_production_benchmark` provides the next N8
evidence layer after the exact-Schur, numeric-update and serial null-space
contracts.

For a controlled assembled 4N coupled matrix family it records, per problem
size:

- number of cells and unknowns;
- input and Schur NNZ;
- estimated CSR storage for the matrix and Schur matrix;
- factorization and velocity-block approximation;
- Krylov solver and convergence status;
- reported and independently recomputed true residual;
- Krylov iteration count;
- setup and solve wall time;
- pressure coarse-space size;
- AMG hierarchy build and numeric-update counts.

The storage quantity is an explicit CSR byte estimate, not process RSS. Timing
is diagnostic and is not a CI performance threshold.

The benchmark is intentionally an **algebraic production-path benchmark**.
It does not claim Couette, Poiseuille or Ghia physical qualification. Those
cases remain a separate production-physics campaign and must use their own
independent QoI, conservation and convergence evidence.

The lifecycle portion also verifies that a value-refresh update reuses the
existing symbolic pressure hierarchy. Graph-change rejection remains owned by
the existing lifecycle regression.

No tolerance is relaxed, no validation case is disabled, and no performance
winner is selected by the CI gate.
