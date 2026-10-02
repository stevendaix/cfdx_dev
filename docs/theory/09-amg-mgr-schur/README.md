# 09 — AMG, MGR and Schur Complements

## 1. Saddle-point structure

In incompressible flow,
\[
\mathcal A=
\begin{bmatrix}A&G\\D&0\end{bmatrix}.
\]
Eliminating velocity gives
\[
S=DA^{-1}G.
\]
A practical solver approximates \(A^{-1}\) and \(S^{-1}\) rather than forming dense inverses.

**CFDX paths:** `src/cfdx/core/linalg/exact_schur.h`, `src/cfdx/core/linalg/block_schur.h`, `src/cfdx/core/linalg/schur_approximation.h`.

## 2. Approximate Schur

A common approximation replaces \(A^{-1}\) with an inexpensive block approximation \(\tilde A^{-1}\):
\[
\tilde S=D\tilde A^{-1}G.
\]
The approximation should be judged by the spectrum/conditioning of the preconditioned coupled operator, not by its name.

## 3. AMG

Algebraic multigrid constructs a hierarchy
\[
A_0,A_1,\ldots,A_L
\]
such that smooth error is reduced on coarse levels. A Galerkin coarse operator is
\[
A_c=R A_f P.
\]
With \(R=P^T\) in a symmetric setting,
\[
A_c=P^TA_fP.
\]

**CFDX paths:** `src/cfdx/core/linalg/hypre_amg.h`, `src/cfdx/core/linalg/hypre_amg.cpp`, `src/cfdx/core/linalg/amg_preconditioner.h`.

## 4. Interpolation

Interpolation maps coarse corrections into fine space:
\[
e_f\approx Pe_c.
\]
Its graph-based construction must preserve the intended low-energy error components. Rank, sparsity and sign properties should be diagnosed.

## 5. V-cycle

A conceptual V-cycle is
\[
x\leftarrow S_{pre}(A,b,x),
\]
\[
r=b-Ax,
\]
\[
r_c=Rr,
\]
\[
A_ce_c=r_c,
\]
\[
x\leftarrow x+Pe_c,
\]
followed by post-smoothing.

## 6. Energy contraction

For SPD problems, a useful measure is
\[
\rho_E=\frac{\|e_{out}\|_A}{\|e_{in}\|_A}.
\]
A qualification gate should specify the tested matrix family and the expected contraction; a residual-only gate is insufficient when the contract is energy-based.

## 7. MGR

Multigrid reduction generalises coarse-grid reduction for coupled/block systems. Variables may be reduced in stages, e.g.
\[
4N\rightarrow3N_{velocity}\rightarrow N_{Schur}.
\]
The strength of coupling and the transfer operators determine whether the reduction is effective.

**CFDX path:** `src/cfdx/core/linalg/` and MGR-specific implementations/tests.

## 8. Update semantics

If matrix values change while topology remains fixed, a preconditioner may require value refresh without rebuilding graph structure. The distinction must be explicit to avoid stale coefficients.

## 9. Verification

Use exact Schur oracles, known SPD matrices, Galerkin identity checks, energy contraction, block-coupled systems and update-values tests.

**Evidence:** `docs/development/NATIVE_AMG_FIELDSPLIT.md`, `docs/validation/AMG_PRECONDITIONER_QUALIFICATION.md`.
