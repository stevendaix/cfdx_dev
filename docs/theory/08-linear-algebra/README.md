# 08 — Linear Algebra

## 1. Discrete problem

Each implicit CFD equation becomes
\[
A x=b.
\]
The residual is
\[
r=b-Ax.
\]
A solver tolerance is meaningful only together with scaling, the initial residual and the true residual recomputed independently.

## 2. Sparse structure

Finite-volume matrices are sparse because a cell interacts with its local stencil. CSR-like storage represents row pointers, column indices and non-zero values without storing zeros.

## 3. Conditioning

For a nonsingular matrix,
\[
\kappa(A)=\|A\|\|A^{-1}\|.
\]
Perturbation theory gives, schematically,
\[
\frac{\|\delta x\|}{\|x\|}
\lesssim
\kappa(A)
\left(
\frac{\|\delta A\|}{\|A\|}
+
\frac{\|\delta b\|}{\|b\|}
\right).
\]
Poor conditioning amplifies numerical perturbations.

## 4. Direct methods

Factorisation writes \(A=LU\) or a related factorisation. It is robust for moderate systems but memory-intensive for large CFD meshes.

## 5. Krylov methods

CG applies to suitable symmetric positive-definite systems. GMRES handles general nonsymmetric systems:
\[
x_m\in x_0+\mathcal K_m(A,r_0),
\]
where
\[
\mathcal K_m=span\{r_0,Ar_0,\ldots,A^{m-1}r_0\}.
\]

## 6. Preconditioning

Solve instead
\[
M^{-1}Ax=M^{-1}b,
\]
where \(M^{-1}\) approximates \(A^{-1}\) cheaply. The objective is not to change the physical equation but to reduce the difficulty seen by the Krylov method.

## 7. Residuals and energy norms

For SPD \(A\), an energy norm is
\[
\|e\|_A=\sqrt{e^TAe}.
\]
A residual contraction
\[
\frac{\|r_{k+1}\|}{\|r_k\|}<1
\]
does not necessarily imply an equivalent error contraction. AMG qualification should therefore use the norm specified by the mathematical contract.

## 8. Breakdown

Failures include singular systems, incompatible RHS, loss of definiteness, zero pivots, stagnation and excessive condition numbers. The solver must report the failure reason rather than silently switching to a different method.

**CFDX paths:** `src/cfdx/core/linalg/`, tests under `tests/unit/`.

## 9. Verification

Use manufactured linear systems, symmetry checks, true residual recomputation, known eigenvalue/conditioning cases and solver iteration histories.
