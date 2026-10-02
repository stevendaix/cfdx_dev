# 08 — Linear Algebra

## 1. Discrete system
Every implicit equation produces
\[
Ax=b,\qquad r=b-Ax.
\]
The true residual must be recomputable independently from the assembled operator.

## 2. Sparse matrix
Finite volume gives local coupling, hence most entries of \(A\) are zero. CSR-style storage uses values, column indices and row offsets. Matrix topology and matrix values must be treated as different update classes.

## 3. Conditioning
\[
\kappa(A)=\|A\|\|A^{-1}\|.
\]
For perturbations,
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
A large condition number amplifies coefficient and round-off errors.

## 4. Krylov methods
\[
\mathcal K_m(A,r_0)=span\{r_0,Ar_0,\ldots,A^{m-1}r_0\}.
\]
CG is appropriate for SPD matrices; GMRES targets general systems. Restarted GMRES limits memory but changes the Krylov process.

## 5. Preconditioning
\[
M^{-1}Ax=M^{-1}b.
\]
A good preconditioner clusters/restructures the spectrum without changing the physical solution. Left, right and split preconditioning have different residual interpretations.

## 6. Energy norm
For SPD \(A\),
\[
\|e\|_A=\sqrt{e^TAe}.
\]
Energy contraction is
\[
\rho_E=\frac{\|e_{out}\|_A}{\|e_{in}\|_A}.
\]
This is distinct from residual contraction \(\|r_{out}\|/\|r_{in}\|\).

## 7. Failure modes
A solver can fail through singularity, incompatible RHS, loss of definiteness, stagnation, breakdown or excessive iteration count. Automatic solver substitution must not hide the original failure.

**CFDX path:** `src/cfdx/core/linalg/`.

## 8. Verification
Use known matrices, exact solutions, independent residuals, symmetry/SPD checks, condition estimates and iteration histories. AMG and Schur are covered in chapter 09.