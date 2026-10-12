# 09 — AMG, MGR and Schur Complements

## 1. Schur reduction
For
\[
\mathcal A=
\begin{bmatrix}A&G\\D&C\end{bmatrix},
\]
eliminating \(u\) gives
\[
S=C-DA^{-1}G.
\]
For incompressible flow \(C=0\), so \(S=-DA^{-1}G\) up to the sign convention of the block equations.

**CFDX paths:** `src/cfdx/core/linalg/exact_schur.h`, `block_schur.h`, `schur_approximation.h`.

## 2. Approximate Schur
\[
\tilde S=C-D\tilde A^{-1}G.
\]
The approximation should preserve the dominant pressure/constraint behaviour. Qualification requires comparison with an exact or independently assembled oracle.

## 3. AMG hierarchy
Given fine matrix \(A_f\), transfer operators \(P\) and \(R\), the Galerkin coarse matrix is
\[
A_c=RA_fP.
\]
The coarse correction is
\[
x_f\leftarrow x_f+Pe_c,
\qquad
A_ce_c=R(b_f-A_fx_f).
\]

## 4. Smoothing
A smoother damps high-frequency/algebraically local error. AMG effectiveness depends on the complementarity between smoothing and coarse correction; neither component is sufficient alone.

## 5. V-cycle
\[
S_{pre}\rightarrow R\rightarrow A_c^{-1}\rightarrow P\rightarrow S_{post}.
\]
A hierarchy is useful only if it reduces the relevant error modes at acceptable setup and application cost.

## 6. Energy analysis
For SPD \(A\),
\[
\|e\|_A^2=e^TAe.
\]
A qualification measure can be
\[
\rho_E=\sqrt{\frac{e_{out}^TAe_{out}}{e_{in}^TAe_{in}}}.
\]
The matrix family, initial error and number of cycles must be recorded.

## 7. MGR
Multigrid reduction can eliminate subsets of variables in stages:
\[
N_{all}\rightarrow N_{velocity}\rightarrow N_{Schur}.
\]
Transfer operators and block ordering determine the resulting coarse operator.

## 8. Update semantics
If graph topology is unchanged but coefficients vary, a preconditioner may support value updates. The contract must distinguish rebuild, refresh and reuse; stale values are a correctness issue, not merely a performance issue.

**CFDX paths:** `src/cfdx/core/linalg/hypre_amg.*`, `amg_preconditioner.h`, `coupled_amg_schur.h`.

## 9. Verification
Use exact Schur matrices, Galerkin identities, energy contraction, block saddle-point systems and explicit update-values tests. Evidence is indexed by `docs/validation/AMG_PRECONDITIONER_QUALIFICATION.md`.