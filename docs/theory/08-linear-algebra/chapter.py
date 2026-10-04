# %% [markdown]
# 08 — Linear Algebra
#
# ## 8.1 Discrete system
#
# CFD assembly produces
# \[
# Ax=b,\qquad r=b-Ax.
# \]
# The residual is an algebraic defect, not a discretisation-error estimate.
#
# ## 8.2 Sparse structure and block systems
#
# A scalar finite-volume row is
# \[
# a_Px_P+\sum_Na_{PN}x_N=b_P.
# \]
# For coupled velocity/pressure/energy variables, the same graph becomes a block matrix. A block preconditioner can exploit this structure without changing the underlying assembled equations.
#
# ## 8.3 Direct and stationary iterations
#
# Gaussian elimination conceptually gives
# \[
# A=LU,\qquad Ly=b,\qquad Ux=y.
# \]
# With \(A=D+L+U\), Jacobi is
# \[
# x^{k+1}=D^{-1}[b-(L+U)x^k],
# \]
# while Gauss–Seidel is
# \[
# (D+L)x^{k+1}=b-Ux^k.
# \]
# Convergence is governed by the spectral radius
# \[
# \rho(T)<1
# \]
# of the corresponding iteration matrix \(T\).
#
# ## 8.4 Conjugate gradients
#
# For SPD \(A\), CG constructs \(A\)-orthogonal search directions and minimises
# \[
# \|e\|_A=\sqrt{e^TAe}.
# \]
# The classical bound is
# \[
# \|e_k\|_A\le
# 2\left(
# \frac{\sqrt{\kappa(A)}-1}{\sqrt{\kappa(A)}+1}
# \right)^k\|e_0\|_A.
# \]
# The SPD requirement matters: applying CG to a genuinely nonsymmetric or indefinite matrix invalidates the standard theory.
#
# ## 8.5 GMRES
#
# \[
# \mathcal K_k(A,r_0)
# =\operatorname{span}\{r_0,Ar_0,\ldots,A^{k-1}r_0\},
# \]
# \[
# x_k=\arg\min_{x\in x_0+\mathcal K_k}\|b-Ax\|_2.
# \]
# The Arnoldi relation can be written
# \[
# AV_k=V_{k+1}\bar H_k.
# \]
# Restarted GMRES limits memory by replacing the Krylov basis after \(m\) iterations, at the cost of potentially losing useful spectral information.
#
# ## 8.6 FGMRES
#
# Flexible GMRES permits a varying preconditioner:
# \[
# z_k=M_k^{-1}v_k,
# \]
# rather than requiring one fixed \(M^{-1}\). This is important when AMG, MGR or inner nonlinear solves change between iterations. It is not interchangeable with ordinary GMRES if the preconditioner is variable.
#
# ## 8.7 BiCGStab
#
# BiCGStab targets nonsymmetric systems using a stabilised bi-Lanczos-type recurrence. It can converge rapidly on some CFD matrices but can also experience breakdown or irregular residual histories. The implementation must distinguish true-residual failure from recurrence breakdown.
#
# ## 8.8 Preconditioning
#
# Given \(M\approx A\),
# \[
# M^{-1}Ax=M^{-1}b.
# \]
# Right preconditioning gives
# \[
# AM^{-1}y=b,\qquad x=M^{-1}y.
# \]
# The solver's stopping metric must therefore identify whether it monitors a preconditioned or true residual.
#
# ## 8.9 True residual, backward error and scaling
#
# The true residual is
# \[
# r=b-Ax.
# \]
# A scale-aware backward-error measure is
# \[
# \eta=
# \frac{\|b-Ax\|}
# {\|b\|+\|A\|\|x\|}.
# \]
# This is often more meaningful than a raw absolute residual for variables with different units.
#
# ## 8.10 Null spaces and singular systems
#
# Pressure Poisson operators with pure Neumann conditions have a constant null mode:
# \[
# A\mathbf1=0.
# \]
# A Krylov solver must either operate in the compatible subspace or receive an explicit null-space/gauge treatment. A generic convergence criterion cannot repair an incompatible right-hand side.
#
# ## 8.11 Conditioning and non-normality
#
# \[
# \kappa(A)=\|A\|\|A^{-1}\|.
# \]
# For nonsymmetric/non-normal matrices, eigenvalues alone do not fully determine transient Krylov behaviour. Pseudospectral sensitivity and field-of-values arguments can be more relevant than the SPD CG estimate.
#
# ## 8.12 Breakdown diagnostics
#
# Failure can result from singularity, indefiniteness, loss of orthogonality, NaN/Inf coefficients, incompatible BCs, stagnation or an ineffective preconditioner. Increasing the iteration limit without identifying the mechanism is not a numerical fix.
#
# ## 8.13 CFDX implementation
#
# See [sparse_matrix.h](../../../src/cfdx/core/linalg/sparse_matrix.h), [linear_system.h](../../../src/cfdx/core/linalg/linear_system.h), [cg_solver.h](../../../src/cfdx/core/linalg/cg_solver.h), [gmres_solver.h](../../../src/cfdx/core/linalg/gmres_solver.h), [bicgstab_solver.h](../../../src/cfdx/core/linalg/bicgstab_solver.h), [linear_solver_models.h](../../../src/cfdx/core/linalg/linear_solver_models.h) and [linear_solver_dispatch.h](../../../src/cfdx/core/linalg/linear_solver_dispatch.h).
#
# Current method-registry evidence distinguishes implemented/verified methods; a solver enum alone is not qualification evidence.
#
# ## 8.14 Executable exact solve
# %%
import numpy as np
A=np.array([[4.,1.],[1.,3.]])
b=np.array([1.,2.])
x=np.linalg.solve(A,b)
r=b-A@x
assert np.linalg.norm(r)<1e-12
