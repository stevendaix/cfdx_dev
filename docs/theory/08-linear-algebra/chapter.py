# %% [markdown]
# 08 — Linear Algebra
#
# ## 8.1 Discrete system
#
# CFD assembly produces
# \[
# Ax=b,\qquad r=b-Ax.
# \]
# The residual is an algebraic defect. It is not by itself a discretisation-error estimate.
#
# ## 8.2 Sparse structure
#
# A finite-volume row contains the owner coefficient and contributions from the stencil:
# \[
# a_Px_P+\sum_Na_{PN}x_N=b_P.
# \]
# Sparse storage records only non-zero coefficients. Matrix sparsity is dictated by mesh connectivity and block coupling.
#
# ## 8.3 Direct methods
#
# Gaussian elimination conceptually performs
# \[
# A=LU,
# \]
# followed by \(Ly=b\) and \(Ux=y\). Sparse direct methods reorder unknowns to control fill-in. Memory growth can make direct factorisation unsuitable for large 3-D CFD systems.
#
# ## 8.4 Jacobi and Gauss–Seidel
#
# Splitting \(A=D+L+U\):
# \[
# x^{k+1}=D^{-1}(b-(L+U)x^k)
# \]
# for Jacobi, while Gauss–Seidel uses updated lower-block values:
# \[
# (D+L)x^{k+1}=b-Ux^k.
# \]
# Their convergence depends on the spectrum of the iteration matrix.
#
# ## 8.5 CG
#
# For symmetric positive-definite \(A\), conjugate gradients minimises the error in
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
# This bound is a theoretical estimate, not a promise for an arbitrary CFD matrix.
#
# ## 8.6 GMRES
#
# \[
# \mathcal K_k(A,r_0)
# =\mathrm{span}\{r_0,Ar_0,\ldots,A^{k-1}r_0\}.
# \]
# GMRES chooses
# \[
# x_k=\arg\min_{x\in x_0+\mathcal K_k}\|b-Ax\|_2.
# \]
# Restarted GMRES replaces the full Krylov basis after a chosen dimension.
#
# ## 8.7 BiCGStab
#
# BiCGStab combines bi-conjugate-gradient ideas with stabilisation. It is useful for nonsymmetric systems but has different breakdown and convergence behaviour from GMRES; it should not be described as a cheaper equivalent.
#
# ## 8.8 Preconditioning
#
# Given \(M\approx A\):
# \[
# M^{-1}Ax=M^{-1}b.
# \]
# Left and right preconditioning produce different residual interpretations:
# \[
# M^{-1}Ax=M^{-1}b
# \quad\text{or}\quad
# AM^{-1}y=b,\;x=M^{-1}y.
# \]
# Diagnostics must state which residual is monitored.
#
# ## 8.9 True residual and scaling
#
# A robust verification periodically recomputes
# \[
# r_{true}=b-Ax
# \]
# from the original matrix. Relative residual:
# \[
# \eta=\frac{\|b-Ax\|}{\|b\|+\|A\|\|x\|}.
# \]
# Scaling is important because raw norms depend on units and variable magnitude.
#
# ## 8.10 Breakdown
#
# Failure can result from singularity, indefiniteness, loss of orthogonality, poor conditioning, NaN/Inf coefficients, incompatible BCs or an ineffective preconditioner. The diagnostic must identify the algebraic failure rather than simply increasing the iteration limit.
#
# ## 8.11 CFDX implementation
#
# See [sparse_matrix.h](../../../src/cfdx/core/linalg/sparse_matrix.h), [linear_system.h](../../../src/cfdx/core/linalg/linear_system.h), [cg_solver.h](../../../src/cfdx/core/linalg/cg_solver.h), [gmres_solver.h](../../../src/cfdx/core/linalg/gmres_solver.h), [bicgstab_solver.h](../../../src/cfdx/core/linalg/bicgstab_solver.h), and [linear_solver_dispatch.h](../../../src/cfdx/core/linalg/linear_solver_dispatch.h).
#
# Tests include the CG, GMRES and Krylov-preconditioning unit suites.
#
# ## 8.12 Executable exact solve
# %%
import numpy as np
A=np.array([[4.,1.],[1.,3.]])
b=np.array([1.,2.])
x=np.linalg.solve(A,b)
r=b-A@x
assert np.linalg.norm(r)<1e-12
