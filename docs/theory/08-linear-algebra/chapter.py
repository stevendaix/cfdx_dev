# %% [markdown]
# Linear Algebra
#
# ## 8.1 Algebraic system
# \[
# Ax=b,\qquad r=b-Ax.
# \]
# A residual is an algebraic diagnostic; it is not a direct estimate of discretisation or modelling error.
#
# ## 8.2 CG
#
# For SPD A, CG uses conjugate directions and the energy norm
# \[
# \|e\|_A=\sqrt{e^TAe}.
# \]
#
# ## 8.3 GMRES
# \[
# \mathcal K_k(A,r_0)=
# \operatorname{span}\{r_0,Ar_0,\ldots,A^{k-1}r_0\}.
# \]
# GMRES minimises the residual over the current Krylov space, subject to restart/preconditioning.
#
# ## 8.4 Preconditioning
# \[
# M^{-1}Ax=M^{-1}b.
# \]
# Preconditioning changes the algebraic representation used for convergence while preserving the target solution when applied consistently.
#
# ## 8.5 Conditioning
# \[
# \kappa(A)=\|A\|\|A^{-1}\|.
# \]
# Conditioning controls sensitivity to perturbations. For non-normal systems, it does not by itself predict Krylov convergence.
#
# ## 8.6 True residual
#
# Recompute
# \[
# r_{\mathrm{true}}=b-Ax
# \]
# from the original operator. A preconditioned residual is not automatically equivalent.
#
# ## 8.7 CFDX traceability
# src/cfdx/core/linalg/sparse_matrix.h  
# src/cfdx/core/linalg/linear_system.h  
# src/cfdx/core/linalg/cg_solver.h  
# src/cfdx/core/linalg/gmres_solver.h  
# src/cfdx/core/linalg/linear_solver_dispatch.h  
# Tests: tests/unit/test_cg_solver.cpp, tests/unit/test_gmres_solver.cpp, tests/unit/test_krylov_preconditioning.cpp
#
# %%
from __future__ import annotations
import numpy as np
A=np.array([[2.,1.],[1.,3.]])
b=np.array([1.,2.])
x=np.linalg.solve(A,b)
assert np.allclose(b-A@x,0.0)
