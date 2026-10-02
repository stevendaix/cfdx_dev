# %% [markdown]
# AMG, MGR and Schur Complement Methods
#
# ## 9.1 Exact Schur
# \[
# A=\begin{bmatrix}A_{uu}&A_{up}\\A_{pu}&A_{pp}\end{bmatrix},
# \qquad
# S=A_{pp}-A_{pu}A_{uu}^{-1}A_{up}.
# \]
#
# ## 9.2 Approximate Schur
# \[
# \widetilde S=A_{pp}-A_{pu}M_u^{-1}A_{up}.
# \]
# The quality of the approximation should be measured against exact/reference Schur action where possible.
#
# ## 9.3 Galerkin AMG
# \[
# A_H=R A_hP.
# \]
# A V-cycle combines smoothing, restriction, coarse correction and prolongation.
#
# ## 9.4 Energy contraction
# \[
# \|e\|_A=\sqrt{e^TAe},
# \qquad
# q_E=\frac{\|e_{out}\|_A}{\|e_{in}\|_A}.
# \]
# Energy contraction is a structural preconditioner diagnostic; residual decrease alone is insufficient.
#
# ## 9.5 MGR
#
# MGR selects coarse variables and reduces the remaining blocks. Variable ordering, transfer operators, relaxation and update semantics are part of the algorithm contract.
#
# ## 9.6 Matrix updates
#
# If coefficient values change while sparsity remains unchanged, numerical hierarchy data may need refreshing. If connectivity changes, the hierarchy structure itself may require rebuilding.
#
# ## 9.7 CFDX traceability
# src/cfdx/core/linalg/exact_schur.h  
# src/cfdx/core/linalg/block_schur.h  
# src/cfdx/core/linalg/schur_approximation.h  
# src/cfdx/core/linalg/amg_preconditioner.h  
# src/cfdx/core/linalg/mgr_preconditioner.h  
# Tests: tests/unit/test_exact_schur.cpp, tests/unit/test_schur_preconditioner.cpp, tests/unit/test_mgr_preconditioner.cpp, tests/validation/test_coupled_block_schur_amg.cpp
#
# %%
from __future__ import annotations
import numpy as np
A=np.array([[4.,1.],[1.,3.]])
e=np.array([1.,-2.])
assert np.sqrt(e@A@e)>0
