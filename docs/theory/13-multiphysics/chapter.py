# %% [markdown]
# Multiphysics
#
# ## 13.1 Coupled residual
# \[
# F(y)=0,\qquad y=(u,p,T,Y_1,\ldots).
# \]
# Newton linearisation:
# \[
# J(y^k)\delta y=-F(y^k),\qquad y^{k+1}=y^k+\delta y.
# \]
#
# ## 13.2 Block Jacobian
# \[
# J=
# \begin{bmatrix}
# J_{uu}&J_{up}&J_{uT}&\cdots\\
# J_{pu}&J_{pp}&J_{pT}&\cdots\\
# J_{Tu}&J_{Tp}&J_{TT}&\cdots
# \end{bmatrix}.
# \]
# Off-diagonal blocks represent coupling.
#
# ## 13.3 Segregated coupling
# \[
# y^{k+1}=S(y^k),\qquad
# y^{k+1}\leftarrow(1-\omega)y^k+\omega S(y^k).
# \]
# Relaxation changes convergence path, not the underlying equations.
#
# ## 13.4 Interface conservation
# \[
# R_\Gamma=\sum_{f\in\Gamma}(q_f^{(1)}+q_f^{(2)}).
# \]
# This defect must be checked independently from field residuals.
#
# ## 13.5 CFDX traceability
# src/cfdx/physics/finite_volume_transport.h  
# src/cfdx/physics/cht_solver.h  
# src/cfdx/physics/radiation_solver.h  
# src/cfdx/physics/steady_incompressible_solver.h  
# Tests: tests/validation/test_level_c_coupled_verification.cpp, tests/validation/test_cht_validation.cpp
#
# %%
from __future__ import annotations
import numpy as np
J=np.array([[3.,1.],[1.,2.]])
F=np.array([1.,-1.])
d=np.linalg.solve(J,-F)
assert np.allclose(J@d,-F)
