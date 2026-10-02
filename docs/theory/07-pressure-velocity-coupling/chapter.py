# %% [markdown]
# Pressure--Velocity Coupling
#
# ## 7.1 Constraint
# \[
# \nabla\cdot\mathbf u=0.
# \]
# Pressure is the Lagrange multiplier enforcing the incompressibility constraint.
#
# ## 7.2 Block system and Schur complement
# \[
# \begin{bmatrix}A_u&G\\D&C\end{bmatrix}
# \begin{bmatrix}u\\p\end{bmatrix}
# =
# \begin{bmatrix}b_u\\b_p\end{bmatrix},
# \]
# \[
# S=C-D A_u^{-1}G.
# \]
# The pressure equation after elimination is \(Sp=b_p-D A_u^{-1}b_u\).
#
# ## 7.3 SIMPLE and SIMPLEC
#
# SIMPLE derives a pressure correction using an approximate inverse of the momentum operator. SIMPLEC modifies the velocity-correction approximation to retain coupling differently. They must be treated as distinct numerical algorithms.
#
# ## 7.4 PISO and PIMPLE
#
# PISO applies multiple pressure corrections in one outer step. PIMPLE combines outer nonlinear iterations with PISO-like corrections. The exact loop ordering, relaxation and stopping criteria define the implemented algorithm.
#
# ## 7.5 Fractional step
# \[
# u^*=u^n+\Delta tR(u^n),
# \]
# \[
# \nabla^2p^{n+1}=\frac{\rho}{\Delta t}\nabla\cdot u^*,
# \]
# \[
# u^{n+1}=u^*-\frac{\Delta t}{\rho}\nabla p^{n+1}.
# \]
#
# ## 7.6 Pressure null space
# Under pure Neumann pressure conditions,
# \[
# p\mapsto p+C
# \]
# leaves velocity unchanged. Gauge fixing or explicit null-space handling is therefore mandatory.
#
# ## 7.7 Rhie--Chow
#
# Collocated meshes require a face mass-flux interpolation that prevents checkerboard pressure modes. The CFDX implementation must be documented from the actual discrete formula.
#
# ## 7.8 CFDX traceability
# src/cfdx/physics/pressure_velocity.h  
# src/cfdx/physics/pressure_velocity_algorithms.h  
# src/cfdx/physics/steady_incompressible_solver.h  
# Tests: tests/unit/test_incompressible.cpp, tests/validation/test_steady_incompressible_solver.cpp
#
# %%
from __future__ import annotations
import numpy as np
A=np.array([[3.,1.],[1.,2.]])
assert np.all(np.linalg.eigvalsh(A)>0)
