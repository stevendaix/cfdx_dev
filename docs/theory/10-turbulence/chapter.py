# %% [markdown]
# Turbulence
#
# ## 10.1 Reynolds decomposition
# \[
# u_i=\overline u_i+u_i',\qquad \overline{u_i'}=0.
# \]
# Averaging nonlinear convection creates
# \[
# -\rho\overline{u_i'u_j'},
# \]
# producing the closure problem.
#
# ## 10.2 Eddy viscosity
# \[
# -\rho\overline{u_i'u_j'}
# =2\mu_tS_{ij}-\frac23\rho k\delta_{ij},
# \qquad
# k=\frac12\overline{u_i'u_i'}.
# \]
# This is a model assumption.
#
# ## 10.3 RANS transport
# A generic turbulence scalar has
# \[
# \frac{Dq}{Dt}=P_q-D_q+\nabla\cdot(D_{eff}\nabla q).
# \]
# Every model-specific production, destruction and cross-diffusion term must be dimensionally and sign consistent.
#
# ## 10.4 SST
# SST combines k-omega near-wall behaviour with k-epsilon-like outer behaviour using blending functions and model-specific stress limitation/cross-diffusion. Qualification must test the exact equations and constants implemented.
#
# ## 10.5 Spalart--Allmaras
# SA transports a modified viscosity variable and includes wall-distance-dependent terms. Wall-distance quality is consequently part of the model numerical chain.
#
# ## 10.6 LES and hybrid models
# Filtering gives
# \[
# \bar u_i=G*u_i,\qquad
# \tau_{ij}^{sgs}=\overline{u_iu_j}-\bar u_i\bar u_j.
# \]
# DES/DDES/IDDES introduce grid/wall-dependent model length scales and therefore require grid-sensitivity evidence.
#
# ## 10.7 CFDX traceability
# src/cfdx/physics/turbulence.h  
# src/cfdx/physics/turbulence_models.h  
# src/cfdx/physics/turbulence_transport.h  
# src/cfdx/physics/sst_solver.h  
# src/cfdx/physics/spalart_allmaras.h  
# src/cfdx/physics/wall_distance.h  
# Tests: tests/unit/test_turbulence_qualification_equations.cpp, tests/unit/test_sst_blending.cpp
#
# %%
from __future__ import annotations
import numpy as np
uprime=np.array([1.,-1.,0.,2.,-2.])
assert np.isclose(np.mean(uprime),0.0)
