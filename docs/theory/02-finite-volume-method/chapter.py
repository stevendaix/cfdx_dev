# %% [markdown]
# Finite-Volume Method
#
# ## 2.1 Integral formulation
#
# For
# \[
# \frac{\partial(\rho\phi)}{\partial t}+\nabla\cdot(\rho\mathbf u\phi)
# =\nabla\cdot(\Gamma\nabla\phi)+S,
# \]
# integration over a control volume gives
# \[
# \frac{d}{dt}\int_{V_P}\rho\phi\,dV+
# \sum_f\int_{A_f}\rho\mathbf u\phi\cdot\mathbf n_f\,dA
# =
# \sum_f\int_{A_f}\Gamma\nabla\phi\cdot\mathbf n_f\,dA+
# \int_{V_P}S\,dV.
# \]
# This is the central FVM principle: conservation is discretised as a balance over a finite volume, not reconstructed from pointwise derivatives.
#
# ## 2.2 Geometry
#
# \[
# \mathbf S_f=A_f\mathbf n_f,\qquad
# \sum_f\mathbf S_f=\mathbf0.
# \]
# The second identity is the closure condition of a closed polyhedron.
#
# ## 2.3 Face fluxes
#
# \[
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f,
# \qquad
# F_{c,f}=\dot m_f\phi_f,
# \]
# \[
# F_{d,f}=-\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
# \]
# One physical internal face must have one conservative flux with opposite orientation in its two neighbouring cells.
#
# ## 2.4 Algebraic equation
#
# After spatial approximation and source linearisation,
# \[
# a_P\phi_P=\sum_Na_{PN}\phi_N+b_P.
# \]
# The coefficients depend on convection, diffusion, temporal discretisation, source treatment and boundary conditions. Properties such as diagonal dominance or an M-matrix structure must therefore be demonstrated for the specific scheme.
#
# ## 2.5 Diffusion and non-orthogonality
#
# With \(\mathbf d_{PN}=\mathbf C_N-\mathbf C_P\),
# \[
# \nabla\phi\cdot\mathbf S_f
# =
# \nabla\phi\cdot\mathbf S_f^\parallel+
# \nabla\phi\cdot\mathbf S_f^\perp.
# \]
# The perpendicular correction is a geometric consequence of non-orthogonality. Its accuracy depends directly on the gradient reconstruction.
#
# ## 2.6 Boundedness
#
# A scheme claiming a local maximum principle must define its admissible bounds and demonstrate them. A generic local requirement is
# \[
# \min_N\phi_N\le\phi_f\le\max_N\phi_N.
# \]
# High-order interpolation can violate this and therefore needs limiting when boundedness is required.
#
# ## 2.7 Verification hierarchy
#
# First verify geometric closure and constant-field invariance. Then verify linear-field exactness, conservation, MMS, refinement order and finally coupled-physics behaviour. Passing a solver benchmark does not replace operator-level verification.
#
# ## 2.8 CFDX traceability
#
# Discrete diffusion: src/cfdx/core/numerics/laplacian.h  
# Divergence: src/cfdx/core/numerics/divergence.h  
# Convection: src/cfdx/core/numerics/convection.h  
# Sources: src/cfdx/core/numerics/source_term.h  
# FVM transport: src/cfdx/physics/finite_volume_transport.h
#
# The relative source links are intentionally kept beside the equations so that the theory remains directly auditable.
#
# %%
from __future__ import annotations
import numpy as np
S=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
assert np.allclose(S.sum(axis=0),0.0)
