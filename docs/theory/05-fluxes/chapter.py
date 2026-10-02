# %% [markdown]
# Fluxes
#
# ## 5.1 Conservative interface quantities
#
# \[
# \mathbf S_f=A_f\mathbf n_f,\qquad
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
# \]
# Convective scalar flux:
# \[
# F_{c,f}=\dot m_f\phi_f.
# \]
# Diffusive flux:
# \[
# F_{d,f}=-\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
# \]
# Pressure force:
# \[
# \mathbf F_{p,f}=-p_f\mathbf S_f.
# \]
#
# ## 5.2 Convection
#
# Upwind uses the upstream state:
# \[
# \phi_f=
# \begin{cases}\phi_P,&\dot m_f>0,\\
# \phi_N,&\dot m_f<0.
# \end{cases}
# \]
# Linear interpolation is
# \[
# \phi_f=(1-w)\phi_P+w\phi_N.
# \]
# Higher-order reconstruction uses gradients and therefore inherits gradient/stencil errors.
#
# ## 5.3 Antisymmetry
#
# For one internal physical face,
# \[
# \mathbf S_{f,N}=-\mathbf S_{f,P},
# \qquad
# F_{N,f}=-F_{P,f}.
# \]
# This identity is a local conservation invariant and should be tested directly.
#
# ## 5.4 Diffusion
#
# With centre vector \(\mathbf d_{PN}\), orthogonal two-point diffusion uses the normal derivative along the centre line. On non-orthogonal meshes,
# \[
# \nabla\phi\cdot\mathbf S_f
# =
# \nabla\phi\cdot\mathbf S_f^\parallel+
# \nabla\phi\cdot\mathbf S_f^\perp,
# \]
# so a correction is required.
#
# ## 5.5 Verification
#
# Test constants, linear fields, face antisymmetry, conservation, boundedness where claimed, and mesh refinement. A global balance cannot replace the local antisymmetry test because compensating errors can cancel globally.
#
# ## 5.6 CFDX traceability
#
# Flux: src/cfdx/core/numerics/flux.h  
# Convection: src/cfdx/core/numerics/convection.h  
# Interpolation: src/cfdx/core/numerics/interpolation.h  
# FVM transport: src/cfdx/physics/finite_volume_transport.h  
# Tests: tests/unit/test_flux.cpp, tests/unit/test_conservation_boundedness.cpp, tests/validation/test_convection_scheme_verification.cpp
#
# %%
from __future__ import annotations
import numpy as np
F=2.5
assert np.isclose(F+(-F),0.0)
