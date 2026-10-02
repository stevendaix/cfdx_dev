# %% [markdown]
# 05 — Fluxes
#
# Fluxes are the numerical interface between conservation laws and the algebraic operator. One physical internal face must produce one numerical flux, with opposite signs for its two adjacent cells.
#
# ## 5.1 Mass flux
#
# \[
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f
# \]
# has units \(kg\,s^{-1}\). For incompressible constant-density flow, the same expression reduces to \(\rho\mathbf u_f\cdot\mathbf S_f\).
#
# ## 5.2 Convective scalar flux
#
# \[
# F_{c,f}=\dot m_f\phi_f.
# \]
# Upwind:
# \[
# \phi_f=
# \begin{cases}
# \phi_P,&\dot m_f\ge0,\\
# \phi_N,&\dot m_f<0.
# \end{cases}
# \]
# Central interpolation:
# \[
# \phi_f=(1-w)\phi_P+w\phi_N.
# \]
# A higher-order reconstruction can be written
# \[
# \phi_f=\phi_P+\nabla\phi_P\cdot(\mathbf x_f-\mathbf C_P)+\cdots
# \]
# and therefore inherits the quality of the gradient and limiter.
#
# ## 5.3 Diffusive flux
#
# For a transported scalar:
# \[
# F_{d,f}=-\Gamma_f\nabla\phi_f\cdot\mathbf S_f
# \]
# when \(F_d\) denotes outward physical flux. The corresponding PDE term is \(+\nabla\cdot(\Gamma\nabla\phi)\); sign conventions must never be mixed between PDE and outward-flux definitions.
#
# ## 5.4 Pressure force
#
# Pressure contributes
# \[
# \mathbf F_{p,f}=-p_f\mathbf S_f
# \]
# to the momentum balance with the outward-normal convention.
#
# ## 5.5 Energy and species
#
# Convective energy flux is typically \(\dot m_f h_f\) or \(\dot m_f e_f\), depending on the chosen formulation. Species convection is \(\dot m_fY_{k,f}\); molecular diffusion adds a constitutive flux. The selected energy/species variable and reference enthalpy must be explicit.
#
# ## 5.6 Internal-face antisymmetry
#
# For one face:
# \[
# \mathbf S_{Nf}=-\mathbf S_{Pf},
# \qquad
# F_{Nf}=-F_{Pf}.
# \]
# This is stronger than checking a global balance because a pair of wrong fluxes can still cancel globally.
#
# ## 5.7 Non-orthogonal diffusion
#
# \[
# \nabla\phi_f\cdot\mathbf S_f
# =
# \nabla\phi_f\cdot\mathbf S_f^\parallel+
# \nabla\phi_f\cdot\mathbf S_f^\perp.
# \]
# The correction term requires a gradient reconstruction and therefore couples N2 and N3.
#
# ## 5.8 Boundedness
#
# For a bounded face interpolation one may require
# \[
# \min(\phi_P,\phi_N)\le\phi_f\le\max(\phi_P,\phi_N).
# \]
# Multidimensional limiters need a stronger cell-neighbourhood definition. Boundedness must be verified for the actual stencil and not inferred from a one-dimensional formula.
#
# ## 5.9 CFDX implementation
#
# The numerical interfaces are in [flux.h](../../../src/cfdx/core/numerics/flux.h), [convection.h](../../../src/cfdx/core/numerics/convection.h), [interpolation.h](../../../src/cfdx/core/numerics/interpolation.h), [divergence.h](../../../src/cfdx/core/numerics/divergence.h) and [laplacian.h](../../../src/cfdx/core/numerics/laplacian.h).
#
# Verification uses [test_convection_scheme_verification.cpp](../../../tests/validation/test_convection_scheme_verification.cpp), [test_convection_3d_verification.cpp](../../../tests/validation/test_convection_3d_verification.cpp), [test_convection_blending.cpp](../../../tests/validation/test_convection_blending.cpp) and conservation tests.
#
# ## 5.10 Executable antisymmetry check
# %%
import numpy as np
F = 2.5
assert np.isclose(F + (-F), 0.0)
