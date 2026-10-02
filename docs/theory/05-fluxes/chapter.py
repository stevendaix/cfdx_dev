# %% [markdown]
# 05 — Fluxes
#
# ## 5.1 Conservation meaning
#
# For a control volume,
# [
# rac{d}{dt}int_Vhophi,dV
# =-oint_{partial V}mathbf F_phicdotmathbf n,dA
# +int_VS_phi,dV.
# ]
# A numerical face flux is therefore a discrete approximation to the surface integral. Its orientation must be explicit.
#
# ## 5.2 Mass and momentum fluxes
#
# [
# dot m_f=ho_fmathbf u_fcdotmathbf S_f.
# ]
# For momentum, the convective flux is
# [
# mathbf F_{m,f}=dot m_fmathbf u_f,
# ]
# while pressure contributes
# [
# mathbf F_{p,f}=-p_fmathbf S_f.
# ]
# Viscous traction for a Newtonian fluid is
# [
# oldsymbol	au=2mumathbf S
# -rac23mu(
ablacdotmathbf u)I,
# qquad
# mathbf F_{	au,f}=oldsymbol	au_fmathbf S_f.
# ]
#
# ## 5.3 Scalar convection
#
# [
# F_{c,f}=dot m_fphi_f.
# ]
# Upwind:
# [
# phi_f=
# egin{cases}
# phi_P,&dot m_fge0,\
# phi_N,&dot m_f<0.
# end{cases}
# ]
# Central:
# [
# phi_f=(1-w)phi_P+wphi_N.
# ]
# A reconstructed scheme has the generic form
# [
# phi_f=phi_P+
# 
ablaphi_Pcdot(mathbf x_f-mathbf C_P)+
# O(h^2).
# ]
#
# ## 5.4 Diffusive flux and non-orthogonal decomposition
#
# [
# F_{d,f}=-Gamma_f
ablaphi_fcdotmathbf S_f.
# ]
# With
# [
# mathbf S_f=mathbf S_f^parallel+mathbf S_f^perp,
# ]
# [
# F_{d,f}
# =-Gamma_f
ablaphi_fcdotmathbf S_f^parallel
# -Gamma_f
ablaphi_fcdotmathbf S_f^perp.
# ]
# The first term can use a two-point centre difference:
# [
# 
ablaphi_fcdotmathbf S_f^parallel
# approx
# rac{phi_N-phi_P}{|mathbf d_{PN}|}
# rac{mathbf d_{PN}cdotmathbf S_f}{|mathbf d_{PN}|}.
# ]
# The correction term depends on a reconstructed gradient and is therefore sensitive to N2.
#
# ## 5.5 Flux antisymmetry
#
# For an internal face:
# [
# mathbf S_{Nf}=-mathbf S_{Pf},
# qquad
# F_{Nf}=-F_{Pf}.
# ]
# This must be enforced at the shared-face level, not recreated independently from two cell rows.
#
# ## 5.6 Boundary fluxes
#
# A prescribed outward flux (q_Gamma) gives
# [
# F_Gamma=q_Gamma A_f.
# ]
# A Dirichlet condition instead determines a boundary state/value and therefore enters through the selected reconstruction and constitutive flux. Robin conditions combine value and flux, for example
# [
# -k
abla Tcdot n=h(T-T_infty).
# ]
#
# ## 5.7 Boundedness and TVD concepts
#
# Pairwise boundedness requires
# [
# min(phi_P,phi_N)lephi_flemax(phi_P,phi_N).
# ]
# For one-dimensional linear advection, a TVD analysis is often expressed using the Courant number
# [
# C=rac{uDelta t}{Delta x}
# ]
# and a flux-limiter function (psi(r)), where
# [
# r=rac{phi_i-phi_{i-1}}{phi_{i+1}-phi_i}.
# ]
# The limiter must satisfy the appropriate Sweby-region constraints for the declared time/space discretisation; a limiter kernel alone is not a complete TVD proof.
#
# ## 5.8 Flux linearisation
#
# After face interpolation, the flux contributes to the algebraic row:
# [
# a_Pphi_P+sum_Na_{PN}phi_N=b_P.
# ]
# The coefficient signs depend on the flux convention and on whether the term is treated explicitly, implicitly, or deferred-correction. The code-level coefficient contract is therefore part of the numerical-method definition.
#
# ## 5.9 CFDX implementation
#
# Core interfaces: [flux.h](../../../src/cfdx/core/numerics/flux.h), [convection.h](../../../src/cfdx/core/numerics/convection.h), [interpolation.h](../../../src/cfdx/core/numerics/interpolation.h), [divergence.h](../../../src/cfdx/core/numerics/divergence.h), [laplacian.h](../../../src/cfdx/core/numerics/laplacian.h).
#
# Verification: [test_convection_scheme_verification.cpp](../../../tests/validation/test_convection_scheme_verification.cpp), [test_convection_3d_verification.cpp](../../../tests/validation/test_convection_3d_verification.cpp), [test_convection_blending.cpp](../../../tests/validation/test_convection_blending.cpp), plus conservation/boundedness campaigns.
#
# ## 5.10 Executable antisymmetry
# %%
import numpy as np
F = 2.5
assert np.isclose(F + (-F), 0.0)
