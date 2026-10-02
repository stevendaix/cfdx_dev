# %% [markdown]
# 01 — Conservation Laws
#
# ## 1.1 Reynolds transport theorem
#
# For an extensive property \(B=\int_V\rho b\,dV\):
# \[
# \frac{dB}{dt}
# =
# \frac{d}{dt}\int_V\rho b\,dV+
# \oint_{\partial V}\rho b\,\mathbf u\cdot\mathbf n\,dA.
# \]
# This theorem connects a material-system balance to a fixed control volume.
#
# ## 1.2 Generic local balance
#
# \[
# \frac{\partial q}{\partial t}
# +\nabla\cdot\mathbf F=s.
# \]
# Integrating:
# \[
# \frac{d}{dt}\int_Vq\,dV+
# \oint_{\partial V}\mathbf F\cdot\mathbf n\,dA
# =\int_Vs\,dV.
# \]
# CFD finite volume starts from this integral statement because it retains conservation at the control-volume level.
#
# ## 1.3 Mass
#
# \[
# \frac{\partial\rho}{\partial t}
# +\nabla\cdot(\rho\mathbf u)=0.
# \]
# The discrete cell balance is
# \[
# \frac{d}{dt}(\rho_PV_P)+\sum_f\dot m_f=0,
# \qquad
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
# \]
# For steady incompressible flow:
# \[
# \sum_f\dot m_f=0.
# \]
#
# ## 1.4 Momentum
#
# \[
# \frac{\partial(\rho\mathbf u)}{\partial t}
# +\nabla\cdot(\rho\mathbf u\otimes\mathbf u)
# =
# -\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
# \]
# Integral momentum contains convective transport, pressure traction
# \[
# -\oint p\mathbf n\,dA,
# \]
# viscous traction
# \[
# \oint\boldsymbol\tau\mathbf n\,dA,
# \]
# and body force.
#
# ## 1.5 Angular momentum
#
# Conservation of angular momentum for a classical continuum implies symmetry of the Cauchy stress in the absence of couple stresses:
# \[
# \boldsymbol\sigma=\boldsymbol\sigma^T.
# \]
# This is a constitutive/continuum assumption that underlies standard Navier–Stokes CFD.
#
# ## 1.6 Energy
#
# \[
# E=e+\frac12|\mathbf u|^2,
# \]
# \[
# \frac{\partial(\rho E)}{\partial t}
# +\nabla\cdot[(\rho E+p)\mathbf u]
# =
# \nabla\cdot(\boldsymbol\tau\mathbf u-\mathbf q)
# +\rho\mathbf f\cdot\mathbf u+\dot q_v.
# \]
# The terms correspond to storage, enthalpy transport, pressure work, viscous work, heat conduction and volumetric sources.
#
# ## 1.7 Species
#
# \[
# \frac{\partial(\rho Y_k)}{\partial t}
# +\nabla\cdot(\rho\mathbf uY_k)
# =-\nabla\cdot\mathbf J_k+\dot\omega_k.
# \]
# Mixture constraints:
# \[
# \sum_kY_k=1,\qquad
# \sum_k\mathbf J_k=0,\qquad
# \sum_k\dot\omega_k=0
# \]
# for a closed conservative mixture formulation.
#
# ## 1.8 Flux decomposition
#
# A generic transported quantity has
# \[
# \mathbf F=
# \underbrace{\rho\mathbf u\phi}_{convection}
# +\underbrace{\mathbf F_d}_{diffusion}
# +\underbrace{\mathbf F_p}_{pressure/constraint}
# +\underbrace{\mathbf F_{other}}_{physics}.
# \]
# The exact decomposition depends on the equation.
#
# ## 1.9 Boundary integrals
#
# A physical-domain balance contains only external boundary fluxes after internal faces cancel. For steady state:
# \[
# \sum_{\partial\Omega}F_f=\int_\Omega s\,dV.
# \]
# This provides a global conservation oracle when all boundary contributions are independently known.
#
# ## 1.10 Internal-face cancellation
#
# For an internal face:
# \[
# \mathbf S_{f,N}=-\mathbf S_{f,P},
# \qquad
# F_{N,f}=-F_{P,f}.
# \]
# Summing over all cells cancels internal contributions exactly. This is the discrete counterpart of locality plus conservation.
#
# ## 1.11 Conservative versus non-conservative forms
#
# The identity
# \[
# \nabla\cdot(\rho\mathbf u\phi)
# =
# \rho\mathbf u\cdot\nabla\phi+
# \phi\nabla\cdot(\rho\mathbf u)
# \]
# shows that conservative and advective forms are equivalent only when the continuity equation is satisfied exactly. In discrete multiphysics, conservative assembly is usually preferred because it exposes the shared face flux.
#
# ## 1.12 Conservation versus boundedness
#
# Conservation:
# \[
# \sum_{cells}R_P=\text{boundary/source balance}.
# \]
# Boundedness concerns admissible values:
# \[
# \phi_{min}\le\phi_h\le\phi_{max}.
# \]
# They are independent numerical properties.
#
# ## 1.13 Verification hierarchy
#
# \[
# geometry\ closure
# \rightarrow face\ antisymmetry
# \rightarrow local\ balance
# \rightarrow global\ balance
# \rightarrow MMS
# \rightarrow refinement
# \rightarrow coupled\ conservation.
# \]
# A small solver residual cannot replace these checks.
#
# ## 1.14 CFDX traceability
#
# Representative operators are [flux.h](../../../src/cfdx/core/numerics/flux.h), [integrate.h](../../../src/cfdx/core/numerics/integrate.h), [divergence.h](../../../src/cfdx/core/numerics/divergence.h), [finite_volume_transport.h](../../../src/cfdx/physics/finite_volume_transport.h) and the physics-specific conservation models.
#
# ## 1.15 Executable antisymmetry
# %%
import numpy as np
F=2.75
assert np.isclose(F+(-F),0.0)
