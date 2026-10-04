# %% [markdown]
# # 05 — Fluxes
#
# ## 5.1 Conservative meaning
#
# For a conserved scalar,
# \[
# \frac{d}{dt}\int_V \rho\phi\,dV
# =-\oint_{\partial V}\mathbf F_\phi\cdot\mathbf n\,dA
# +\int_V S_\phi\,dV.
# \]
#
# The finite-volume method therefore acts on **integrated face fluxes**. A shared internal face must contribute equal and opposite amounts to its owner and neighbour equations.
#
# ## 5.2 Mass, momentum and scalar fluxes
#
# With \(\mathbf S_f=A_f\mathbf n_f\) oriented owner-to-neighbour,
# \[
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
# \]
# Momentum convection is \(\dot m_f\mathbf u_f\). Pressure contributes \(p_f\mathbf S_f\) to the momentum balance with the sign determined by the adopted residual convention. For a Newtonian fluid,
# \[
# \boldsymbol\tau
# =2\mu\mathbf S-\frac23\mu(\nabla\cdot\mathbf u)\mathbf I,
# \qquad
# \mathbf S=\frac12(\nabla\mathbf u+\nabla\mathbf u^T).
# \]
#
# For a transported scalar,
# \[
# F_{c,f}=\dot m_f\phi_f.
# \]
# The same face mass flux must be used consistently by all transported equations sharing that flux.
#
# ## 5.3 Face interpolation
#
# Linear interpolation is
# \[
# \phi_f=(1-w_f)\phi_P+w_f\phi_N,
# \]
# while a reconstructed face value has the generic form
# \[
# \phi_f=\phi_P+\nabla\phi_P\cdot(\mathbf x_f-\mathbf C_P)+O(h^2).
# \]
# The face reconstruction contract is deliberately independent of the gradient implementation: Green–Gauss, LS and WLS gradients must be replaceable without changing the definition of the face-flux interface.
#
# ## 5.4 Upwind, central and high-resolution convection
#
# Upwind selects the upwind state:
# \[
# \phi_f=
# \begin{cases}
# \phi_P,&\dot m_f\ge0,\\
# \phi_N,&\dot m_f<0.
# \end{cases}
# \]
# Central interpolation is second-order only under the geometric/solution assumptions for which its truncation-error derivation holds; arbitrary polyhedral meshes require explicit verification rather than a blanket order claim.
#
# A high-resolution flux can be expressed schematically as
# \[
# \phi_f=\phi_f^{LO}+\psi(r)\left(\phi_f^{HO}-\phi_f^{LO}\right).
# \]
# The limiter, smoothness ratio, boundary treatment and multidimensional extension are part of the scheme definition.
#
# ## 5.5 Diffusive flux
#
# \[
# F_{d,f}=-\Gamma_f\nabla\phi_f\cdot\mathbf S_f.
# \]
# Decompose
# \[
# \mathbf S_f=\mathbf S_f^{\parallel}+\mathbf S_f^{\perp}
# \]
# to separate the centre-to-centre contribution from the non-orthogonal correction. The correction depends on reconstructed gradients and therefore couples N3 diffusion behaviour to N2 gradient accuracy.
#
# ## 5.6 Internal-face antisymmetry
#
# For an internal face,
# \[
# \mathbf S_{Nf}=-\mathbf S_{Pf},
# \qquad
# F_{Nf}=-F_{Pf}.
# \]
# The implementation should assemble a single physical face flux and apply opposite algebraic signs rather than independently recomputing two nominally equivalent fluxes.
#
# This invariant is stronger than global convergence: a converged solution can still hide a local assembly defect.
#
# ## 5.7 Boundary fluxes
#
# A prescribed normal flux satisfies
# \[
# F_\Gamma=q_nA_f.
# \]
# Dirichlet conditions determine a boundary state/value and therefore affect both convection and diffusion fluxes. Robin conditions combine value and flux, e.g.
# \[
# -k\nabla T\cdot\mathbf n=h(T_s-T_\infty).
# \]
# Boundary signs must follow the owner-cell outward-normal convention.
#
# ## 5.8 Boundedness and TVD
#
# Pairwise bounded interpolation satisfies
# \[
# \min(\phi_P,\phi_N)\le\phi_f\le\max(\phi_P,\phi_N).
# \]
# For 1-D linear advection,
# \[
# C=\frac{u\Delta t}{\Delta x},
# \qquad
# r=\frac{\phi_i-\phi_{i-1}}{\phi_{i+1}-\phi_i}.
# \]
# A TVD limiter must be analysed together with the underlying temporal/spatial discretisation and boundary treatment. Satisfying a limiter formula in isolation is not a proof that the multidimensional CFDX scheme is globally bounded.
#
# ## 5.9 Algebraic linearisation
#
# A linearised scalar equation has
# \[
# a_P\phi_P+\sum_Na_{PN}\phi_N=b_P.
# \]
# Explicit, implicit and deferred-correction pieces produce different coefficient structures. The method documentation must therefore state what is assembled implicitly, what is lagged, and what enters the source term.
#
# ## 5.10 Failure modes
#
# Important failure signatures include:
#
# - owner/neighbour sign reversal;
# - inconsistent mass and scalar face fluxes;
# - non-conservative boundary signs;
# - unbounded high-order reconstruction;
# - excessive numerical diffusion;
# - oscillations on steep gradients;
# - loss of accuracy from poor gradients;
# - non-orthogonal correction instability;
# - inconsistent pressure/velocity face fluxes.
#
# Diagnostics should identify the first failed invariant before changing solver relaxation or tolerances.
#
# ## 5.11 CFDX implementation and evidence
#
# Core interfaces: [flux.h](../../../src/cfdx/core/numerics/flux.h), [convection.h](../../../src/cfdx/core/numerics/convection.h), [interpolation.h](../../../src/cfdx/core/numerics/interpolation.h), [divergence.h](../../../src/cfdx/core/numerics/divergence.h), [laplacian.h](../../../src/cfdx/core/numerics/laplacian.h).
#
# Verification references include [test_convection_scheme_verification.cpp](../../../tests/validation/test_convection_scheme_verification.cpp), [test_convection_3d_verification.cpp](../../../tests/validation/test_convection_3d_verification.cpp), [test_convection_blending.cpp](../../../tests/validation/test_convection_blending.cpp), plus conservation/boundedness campaigns.
#
# ## 5.12 Executable antisymmetry
# %%
import numpy as np
F = 2.5
assert np.isclose(F + (-F), 0.0)
