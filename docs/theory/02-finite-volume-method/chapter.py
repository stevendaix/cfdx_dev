# %% [markdown]
# 02 — Finite-Volume Method
#
# This chapter derives the finite-volume discretisation from the conservation law to the algebraic row. The key principle is that each control volume owns a balance: internal-face fluxes are shared and cancel exactly when the same numerical face flux is used with opposite orientation.
#
# ## 2.1 Generic transport equation
#
# We use
# \[
# \frac{\partial(\rho\phi)}{\partial t}
# +\nabla\cdot(\rho\mathbf u\phi)
# =\nabla\cdot(\Gamma\nabla\phi)+S_\phi .
# \]
# Here \(\rho\,[kg\,m^{-3}]\), \(\mathbf u\,[m\,s^{-1}]\), \(\Gamma\) is the diffusion coefficient and \(S_\phi\) is a volumetric source. The dimensions of \(\phi\) depend on the transported quantity.
#
# ## 2.2 Control-volume integration
#
# Integrating over \(V_P\) and applying Gauss' theorem:
# \[
# \frac{d}{dt}\int_{V_P}\rho\phi\,dV+
# \sum_f\int_{A_f}\rho\mathbf u\phi\cdot\mathbf n_f\,dA
# =
# \sum_f\int_{A_f}\Gamma\nabla\phi\cdot\mathbf n_f\,dA+
# \int_{V_P}S_\phi\,dV.
# \]
# Define \(\mathbf S_f=A_f\mathbf n_f\). For a closed polyhedron,
# \[
# \sum_f\mathbf S_f=0.
# \]
#
# ## 2.3 Spatial approximation
#
# With cell-centred values:
# \[
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f,\qquad
# F_{c,f}=\dot m_f\phi_f,
# \]
# and
# \[
# F_{d,f}=\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
# \]
# The sign convention must be stated once: here positive \(F_f\) is outward from owner cell \(P\). The neighbour sees the same physical flux with the opposite orientation.
#
# ## 2.4 Semi-discrete row
#
# For implicit Euler:
# \[
# \frac{(\rho\phi V)_P^{n+1}-(\rho\phi V)_P^n}{\Delta t}
# +\sum_fF_{c,f}^{n+1}
# =\sum_fF_{d,f}^{n+1}+S_PV_P .
# \]
# After linearisation:
# \[
# a_P\phi_P=\sum_Na_{PN}\phi_N+b_P.
# \]
# The coefficient row is not universal: it depends on convection scheme, diffusion reconstruction, source linearisation, temporal scheme and BC discretisation.
#
# ## 2.5 Diffusion on an orthogonal face
#
# Let \(\mathbf d_{PN}=\mathbf C_N-\mathbf C_P\). If \(\mathbf S_f\parallel\mathbf d_{PN}\),
# \[
# (\nabla\phi)_f\cdot\mathbf S_f
# \approx
# \frac{\phi_N-\phi_P}{|\mathbf d_{PN}|}
# \frac{\mathbf d_{PN}\cdot\mathbf S_f}{|\mathbf d_{PN}|}.
# \]
# More generally decompose
# \[
# \mathbf S_f=\mathbf S_f^\parallel+\mathbf S_f^\perp
# \]
# and reconstruct the correction using a verified gradient.
#
# ## 2.6 Source linearisation
#
# A common linearisation is
# \[
# S_\phi=S_C+S_P\phi_P,
# \]
# so that the diagonal contribution is modified by the chosen source treatment. Stability and boundedness depend on the signs and treatment of these terms; they must be verified rather than assumed.
#
# ## 2.7 Boundary faces
#
# A boundary face has no ordinary neighbour. Dirichlet, Neumann and Robin conditions therefore become algebraic contributions to the owner row. For a prescribed diffusive flux \(q_f\),
# \[
# F_{d,f}=q_fA_f.
# \]
# For a prescribed value, the face value and gradient reconstruction must be derived from the actual boundary discretisation.
#
# ## 2.8 Conservation and boundedness are different
#
# Conservation requires shared internal fluxes to be antisymmetric:
# \[
# F_{P,f}+F_{N,f}=0.
# \]
# Boundedness instead concerns admissible solution values, e.g.
# \[
# \min_{N\in\mathcal N(P)}\phi_N
# \le\phi_f\le
# \max_{N\in\mathcal N(P)}\phi_N .
# \]
# A scheme may be conservative but unbounded, or bounded but non-conservative if implemented incorrectly.
#
# ## 2.9 Algebraic properties
#
# For diffusion-dominated scalar equations, an M-matrix-like structure is often desirable:
# \[
# a_P>0,\qquad a_{PN}\le0
# \]
# under the relevant sign convention, together with appropriate diagonal dominance. Convection, non-orthogonal correction and source terms can invalidate these properties; therefore the property belongs to a declared scheme and mesh regime.
#
# ## 2.10 Polyhedral FVM
#
# The method does not require hexahedra. For arbitrary polyhedra the essential objects are oriented faces, owner/neighbour connectivity, face area vectors, cell volumes and consistent centres. The topological and geometric invariants are documented in the mesh chapter.
#
# ## 2.11 Verification sequence
#
# 1. closed-cell area-vector closure;
# 2. constant-field invariance;
# 3. exact cancellation of internal fluxes;
# 4. linear-field diffusion;
# 5. MMS;
# 6. mesh refinement and observed order;
# 7. coupled solver benchmarks.
#
# A successful global benchmark is not sufficient to diagnose an operator-level error.
#
# ## 2.12 CFDX implementation
#
# Core numerical operators are in [divergence.h](../../../src/cfdx/core/numerics/divergence.h), [laplacian.h](../../../src/cfdx/core/numerics/laplacian.h), [convection.h](../../../src/cfdx/core/numerics/convection.h), [source_term.h](../../../src/cfdx/core/numerics/source_term.h), and [flux.h](../../../src/cfdx/core/numerics/flux.h). The transport assembly is in [finite_volume_transport.h](../../../src/cfdx/physics/finite_volume_transport.h).
#
# ## 2.13 Executable analytical check
# %%
import numpy as np
S = np.array([[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.]])
assert np.allclose(S.sum(axis=0), 0.0)
# For a constant field, the convective contribution is zero when the mass fluxes close.
phi = 7.0
assert np.isclose(phi * S.sum(axis=0)[0], 0.0)


# %% [markdown]
# ## 2.14 Discrete conservation and algebraic assembly
#
# For a cell P, the conservative semi-discrete balance can be written
#
# \[
# \frac{d}{dt}(\rho_P\phi_PV_P)+\sum_f F_{c,f}-\sum_f F_{d,f}-S_PV_P=0.
# \]
#
# The assembly contract is that every internal-face contribution is generated once with owner orientation and once with the opposite neighbour orientation. The two contributions must represent the same physical numerical flux with opposite signs.
#
# For a linearised row,
#
# \[
# a_P\phi_P+\sum_N a_{PN}\phi_N=b_P.
# \]
#
# The coefficient pattern depends on the selected interpolation, diffusion reconstruction, source treatment, boundary condition and time scheme. It must therefore be inspected as an assembled operator rather than inferred from the continuous equation alone.
#
# ## 2.15 Verification ladder
#
# The recommended sequence is:
#
# 1. geometric face-vector closure;
# 2. constant-field invariance;
# 3. internal-face flux antisymmetry;
# 4. exact simple-field diffusion;
# 5. isolated operator MMS;
# 6. mesh convergence;
# 7. coupled benchmark.
#
# This ordering localises defects before nonlinear coupling can mask them.
#
# ## 2.16 Limitations and improvement paths
#
# Non-orthogonal and skewed meshes introduce reconstruction error. Higher-order convection introduces boundedness/oscillation trade-offs. Strong source terms can alter diagonal dominance. Future improvements should therefore be qualified per operator, mesh family and scheme rather than by a single global benchmark.
