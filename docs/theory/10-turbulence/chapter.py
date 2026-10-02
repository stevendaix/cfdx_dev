# %% [markdown]
# 10 — Turbulence
#
# This chapter separates the exact averaged/filter equations from closure assumptions. A turbulence model is not a numerical discretisation: the model equations, wall-distance definition, boundary conditions and discretisation must all be verified separately.
#
# ## 10.1 Reynolds decomposition
#
# \[
# u_i=\overline u_i+u_i',\qquad \overline{u_i'}=0.
# \]
# Averaging the nonlinear term produces the Reynolds stress:
# \[
# -\rho\overline{u_i'u_j'}.
# \]
# The closure problem is to express this unknown tensor using resolved mean quantities and additional model variables.
#
# ## 10.2 Reynolds stress and \(k\)
#
# \[
# k=\frac12\overline{u_i'u_i'}.
# \]
# The Boussinesq hypothesis writes
# \[
# -\rho\overline{u_i'u_j'}
# =2\mu_tS_{ij}-\frac23\rho k\delta_{ij}.
# \]
# This is a constitutive closure, not an identity.
#
# ## 10.3 Generic RANS transport
#
# A transported turbulence variable \(q\) has the form
# \[
# \frac{\partial(\rho q)}{\partial t}
# +\nabla\cdot(\rho\mathbf uq)
# =
# \nabla\cdot(D_q\nabla q)+P_q-D_q^{loss}+S_q.
# \]
# Every production/destruction/cross-diffusion term must be dimensionally consistent and sign-consistent with the implemented model.
#
# ## 10.4 \(k\)-\(\varepsilon\)
#
# A standard two-equation form is
# \[
# \frac{Dk}{Dt}=P_k-\varepsilon+
# \nabla\cdot\left[
# \left(\nu+\frac{\nu_t}{\sigma_k}\right)\nabla k
# \right],
# \]
# \[
# \frac{D\varepsilon}{Dt}
# =C_{\varepsilon1}\frac{\varepsilon}{k}P_k
# -C_{\varepsilon2}\frac{\varepsilon^2}{k}
# +\nabla\cdot\left[
# \left(\nu+\frac{\nu_t}{\sigma_\varepsilon}\right)\nabla\varepsilon
# \right].
# \]
# Exact constants and any low-Re/compressibility corrections are implementation-specific.
#
# ## 10.5 \(k\)-\(\omega\) and SST
#
# A generic \(k\)-\(\omega\) model contains
# \[
# \frac{Dk}{Dt}=P_k-\beta^*k\omega+\nabla\cdot(D_k\nabla k),
# \]
# \[
# \frac{D\omega}{Dt}
# =\alpha\frac{\omega}{k}P_k-\beta\omega^2
# +\nabla\cdot(D_\omega\nabla\omega)+CD.
# \]
# SST blends near-wall \(k\)-\(\omega\) and outer-layer behaviour:
# \[
# \Phi=F_1\Phi_1+(1-F_1)\Phi_2.
# \]
# The cross-diffusion and eddy-viscosity limiter are part of the model contract.
#
# ## 10.6 Spalart–Allmaras
#
# SA transports a modified viscosity variable \(\tilde\nu\):
# \[
# \frac{D\tilde\nu}{Dt}
# =P_{\tilde\nu}-D_{\tilde\nu}
# +\frac1{\sigma}\left[
# \nabla\cdot((\nu+\tilde\nu)\nabla\tilde\nu)
# +C_{b2}|\nabla\tilde\nu|^2
# \right].
# \]
# Production/destruction functions depend on wall distance and model functions. Wall-distance quality is therefore a direct model dependency.
#
# ## 10.7 LES
#
# Filtering gives
# \[
# \bar u_i=G*u_i,\qquad
# \tau_{ij}^{sgs}=\overline{u_iu_j}-\bar u_i\bar u_j.
# \]
# A Smagorinsky-type model uses
# \[
# \nu_t=(C_s\Delta)^2|S|,
# \qquad
# |S|=\sqrt{2S_{ij}S_{ij}}.
# \]
# WALE and dynamic models use different invariants and coefficients.
#
# ## 10.8 DES/DDES/IDDES
#
# Hybrid models replace or modify the turbulence length scale using grid and wall-distance information. Therefore model qualification must include grid sensitivity and shielding/wall-treatment cases; a generic LES result cannot qualify a hybrid model.
#
# ## 10.9 Boundary conditions and positivity
#
# Turbulence variables often have physical admissibility constraints:
# \[
# k\ge0,\qquad\varepsilon>0,\qquad\omega>0,\qquad\nu_t\ge0.
# \]
# Boundary conditions must be compatible with the wall treatment and far-field/inlet specification.
#
# ## 10.10 Verification ladder
#
# Verify individual algebraic terms, manufactured transport equations, positivity, wall-distance coupling, model limits and benchmark flows. Validation against experiments is a separate stage.
#
# ## 10.11 CFDX implementation
#
# See [turbulence.h](../../../src/cfdx/physics/turbulence.h), [turbulence_models.h](../../../src/cfdx/physics/turbulence_models.h), [turbulence_transport.h](../../../src/cfdx/physics/turbulence_transport.h), [turbulence_solver.h](../../../src/cfdx/physics/turbulence_solver.h), [spalart_allmaras.h](../../../src/cfdx/physics/spalart_allmaras.h) and [wall_distance.h](../../../src/cfdx/physics/wall_distance.h).
#
# Tests include [test_turbulence_qualification_equations.cpp](../../../tests/unit/test_turbulence_qualification_equations.cpp), [test_sst_blending.cpp](../../../tests/unit/test_sst_blending.cpp) and the turbulence hardening/qualification matrix.
#
# ## 10.12 Executable Reynolds decomposition check
# %%
import numpy as np
uprime=np.array([-1.,1.,0.,2.,-2.])
assert np.isclose(np.mean(uprime),0.0)
