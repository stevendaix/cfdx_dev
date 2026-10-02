# %% [markdown]
# 12 — Radiation
#
# ## 12.1 Black-body emission
#
# \[
# E_b=\sigma T^4,\qquad
# \sigma=5.670374419\times10^{-8}\;W\,m^{-2}K^{-4}.
# \]
# A grey diffuse surface emits
# \[
# E=\varepsilon\sigma T^4.
# \]
#
# ## 12.2 Irradiation and radiosity
#
# Let \(G\) be irradiation and \(J\) radiosity:
# \[
# J=\varepsilon\sigma T^4+(1-\varepsilon)G.
# \]
# Net outward radiative flux:
# \[
# q''=J-G.
# \]
# For opaque grey surfaces, \(\alpha+\rho=1\) and, at thermal equilibrium, Kirchhoff's law gives \(\alpha=\varepsilon\).
#
# ## 12.3 View factors
#
# The view factor \(F_{ij}\) is the fraction of radiation leaving surface \(i\) that reaches \(j\). For a closed enclosure:
# \[
# \sum_jF_{ij}=1.
# \]
# Reciprocity:
# \[
# A_iF_{ij}=A_jF_{ji}.
# \]
# These are geometric identities and therefore independent S2S verification targets.
#
# ## 12.4 S2S radiosity system
#
# Irradiation is
# \[
# G_i=\sum_jF_{ij}J_j.
# \]
# Substitution gives
# \[
# J_i=\varepsilon_i\sigma T_i^4+
# (1-\varepsilon_i)\sum_jF_{ij}J_j.
# \]
# Hence
# \[
# \left[I-\operatorname{diag}(1-\varepsilon)F\right]J
# =\varepsilon\sigma T^4.
# \]
# The surface net heat rate is \(Q_i=A_i(J_i-G_i)\).
#
# ## 12.5 Radiation conservation
#
# For an isolated enclosure, the net exchange satisfies
# \[
# \sum_iQ_i=0
# \]
# when there are no volumetric sources/sinks and the view-factor system is consistent. A non-zero defect can originate from geometry, reciprocity, matrix solution or boundary data.
#
# ## 12.6 P1 approximation
#
# A common grey, isotropic-scattering P1 form is
# \[
# -\nabla\cdot\left(\frac{1}{3\beta}\nabla G\right)
# +\beta G
# =4\beta\sigma T^4,
# \]
# where \(\beta\) is an extinction coefficient under the adopted convention. Absorption and scattering alter the exact coefficients/source terms.
#
# ## 12.7 DOM
#
# The radiative transfer equation may be written
# \[
# \mathbf s\cdot\nabla I(\mathbf x,\mathbf s)
# +\beta I
# =\kappa I_b+\sigma_s\int_{4\pi}\Phi(\mathbf s,\mathbf s')I(\mathbf s')\,d\Omega'.
# \]
# Discrete ordinates replace the angular integral by directions \(\mathbf s_m\) and weights \(w_m\):
# \[
# \mathbf s_m\cdot\nabla I_m+\beta I_m
# =\kappa I_b+\sigma_s\sum_n w_n\Phi_{mn}I_n.
# \]
# Angular quadrature accuracy and ray/discretisation effects must be part of verification.
#
# ## 12.8 Coupling to energy
#
# Radiation contributes an energy source/sink:
# \[
# S_T=-\nabla\cdot\mathbf q_{rad}
# \]
# in a volume formulation, or an equivalent surface heat flux at boundaries. The sign must match the thermal equation convention.
#
# ## 12.9 Verification
#
# S2S: view-factor closure → reciprocity → analytical two-surface exchange → enclosure energy conservation.
#
# P1/DOM: manufactured optical transport → limiting/positivity → angular/refinement studies → coupled thermal benchmark.
#
# ## 12.10 CFDX implementation
#
# See [radiation.h](../../../src/cfdx/physics/radiation.h), [radiation_models.h](../../../src/cfdx/physics/radiation_models.h), [radiation_s2s.h](../../../src/cfdx/physics/radiation_s2s.h), [radiation_solver.h](../../../src/cfdx/physics/radiation_solver.h) and [radiation_streaming.h](../../../src/cfdx/physics/radiation_streaming.h).
#
# Tests: [test_radiation_vv.cpp](../../../tests/validation/test_radiation_vv.cpp), [test_s2s_radiation_vv.cpp](../../../tests/validation/test_s2s_radiation_vv.cpp), [test_phase12_radiation_vv.cpp](../../../tests/validation/test_phase12_radiation_vv.cpp).
#
# ## 12.11 Executable Stefan–Boltzmann check
# %%
import numpy as np
sigma=5.670374419e-8
T=300.0
assert np.isclose(sigma*T**4,459.300326939)
