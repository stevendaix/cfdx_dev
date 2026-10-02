# %% [markdown]
# 12 — Radiation
#
# ## 12.1 Blackbody
#
# \[
# E_b=\sigma T^4,\qquad
# \sigma=5.670374419\times10^{-8}\;W\,m^{-2}K^{-4}.
# \]
# Grey diffuse emission:
# \[
# E=\varepsilon E_b.
# \]
#
# ## 12.2 Irradiation and radiosity
#
# \[
# J=E+(1-\varepsilon)G,
# \]
# \[
# q''=J-G.
# \]
# Thus
# \[
# q''=\varepsilon(E_b-G).
# \]
#
# ## 12.3 Two-surface exchange
#
# For black surfaces:
# \[
# Q_{1\rightarrow2}
# =
# A_1F_{12}\sigma(T_1^4-T_2^4).
# \]
# For grey surfaces, the resistance-network form is
# \[
# Q=
# \frac{\sigma(T_1^4-T_2^4)}
# {\frac{1-\varepsilon_1}{A_1\varepsilon_1}
# +\frac1{A_1F_{12}}
# +\frac{1-\varepsilon_2}{A_2\varepsilon_2}}.
# \]
#
# ## 12.4 View factors
#
# \[
# \sum_jF_{ij}=1,
# \qquad
# A_iF_{ij}=A_jF_{ji},
# \qquad
# 0\le F_{ij}\le1.
# \]
# Closure, reciprocity and bounds are independent geometry checks.
#
# ## 12.5 S2S radiosity system
#
# \[
# G_i=\sum_jF_{ij}J_j,
# \]
# therefore
# \[
# \left[I-\operatorname{diag}(1-\varepsilon)F\right]J
# =\varepsilon\sigma T^4.
# \]
# Net surface power:
# \[
# Q_i=A_i(J_i-G_i).
# \]
#
# ## 12.6 P1
#
# A common participating-media P1 equation is
# \[
# -\nabla\cdot\left(\frac{1}{3\beta}\nabla G\right)
# +\beta G
# =4\beta\sigma T^4
# \]
# under a particular grey, isotropic convention. The exact implementation includes absorption/scattering and source conventions that must be checked against the source.
#
# ## 12.7 DOM
#
# The radiative transfer equation:
# \[
# \mathbf s\cdot\nabla I+\beta I
# =
# \kappa I_b+
# \sigma_s\int_{4\pi}\Phi(\mathbf s,\mathbf s')I(\mathbf s')d\Omega'.
# \]
# Discrete ordinates:
# \[
# \mathbf s_m\cdot\nabla I_m+\beta I_m
# =
# \kappa I_b+
# \sigma_s\sum_nw_n\Phi_{mn}I_n.
# \]
# Angular quadrature must reproduce declared moments:
# \[
# \sum_mw_m,\quad
# \sum_mw_m\mathbf s_m,\quad
# \sum_mw_m\mathbf s_m\otimes\mathbf s_m.
# \]
#
# ## 12.8 Participating media
#
# Absorption \(\kappa\), scattering \(\sigma_s\) and extinction
# \[
# \beta=\kappa+\sigma_s
# \]
# define optical thickness
# \[
# \tau=\beta L.
# \]
# Thin/intermediate/thick regimes should be tested separately.
#
# ## 12.9 Rosseland diffusion
#
# In optically thick media:
# \[
# \mathbf q_{rad}
# =-\frac{16\sigma T^3}{3\kappa_R}\nabla T.
# \]
# Hence an effective radiative conductivity is
# \[
# k_{rad}=\frac{16\sigma T^3}{3\kappa_R}.
# \]
# This is a diffusion approximation, not a replacement for general DOM/S2S physics.
#
# ## 12.10 Nonlinear radiation-energy coupling
#
# For a volumetric radiative loss \(q_{rad}(T)\), linearisation gives
# \[
# S(T)\approx S(T_*)+S_p(T-T_*),
# \qquad
# S_p=\frac{dS}{dT}(T_*).
# \]
# For \(q_{rad}=4\kappa\sigma T^4\):
# \[
# \frac{dq_{rad}}{dT}=16\kappa\sigma T^3.
# \]
# Correct source sign and linearisation are essential to energy stability.
#
# ## 12.11 Verification
#
# \[
# blackbody\rightarrow grey\ wall\rightarrow
# view\ factor\ closure/reciprocity\rightarrow
# S2S\ exchange\rightarrow
# P1\ equilibrium\rightarrow
# DOM\ angular\ moments\rightarrow
# coupled\ thermal\ energy.
# \]
#
# ## 12.12 CFDX implementation
#
# See [radiation.h](../../../src/cfdx/physics/radiation.h), [radiation_models.h](../../../src/cfdx/physics/radiation_models.h), [radiation_s2s.h](../../../src/cfdx/physics/radiation_s2s.h), [radiation_solver.h](../../../src/cfdx/physics/radiation_solver.h) and [radiation_streaming.h](../../../src/cfdx/physics/radiation_streaming.h).
#
# Tests: [test_radiation_vv.cpp](../../../tests/validation/test_radiation_vv.cpp), [test_s2s_radiation_vv.cpp](../../../tests/validation/test_s2s_radiation_vv.cpp), [test_phase12_radiation_vv.cpp](../../../tests/validation/test_phase12_radiation_vv.cpp), plus thermal/radiation regression.
#
# ## 12.13 Executable Stefan–Boltzmann check
# %%
import numpy as np
sigma=5.670374419e-8
assert np.isclose(sigma*300.0**4,459.300326939)
