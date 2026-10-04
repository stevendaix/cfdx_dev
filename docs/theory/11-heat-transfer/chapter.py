# %% [markdown]
# 11 — Heat Transfer
#
# ## 11.1 Temperature and enthalpy formulations
#
# Constant-property temperature transport:
# $
# \rho c_p\frac{DT}{Dt}
# =\nabla\cdot(k\nabla T)+S_T.
# $
# Enthalpy satisfies
# $
# \rho\frac{Dh}{Dt}
# =\frac{Dp}{Dt}+\nabla\cdot(k\nabla T)+\Phi+S_h
# $
# in a representative single-phase formulation, where $\Phi$ is viscous dissipation. The exact energy formulation depends on the selected thermophysical model.
#
# ## 11.2 Fourier law
#
# $
# \mathbf q=-k\nabla T.
# $
# For a two-point owner-neighbour approximation:
# $
# q_n\approx-k\frac{T_N-T_P}{d_{PN}}.
# $
# CFDX's thermal primitive exposes this relation explicitly.
#
# ## 11.3 Thermal diffusivity
#
# $
# \alpha=\frac{k}{\rho c_p}.
# $
# For constant properties and no source:
# $
# \partial_tT=\alpha\nabla^2T.
# $
#
# ## 11.4 Convection
#
# At a convective boundary:
# $
# -k\nabla T\cdot\mathbf n=h(T_s-T_\infty).
# $
# The convective heat-transfer coefficient is a model/input quantity unless derived by a coupled fluid solution.
#
# ## 11.5 Radiation boundary
#
# $
# q''_{rad}=\varepsilon\sigma(T_s^4-T_{sur}^4).
# $
# Thermal/radiation coupling is nonlinear because
# $
# \frac{dq''_{rad}}{dT_s}=4\varepsilon\sigma T_s^3.
# $
#
# ## 11.6 Sources and linearisation
#
# A source can be linearised about $T_0$:
# $
# S(T)\approx S(T_0)+
# \left.\frac{dS}{dT}\right|_{T_0}(T-T_0)
# =S_u+S_pT.
# $
# Correct reconstruction requires
# $
# S_u+S_pT_0=S(T_0).
# $
# This identity is a useful code-verification test.
#
# ## 11.7 CHT interface
#
# Perfect contact:
# $
# T_1=T_2,\qquad
# q_1''+q_2''=0.
# $
# With a finite contact conductance $G$:
# $
# q''=G(T_1-T_2).
# $
# Series thermal resistance:
# $
# G=\frac{A}{d_1/k_1+d_2/k_2}.
# $
#
# ## 11.8 Dimensionless groups
#
# $
# Pe=\frac{UL}{\alpha},\quad
# Pr=\frac{\nu}{\alpha},\quad
# Nu=\frac{hL}{k},
# $
# $
# Bi=\frac{hL_c}{k_s},\quad
# Fo=\frac{\alpha t}{L_c^2}.
# $
#
# ## 11.9 Analytical benchmarks
#
# 1-D steady conduction:
# $
# T(x)=T_0+\frac{T_L-T_0}{L}x.
# $
# Uniform generation:
# $
# kT''+q'''=0
# \Rightarrow
# T(x)=-\frac{q'''}{2k}x^2+C_1x+C_2.
# $
# The constants follow from the two boundary conditions.
#
# Fully developed heated-channel balances can be checked through
# $
# \dot m c_p(T_b-T_{b,0})=Q_{wall}.
# $
#
# ## 11.10 Thermophysical closures
#
# CFDX supports density/enthalpy/viscosity/conductivity model families including constant and temperature-dependent laws and ideal-gas density. A general EOS relation is
# $
# \rho=\rho(p,T,\ldots),
# $
# while ideal gas gives
# $
# \rho=\frac{p}{RT}.
# $
# Property-model verification is separate from solver-level thermal V&V.
#
# ## 11.11 Verification ladder
#
# $
# property\ law\rightarrow
# conduction\ kernel\rightarrow
# source\ linearisation\rightarrow
# transient\ MMS\rightarrow
# CHT\ interface\rightarrow
# radiation/energy\ coupling.
# $
#
# ## 11.12 CFDX implementation
#
# Thermal primitive: [thermal.h](../../../src/cfdx/physics/thermal.h). Energy/CHT: [energy_solver.h](../../../src/cfdx/physics/energy_solver.h), [cht_solver.h](../../../src/cfdx/physics/cht_solver.h). Thermophysical closures: [thermophysical_models.h](../../../src/cfdx/physics/thermophysical_models.h).
#
# Tests include [test_thermal_vv.cpp](../../../tests/validation/test_thermal_vv.cpp), [test_cht_validation.cpp](../../../tests/validation/test_cht_validation.cpp), [test_thermophysical_models_vv.cpp](../../../tests/validation/test_thermophysical_models_vv.cpp) and the thermal/radiation regression matrix.
#
# ## 11.13 Executable diffusivity check
# %%
import numpy as np
k,rho,cp=10.,1000.,4200.
assert np.isclose(k/(rho*cp),10./4.2e6)
