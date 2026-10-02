# %% [markdown]
# 11 — Heat Transfer
#
# ## 11.1 Energy balance
#
# A simple constant-property temperature equation is
# \[
# \rho c_p\frac{DT}{Dt}
# =\nabla\cdot(k\nabla T)+S_T.
# \]
# More generally, CFDX may solve an enthalpy/total-energy formulation:
# \[
# \frac{\partial(\rho h)}{\partial t}
# +\nabla\cdot(\rho\mathbf uh)
# =\nabla\cdot(k\nabla T)+S_h+\text{pressure/work/species terms},
# \]
# depending on the physics model. The exact solved variable must therefore be identified from the source.
#
# ## 11.2 Fourier law
#
# \[
# \mathbf q=-k\nabla T,\qquad
# q_n=-k\nabla T\cdot\mathbf n.
# \]
# \(k\,[W\,m^{-1}K^{-1}]\), so \(\mathbf q\,[W\,m^{-2}]\).
#
# ## 11.3 Thermal diffusivity
#
# \[
# \alpha=\frac{k}{\rho c_p}.
# \]
# The transient conduction equation becomes
# \[
# \frac{\partial T}{\partial t}=\alpha\nabla^2T
# \]
# for constant properties and no source.
#
# ## 11.4 Convection
#
# At a fluid wall:
# \[
# q''=h(T_s-T_\infty).
# \]
# The boundary condition is
# \[
# -k\nabla T\cdot\mathbf n=h(T_s-T_\infty).
# \]
# The sign depends on the outward-normal convention and must match the heat-flux definition used in the code.
#
# ## 11.5 Thermal boundary conditions
#
# Dirichlet:
# \[
# T=T_b.
# \]
# Neumann:
# \[
# -k\nabla T\cdot\mathbf n=q_b''.
# \]
# Robin:
# \[
# -k\nabla T\cdot\mathbf n=h(T-T_\infty).
# \]
# Radiation can add
# \[
# q''_{rad}=\varepsilon\sigma(T_s^4-T_{sur}^4).
# \]
#
# ## 11.6 Conjugate heat transfer
#
# At an ideal fluid-solid interface:
# \[
# T_f=T_s,
# \]
# \[
# \mathbf q_f\cdot\mathbf n_f+
# \mathbf q_s\cdot\mathbf n_s=0.
# \]
# Temperature continuity and flux conservation are independent checks.
#
# ## 11.7 Dimensionless groups
#
# \[
# Pe=\frac{UL}{\alpha},\qquad
# Re=\frac{\rho UL}{\mu},\qquad
# Pr=\frac{\nu}{\alpha},
# \]
# \[
# Nu=\frac{hL}{k},\qquad
# Bi=\frac{hL_c}{k_s},\qquad
# Fo=\frac{\alpha t}{L_c^2}.
# \]
# These groups identify the relative importance of transport mechanisms and provide benchmark-independent scaling checks.
#
# ## 11.8 Analytical benchmarks
#
# 1-D steady conduction with no source:
# \[
# \frac{d^2T}{dx^2}=0
# \Rightarrow T(x)=T_0+\frac{T_L-T_0}{L}x.
# \]
# Uniform volumetric generation:
# \[
# k\frac{d^2T}{dx^2}+q'''=0
# \Rightarrow
# T(x)=T_0+\frac{C_1}{ }x-\frac{q'''}{2k}x^2,
# \]
# with constants fixed by the two boundary conditions.
#
# ## 11.9 Verification
#
# The ladder is analytical 1-D conduction → MMS → transient conduction → interface conservation → coupled thermal-radiation/flow benchmarks. Thermal residual convergence alone is insufficient for interface conservation.
#
# ## 11.10 CFDX implementation
#
# See [thermal.h](../../../src/cfdx/physics/thermal.h), [energy_solver.h](../../../src/cfdx/physics/energy_solver.h), [cht_solver.h](../../../src/cfdx/physics/cht_solver.h) and the conductivity/thermophysical model families.
#
# Verification includes [test_thermal_vv.cpp](../../../tests/validation/test_thermal_vv.cpp), [test_cht_validation.cpp](../../../tests/validation/test_cht_validation.cpp) and the thermal-radiation matrix.
#
# ## 11.11 Executable property check
# %%
import numpy as np
k,rho,cp=10.,1000.,4200.
alpha=k/(rho*cp)
assert np.isclose(alpha, 10./4.2e6)
