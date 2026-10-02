# %% [markdown]
# Heat Transfer
#
# ## 11.1 Energy equation
# \[
# \rho c_p\frac{DT}{Dt}=\nabla\cdot(k\nabla T)+S_T.
# \]
# Thermal diffusivity:
# \[
# \alpha=\frac{k}{\rho c_p}.
# \]
#
# ## 11.2 Fourier law
# \[
# \mathbf q=-k\nabla T,\qquad
# q_n=-k\nabla T\cdot\mathbf n.
# \]
#
# ## 11.3 Boundary conditions
# Fixed temperature: \(T=T_b\).  
# Fixed flux:
# \[
# -k\nabla T\cdot\mathbf n=q_b''.
# \]
# Convection:
# \[
# -k\nabla T\cdot\mathbf n=h(T_s-T_\infty).
# \]
#
# ## 11.4 Conjugate heat transfer
# At a perfect interface:
# \[
# T_1=T_2,\qquad
# \mathbf q_1\cdot\mathbf n_1+\mathbf q_2\cdot\mathbf n_2=0.
# \]
# Temperature continuity and heat-flux conservation are separate requirements.
#
# ## 11.5 Dimensionless groups
# \[
# Pe=\frac{UL}{\alpha},\qquad
# Bi=\frac{hL_c}{k_s},\qquad
# Fo=\frac{\alpha t}{L_c^2}.
# \]
#
# ## 11.6 Verification
# One-dimensional conduction, manufactured solutions, transient analytical solutions and CHT interface balances provide progressively stronger evidence.
#
# ## 11.7 CFDX traceability
# src/cfdx/physics/thermal.h  
# src/cfdx/physics/energy_solver.h  
# src/cfdx/physics/cht_solver.h  
# src/cfdx/transport/conductivity/conductivity.h  
# Tests: tests/validation/test_thermal_vv.cpp, tests/validation/test_cht_validation.cpp, tests/validation/test_thermal_radiation_model_matrix.cpp
#
# %%
from __future__ import annotations
k,rho,cp=10.,1000.,4200.
assert np.isclose(k/(rho*cp),10./4.2e6)
