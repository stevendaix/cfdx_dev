# %% [markdown]
# 13 — Multiphysics Coupling
#
# ## 13.1 Coupled problem
#
# Let
# \[
# y=(u,p,T,Y_1,\ldots,Y_N).
# \]
# The complete coupled model is
# \[
# F(y)=0.
# \]
# Each block contains one physical subsystem, while off-diagonal blocks represent feedback between subsystems.
#
# ## 13.2 Newton linearisation
#
# Around \(y^k\):
# \[
# F(y^k+\delta y)\approx F(y^k)+J(y^k)\delta y.
# \]
# Solve
# \[
# J(y^k)\delta y=-F(y^k),
# \qquad
# y^{k+1}=y^k+\delta y.
# \]
# A damped update uses
# \[
# y^{k+1}=y^k+\omega\delta y,\qquad0<\omega\le1.
# \]
#
# ## 13.3 Block Jacobian
#
# \[
# J=
# \begin{bmatrix}
# J_{uu}&J_{up}&J_{uT}&\cdots\\
# J_{pu}&J_{pp}&J_{pT}&\cdots\\
# J_{Tu}&J_{Tp}&J_{TT}&\cdots
# \end{bmatrix}.
# \]
# For example \(J_{uT}\) represents thermal feedback on momentum through density/viscosity/buoyancy, while \(J_{Tu}\) represents convective heat transport by velocity.
#
# ## 13.4 Segregated fixed-point coupling
#
# \[
# y^{k+1}=S(y^k).
# \]
# Relaxation:
# \[
# y^{k+1}\leftarrow(1-\omega)y^k+\omega S(y^k).
# \]
# If the local linearisation is \(S'(y^\*)\), convergence requires the spectral radius to satisfy
# \[
# \rho(S'(y^\*))<1.
# \]
# This explains why strongly coupled physics can require under-relaxation or stronger block coupling.
#
# ## 13.5 Interface coupling
#
# For an interface \(\Gamma\), define a defect
# \[
# R_\Gamma=
# \sum_{f\in\Gamma}
# \left(q_f^{(1)}+q_f^{(2)}\right).
# \]
# Correct coupling requires \(R_\Gamma\rightarrow0\) independently of field residuals.
#
# ## 13.6 Thermal-fluid coupling
#
# Density/viscosity dependence:
# \[
# \rho=\rho(T),\qquad\mu=\mu(T)
# \]
# changes momentum and continuity. Buoyancy may add
# \[
# \mathbf f_b=\rho\mathbf g
# \]
# or a Boussinesq approximation
# \[
# \rho\mathbf g\approx\rho_0[1-\beta(T-T_0)]\mathbf g.
# \]
# The exact model and sign convention must be explicit.
#
# ## 13.7 Radiation-thermal coupling
#
# Thermal balance receives
# \[
# S_{rad}=-\nabla\cdot\mathbf q_{rad}
# \]
# or a boundary radiative flux. Because \(q_{rad}\sim T^4\), the coupling is strongly nonlinear at high temperature.
#
# ## 13.8 Turbulence coupling
#
# Turbulence changes effective transport:
# \[
# \mu_{eff}=\mu+\mu_t,
# \]
# and often
# \[
# k_{eff}=k+k_t.
# \]
# Thus turbulence residuals, wall distance and thermal equations can form a coupled nonlinear loop.
#
# ## 13.9 Convergence hierarchy
#
# Monitor separately:
# \[
# \|F_u\|,\;\|F_p\|,\;\|F_T\|,\;\|F_Y\|,
# \]
# interface defects, conservation defects and linear-solver residuals. A small linear residual does not prove nonlinear convergence; a small nonlinear residual does not prove conservation if the residual itself is assembled incorrectly.
#
# ## 13.10 Verification
#
# Start with decoupled analytical limits, then one-way coupling, then two-way coupling, then full multiphysics. Every added coupling should have an independent invariant or limiting case.
#
# ## 13.11 CFDX implementation
#
# Coupling is distributed through [finite_volume_transport.h](../../../src/cfdx/physics/finite_volume_transport.h), [cht_solver.h](../../../src/cfdx/physics/cht_solver.h), [radiation_solver.h](../../../src/cfdx/physics/radiation_solver.h), [steady_incompressible_solver.h](../../../src/cfdx/physics/steady_incompressible_solver.h) and the associated physics families.
#
# Tests include [test_level_c_coupled_verification.cpp](../../../tests/validation/test_level_c_coupled_verification.cpp) and CHT/radiation validation.
#
# ## 13.12 Executable Newton step
# %%
import numpy as np
J=np.array([[3.,1.],[1.,2.]])
F=np.array([1.,-1.])
delta=np.linalg.solve(J,-F)
assert np.allclose(J@delta,-F)
