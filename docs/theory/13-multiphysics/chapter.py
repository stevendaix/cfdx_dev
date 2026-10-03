# %% [markdown]
# # 13 — Multiphysics Coupling
#
# ## 13.1 Coupled nonlinear problem
#
# Let
# \[
# y=(\mathbf u,p,T,Y_1,\ldots,Y_N).
# \]
# The coupled model is
# \[
# F(y)=0.
# \]
# Diagonal Jacobian blocks describe self-physics; off-diagonal blocks describe feedback between physics.
#
# A coupling claim must therefore specify not only which models exist, but which feedback paths are actually assembled.
#
# ## 13.2 Newton and fixed-point formulations
#
# Newton linearisation gives
# \[
# J(y^k)\delta y=-F(y^k),
# \qquad
# y^{k+1}=y^k+\omega\delta y.
# \]
# A segregated fixed-point iteration is
# \[
# y^{k+1}=S(y^k),
# \]
# possibly relaxed as
# \[
# y^{k+1}\leftarrow(1-\omega)y^k+\omega S(y^k).
# \]
# Local fixed-point convergence requires the spectral radius
# \[
# \rho(S'(y^*))<1.
# \]
#
# ## 13.3 Block Jacobian and Schur structure
#
# \[
# J=
# \begin{bmatrix}
# J_{uu}&J_{up}&J_{uT}&J_{uY}\\
# J_{pu}&J_{pp}&J_{pT}&J_{pY}\\
# J_{Tu}&J_{Tp}&J_{TT}&J_{TY}\\
# J_{Yu}&J_{Yp}&J_{YT}&J_{YY}
# \end{bmatrix}.
# \]
# Strong off-diagonal coupling can make a fully or partially coupled block strategy preferable to completely segregated iterations. The numerical method must document which blocks are assembled exactly, approximated, lagged or omitted.
#
# ## 13.4 Interface conservation
#
# For an interface \(\Gamma\),
# \[
# R_\Gamma=\sum_{f\in\Gamma}\left(q_f^{(1)}+q_f^{(2)}\right).
# \]
# A conservative coupling requires \(R_\Gamma\to0\) independently of the field-value error.
#
# For perfect thermal contact,
# \[
# T_1=T_2,
# \qquad
# q_1''+q_2''=0.
# \]
# With finite contact conductance,
# \[
# q''=G(T_1-T_2).
# \]
# Both temperature jump and heat-flux conservation must be checked.
#
# ## 13.5 Thermal-fluid coupling
#
# Temperature-dependent properties give
# \[
# \rho=\rho(p,T),\qquad \mu=\mu(T),\qquad k=k(T).
# \]
# Under Boussinesq,
# \[
# \rho\mathbf g\approx
# \rho_0[1-\beta(T-T_0)]\mathbf g.
# \]
# The approximation assumes density variations are sufficiently small for the selected formulation; it must not be confused with a general variable-density equation of state.
#
# ## 13.6 Radiation and turbulence feedback
#
# Radiation enters through a volumetric source or boundary heat flux:
# \[
# S_{rad}=-\nabla\cdot\mathbf q_{rad},
# \]
# with nonlinear dependence such as \(q_{rad}\propto T^4\).
#
# Turbulence changes effective transport:
# \[
# \mu_{eff}=\mu+\mu_t,
# \qquad
# k_{eff}=k+k_t.
# \]
# Consequently, a fully coupled thermal-turbulent-radiative calculation contains feedback loops rather than independent post-processing steps.
#
# ## 13.7 Coupling strength and dimensionless groups
#
# \[
# Re=\frac{\rho UL}{\mu},
# \quad
# Pr=\frac{\mu c_p}{k},
# \quad
# Ra=\frac{g\beta\Delta T L^3}{\nu\alpha},
# \quad
# Bi=\frac{hL_c}{k_s}.
# \]
# These identify relevant regimes but do not constitute verification criteria by themselves.
#
# ## 13.8 Convergence hierarchy
#
# Monitor independently:
# \[
# \|F_u\|,\ \|F_p\|,\ \|F_T\|,\ \|F_Y\|,
# \]
# interface defects, global conservation defects, linear residuals and physical quantities of interest.
#
# A small linear residual proves only that the current linearised system was solved accurately. It does not prove convergence of the nonlinear coupled problem or correctness of the physical model.
#
# ## 13.9 Coupling verification ladder
#
# Use the sequence
# \[
# single\ physics
# \rightarrow one\!-!way\ coupling
# \rightarrow two\!-!way\ coupling
# \rightarrow fully\ coupled.
# \]
# Each stage requires a limiting case in which the coupling disappears and an independent interface/conservation diagnostic.
#
# Recommended examples:
#
# - fluid only → thermal as passive scalar;
# - thermal → Boussinesq feedback;
# - thermal → radiation boundary exchange;
# - turbulence → effective transport;
# - CHT → finite interface resistance;
# - fully coupled thermal/radiation/flow.
#
# ## 13.10 Failure modes
#
# Typical coupling failures include:
#
# - interface flux sign mismatch;
# - double-counted source terms;
# - stale property values;
# - inconsistent linearisation;
# - under-relaxation hiding a non-convergent fixed point;
# - converged individual fields with a non-converged interface balance;
# - incorrect Boussinesq model selection;
# - radiation/energy source sign errors.
#
# A coupling diagnostic must identify which block or interface first violates its contract.
#
# ## 13.11 CFDX implementation and evidence
#
# Coupling is distributed through [finite_volume_transport.h](../../../src/cfdx/physics/finite_volume_transport.h), [cht_solver.h](../../../src/cfdx/physics/cht_solver.h), [radiation_solver.h](../../../src/cfdx/physics/radiation_solver.h), [steady_incompressible_solver.h](../../../src/cfdx/physics/steady_incompressible_solver.h) and associated physics families.
#
# Tests include [test_level_c_coupled_verification.cpp](../../../tests/validation/test_level_c_coupled_verification.cpp) and the CHT/radiation validation matrix. The exact tests and source paths must be rechecked when implementation changes.
#
# ## 13.12 Executable Newton step
# %%
import numpy as np
J=np.array([[3.,1.],[1.,2.]])
F=np.array([1.,-1.])
delta=np.linalg.solve(J,-F)
assert np.allclose(J@delta,-F)
