# %% [markdown]
# 13 — Multiphysics Coupling
#
# ## 13.1 Coupled problem
#
# Let
# [
# y=(u,p,T,Y_1,ldots,Y_N).
# ]
# The complete model is
# [
# F(y)=0.
# ]
# Diagonal Jacobian blocks represent self-physics; off-diagonal blocks represent feedback.
#
# ## 13.2 Newton and fixed-point formulations
#
# Newton linearisation gives
# [
# F(y^k+delta y)approx F(y^k)+J(y^k)delta y,
# ]
# [
# J(y^k)delta y=-F(y^k),
# qquad
# y^{k+1}=y^k+omegadelta y.
# ]
# A segregated scheme instead applies
# [
# y^{k+1}=S(y^k),
# ]
# or relaxed
# [
# y^{k+1}leftarrow(1-omega)y^k+omega S(y^k).
# ]
# Local fixed-point convergence requires
# [
# ho(S'(y^*))<1.
# ]
#
# ## 13.3 Block Jacobian
#
# [
# J=
# egin{bmatrix}
# J_{uu}&J_{up}&J_{uT}&J_{uY}\
# J_{pu}&J_{pp}&J_{pT}&J_{pY}\
# J_{Tu}&J_{Tp}&J_{TT}&J_{TY}\
# J_{Yu}&J_{Yp}&J_{YT}&J_{YY}
# end{bmatrix}.
# ]
# This structure explains why monolithic block preconditioners can outperform completely segregated iterations when feedback is strong.
#
# ## 13.4 Interface conservation
#
# For an interface (Gamma),
# [
# R_Gamma=
# sum_{finGamma}
# left(q_f^{(1)}+q_f^{(2)}ight).
# ]
# Correct coupling requires (R_Gamma	o0). For perfect thermal contact,
# [
# T_1=T_2,qquad q_1''+q_2''=0.
# ]
# For finite contact conductance,
# [
# q''=G(T_1-T_2).
# ]
#
# ## 13.5 Thermal-fluid coupling
#
# [
# ho=ho(T),qquadmu=mu(T),
# ]
# modifies momentum and continuity. Buoyancy can be represented by
# [
# mathbf f_b=homathbf g
# ]
# or, under the Boussinesq approximation,
# [
# homathbf gapprox
# ho_0[1-eta(T-T_0)]mathbf g.
# ]
# The approximation is valid only within its stated density-variation assumptions.
#
# ## 13.6 Radiation and turbulence feedback
#
# Radiation enters the energy equation through
# [
# S_{rad}=-
ablacdotmathbf q_{rad}
# ]
# or a boundary flux, with the nonlinear blackbody dependence
# [
# q_{rad}propto T^4.
# ]
# Turbulence modifies transport:
# [
# mu_{eff}=mu+mu_t,qquad
# k_{eff}=k+k_t.
# ]
# These dependencies create nonlinear feedback loops between momentum, turbulence, energy and radiation.
#
# ## 13.7 Dimensionless coupling indicators
#
# Useful indicators include
# [
# Re=rac{ho UL}{mu},qquad
# Pr=rac{mu c_p}{k},qquad
# Ra=rac{getaDelta TL^3}{
ualpha},
# ]
# and, for conjugate heat transfer,
# [
# Bi=rac{hL_c}{k_s}.
# ]
# They help identify whether a coupling is weak, diffusion dominated or strongly interactive; they do not replace a verification case.
#
# ## 13.8 Convergence hierarchy
#
# Monitor separately
# [
# |F_u|,|F_p|,|F_T|,|F_Y|,
# ]
# interface defects, global conservation defects, linear residuals and physically relevant integral quantities. A small linear residual proves only that the current linearised system was solved accurately.
#
# ## 13.9 Verification ladder
#
# [
# 	ext{single physics}
# ightarrow	ext{one-way coupling}
# ightarrow	ext{two-way coupling}
# ightarrow	ext{fully coupled}.
# ]
# Every coupling should have a limiting case in which its influence vanishes and an independent interface/conservation diagnostic.
#
# ## 13.10 CFDX implementation
#
# Coupling is distributed through [finite_volume_transport.h](../../../src/cfdx/physics/finite_volume_transport.h), [cht_solver.h](../../../src/cfdx/physics/cht_solver.h), [radiation_solver.h](../../../src/cfdx/physics/radiation_solver.h), [steady_incompressible_solver.h](../../../src/cfdx/physics/steady_incompressible_solver.h) and the associated physics families.
#
# Tests include [test_level_c_coupled_verification.cpp](../../../tests/validation/test_level_c_coupled_verification.cpp) and CHT/radiation validation.
#
# ## 13.11 Executable Newton step
# %%
import numpy as np
J=np.array([[3.,1.],[1.,2.]])
F=np.array([1.,-1.])
delta=np.linalg.solve(J,-F)
assert np.allclose(J@delta,-F)
