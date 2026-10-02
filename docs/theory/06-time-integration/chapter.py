# %% [markdown]
# 06 — Time Integration
#
# ## 6.1 Semi-discrete problem
#
# After spatial discretisation,
# [
# Mrac{dU}{dt}=R(U,t).
# ]
# (M) contains cell volumes, densities or other storage coefficients. (R) is the assembled physical residual. Temporal order can only be measured when spatial error and nonlinear/linear iteration error are smaller than the temporal error being studied.
#
# ## 6.2 Explicit Euler
#
# [
# U^{n+1}=U^n+Delta t,M^{-1}R(U^n).
# ]
# For (u'=lambda u),
# [
# G_E(z)=1+z,qquad z=lambdaDelta t.
# ]
# For the decay equation (u'=-lambda u), stability requires
# [
# |1-lambdaDelta t|le1
# quadLongrightarrowquad
# 0lelambdaDelta tle2.
# ]
#
# ## 6.3 Implicit Euler
#
# [
# rac{U^{n+1}-U^n}{Delta t}
# =R(U^{n+1}),
# ]
# and for (u'=-lambda u),
# [
# G_I=rac1{1+lambdaDelta t}.
# ]
# The method is first-order and A-stable for the linear test equation. It is strongly dissipative for large negative (z).
#
# ## 6.4 Crank–Nicolson
#
# [
# rac{U^{n+1}-U^n}{Delta t}
# =rac12[R(U^{n+1})+R(U^n)].
# ]
# For (u'=-lambda u),
# [
# G_{CN}=
# rac{1-rac12lambdaDelta t}
# {1+rac12lambdaDelta t}.
# ]
# It is second-order and A-stable but not L-stable; stiff components can therefore remain weakly damped.
#
# ## 6.5 BDF2
#
# Constant-step BDF2 is
# [
# rac{3U^{n+1}-4U^n+U^{n-1}}{2Delta t}
# =R(U^{n+1}).
# ]
# For variable steps, let
# [
# omega=rac{Delta t_n}{Delta t_{n-1}}.
# ]
# The derivative coefficients must be recomputed from the three actual time levels. A useful coefficient representation is
# [
# rac{dU}{dt}igg|_{n+1}
# approx
# rac{1}{Delta t_n}
# left(
# a_0U^{n+1}+a_1U^n+a_2U^{n-1}
# ight),
# ]
# where (a_0,a_1,a_2) are functions of (omega), reducing to (3/2,-2,1/2) at (omega=1). Using constant-step coefficients at variable (Delta t) is inconsistent.
#
# ## 6.6 Runge–Kutta
#
# An explicit (s)-stage Runge–Kutta method has
# [
# U^{(i)}=U^n+Delta tsum_{j<i}a_{ij}k_j,
# qquad
# k_i=M^{-1}R(U^{(i)}),
# ]
# [
# U^{n+1}=U^n+Delta tsum_i b_i k_i.
# ]
# Second-order consistency requires, among other conditions,
# [
# sum_i b_i=1,qquad
# sum_i b_ic_i=rac12.
# ]
# For third order additional rooted-tree/Butcher conditions are required. Thus a two-stage implementation is not automatically second order merely because it has two stages.
#
# ## 6.7 CFL, Fourier and acoustic scales
#
# For advection,
# [
# CFL=rac{|u|Delta t}{Delta x}.
# ]
# For diffusion,
# [
# Fo=rac{alphaDelta t}{Delta x^2}.
# ]
# For compressible/acoustic transport a characteristic scale is
# [
# CFL_a=rac{(|u|+a)Delta t}{Delta x},
# ]
# where (a) is the speed of sound. Bounds depend on dimension, stencil and time integrator.
#
# ## 6.8 Local time stepping and pseudo-time
#
# A pseudo-transient steady solve is
# [
# Mrac{U^{k+1}-U^k}{Delta	au}+R(U^{k+1})=0.
# ]
# (	au) is numerical time and must remain distinct from physical (t).
#
# ## 6.9 Adaptive stepping
#
# A basic controller based on an estimated order-(p) local error is
# [
# Delta t_{new}
# =
# Delta t_{old}
# left(rac{tol}{err}ight)^{1/(p+1)}
# ]
# with a declared safety factor and min/max bounds. Acceptance/rejection of a step is part of the algorithm state and therefore belongs in restart metadata.
#
# ## 6.10 Temporal/spatial error separation
#
# If
# [
# E(h,Delta t)approx C_s h^q+C_tDelta t^p,
# ]
# then refining (Delta t) while leaving (h) fixed eventually produces a spatial-error plateau. Conversely, spatial refinement with fixed (Delta t) eventually measures temporal error. A credible temporal-order study must demonstrate the appropriate asymptotic regime.
#
# ## 6.11 Restart semantics
#
# A one-step method needs (U^n). BDF2 additionally needs (U^{n-1}) and the previous time-step information. An adaptive integrator may also require controller state. Therefore a restart format that stores only the instantaneous field cannot, in general, reproduce a multistep trajectory bit-for-bit.
#
# ## 6.12 CFDX implementation
#
# Temporal code is in [temporal.h](../../../src/cfdx/physics/temporal.h), [low_storage_time_integration.h](../../../src/cfdx/physics/low_storage_time_integration.h) and [adaptive_cfl.h](../../../src/cfdx/physics/adaptive_cfl.h). The method registry also records the declared temporal methods and their verification status.
#
# The validation matrix is [test_temporal_order_matrix.cpp](../../../tests/validation/test_temporal_order_matrix.cpp); unit tests cover temporal history and variable-step coefficients.
#
# ## 6.13 Executable amplification check
# %%
import numpy as np
lam, dt = 2.0, 0.1
GE = 1-lam*dt
GI = 1/(1+lam*dt)
GCN = (1-0.5*lam*dt)/(1+0.5*lam*dt)
assert abs(GE) <= 1 and abs(GI) <= 1 and abs(GCN) <= 1
