# %% [markdown]
# 06 — Time Integration
#
# ## 6.1 Semi-discrete equation
#
# After spatial discretisation:
# \[
# M\frac{dU}{dt}=R(U,t).
# \]
# \(M\) contains the discrete storage measure; \(R\) is the assembled residual/flux balance. A temporal convergence test must keep spatial and nonlinear errors sufficiently small.
#
# ## 6.2 Explicit Euler
#
# \[
# U^{n+1}=U^n+\Delta t\,M^{-1}R(U^n).
# \]
# For \(u'=-\lambda u\), \(\lambda>0\):
# \[
# G_E=1-\lambda\Delta t,\qquad |G_E|\le1
# \]
# gives \(0\le\lambda\Delta t\le2\) for this scalar decay problem.
#
# ## 6.3 Implicit Euler
#
# \[
# \frac{U^{n+1}-U^n}{\Delta t}
# =-\lambda U^{n+1},
# \qquad
# G_I=\frac1{1+\lambda\Delta t}.
# \]
# It is first-order accurate but A-stable for the linear test equation.
#
# ## 6.4 Crank–Nicolson
#
# \[
# \frac{U^{n+1}-U^n}{\Delta t}
# =\frac12[R(U^{n+1})+R(U^n)].
# \]
# For \(u'=-\lambda u\),
# \[
# G_{CN}=\frac{1-\frac12\lambda\Delta t}
# {1+\frac12\lambda\Delta t}.
# \]
# The method is second-order but can preserve oscillatory behaviour for stiff transients.
#
# ## 6.5 BDF2
#
# \[
# \frac{3U^{n+1}-4U^n+U^{n-1}}{2\Delta t}
# =R(U^{n+1}).
# \]
# For variable time steps the coefficients must be derived from the actual \(t_{n-1},t_n,t_{n+1}\); copying constant-step coefficients is not mathematically valid.
#
# ## 6.6 CFL and diffusion restrictions
#
# For advection:
# \[
# CFL=\frac{|u|\Delta t}{\Delta x}.
# \]
# For a diffusion equation a characteristic explicit restriction is
# \[
# Fo=\frac{\alpha\Delta t}{\Delta x^2},
# \]
# with the allowable bound depending on dimension and spatial stencil. These are stability constraints, not universal accuracy criteria.
#
# ## 6.7 Pseudo-transient continuation
#
# A steady residual equation
# \[
# R(U)=0
# \]
# can be advanced using
# \[
# M\frac{U^{n+1}-U^n}{\Delta\tau}+R(U^{n+1})=0.
# \]
# \(\tau\) is numerical pseudo-time. It must not be reported as physical time.
#
# ## 6.8 Adaptive time step
#
# A controller may target a CFL or local truncation/error indicator:
# \[
# \Delta t_{new}
# =\Delta t_{old}
# \left(\frac{tol}{err}\right)^{1/(p+1)}
# \]
# for a basic order-\(p\) error model, with safety factors and declared min/max limits.
#
# ## 6.9 Restart invariants
#
# A restart must restore all history required by the scheme: for BDF2 at least \(U^n,U^{n-1}\), while one-step methods need only the current state plus controller state if adaptive stepping is used.
#
# ## 6.10 CFDX implementation
#
# Temporal code is in [temporal.h](../../../src/cfdx/physics/temporal.h), [low_storage_time_integration.h](../../../src/cfdx/physics/low_storage_time_integration.h) and [adaptive_cfl.h](../../../src/cfdx/physics/adaptive_cfl.h). The exact supported schemes and coefficients must be kept synchronized with this chapter.
#
# Verification: [test_temporal_order_matrix.cpp](../../../tests/validation/test_temporal_order_matrix.cpp) and the temporal unit tests.
#
# ## 6.11 Executable stability check
# %%
import numpy as np
lam, dt = 2.0, 0.1
GE = 1-lam*dt
GI = 1/(1+lam*dt)
GCN = (1-0.5*lam*dt)/(1+0.5*lam*dt)
assert abs(GE) <= 1 and abs(GI) <= 1 and abs(GCN) <= 1
