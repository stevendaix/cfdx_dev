# %% [markdown]
# Time Integration
#
# ## 6.1 Semi-discrete system
# \[
# M\frac{dU}{dt}=R(U,t).
# \]
# Temporal accuracy must be separated from spatial and iterative errors.
#
# ## 6.2 Explicit Euler
# \[
# U^{n+1}=U^n+\Delta t\,M^{-1}R(U^n).
# \]
# For \(u'=-\lambda u\), \(G=1-\lambda\Delta t\); stability of this scalar mode requires \(|G|\le1\).
#
# ## 6.3 Implicit Euler
# \[
# M\frac{U^{n+1}-U^n}{\Delta t}=R(U^{n+1}),
# \qquad
# G=\frac1{1+\lambda\Delta t}
# \]
# for the decay test equation.
#
# ## 6.4 Crank--Nicolson
# \[
# M\frac{U^{n+1}-U^n}{\Delta t}
# =\frac12[R(U^{n+1})+R(U^n)].
# \]
# It is formally second-order for smooth problems but may oscillate for stiff modes.
#
# ## 6.5 BDF2
# \[
# M\frac{3U^{n+1}-4U^n+U^{n-1}}{2\Delta t}=R(U^{n+1}).
# \]
# Variable-step coefficients must be derived from actual time nodes.
#
# ## 6.6 CFL and pseudo-time
# \[
# CFL=\frac{|u|\Delta t}{\Delta x}.
# \]
# A pseudo-transient method changes the path to steady state; it must not be confused with physical time integration.
#
# ## 6.7 Restart
# Restart must restore every history variable required by the selected temporal scheme.
#
# ## 6.8 CFDX traceability
# src/cfdx/physics/temporal.h  
# src/cfdx/physics/low_storage_time_integration.h  
# src/cfdx/physics/adaptive_cfl.h  
# src/cfdx/core/numerics/temporal.h  
# Tests: tests/unit/test_temporal.cpp, tests/validation/test_temporal_order_matrix.cpp
#
# %%
from __future__ import annotations
import numpy as np
lam,dt=2.0,0.1
G=1/(1+lam*dt)
assert 0<G<1
