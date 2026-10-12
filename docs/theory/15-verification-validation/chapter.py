# %% [markdown]
# 15 — Verification, Validation and Qualification
#
# ## 15.1 Vocabulary
#
# **Code verification:** is the implementation faithful to the specified mathematical model?
#
# **Solution verification:** is the computed numerical solution sufficiently converged with respect to discretisation and iteration?
#
# **Validation:** does the mathematical model reproduce the target physical behaviour within stated uncertainty?
#
# **Qualification:** has a declared capability population met explicit engineering acceptance criteria with retained evidence?
#
# These labels must never be collapsed into one generic “pass”.
#
# ## 15.2 Manufactured solutions
#
# Choose a smooth exact field $u_e$. If
# $
# Lu=f,
# $
# construct
# $
# f_{MMS}=L(u_e).
# $
# Solve the manufactured problem and measure
# $
# E_2=
# \left(
# \frac{\sum_PV_P|u_P-u_e(\mathbf x_P)|^2}
# {\sum_PV_P|u_e(\mathbf x_P)|^2}
# \right)^{1/2}.
# $
# MMS isolates discretisation implementation from uncertain physical data.
#
# ## 15.3 Spatial convergence
#
# If
# $
# E_h=Ch^p,
# $
# then for two levels:
# $
# p_{obs}=
# \frac{\ln(E_H/E_h)}{\ln(H/h)}.
# $
# Three or more levels are preferred because pre-asymptotic behaviour can otherwise mimic an incorrect order.
#
# ## 15.4 Temporal convergence
#
# Hold spatial error sufficiently below temporal error and use
# $
# E_{\Delta t}=C\Delta t^p.
# $
# Then repeat with at least two refinement ratios and keep nonlinear/linear tolerances sufficiently tight that iteration error does not dominate.
#
# ## 15.5 Iterative error
#
# Let $x^\*$ be the converged discrete solution and $x_k$ the current iterate. The iterative error is
# $
# e_{iter}=x_k-x^\*.
# $
# The residual
# $
# r_k=b-Ax_k
# $
# is related to this error by $r_k=-Ae_{iter}$. A small residual therefore implies small error only under conditioning assumptions.
#
# ## 15.6 Conservation verification
#
# Define
# $
# R_C=
# \frac{d}{dt}\sum_Pq_PV_P+
# \sum_{f\in\partial\Omega}F_f-
# \sum_Ps_PV_P.
# $
# Normalise with a declared scale $Q_{ref}$:
# $
# \epsilon_C=\frac{|R_C|}{Q_{ref}}.
# $
# Report both signed and absolute defects when cancellation matters.
#
# ## 15.7 Richardson/GCI
#
# For
# $
# Q_h=Q+C h^p+O(h^{p+1}),
# $
# Richardson extrapolation is
# $
# Q_{ext}=Q_h+\frac{Q_h-Q_H}{r^p-1}.
# $
# A common GCI form is
# $
# GCI_{fine}
# =F_s\frac{|Q_{fine}-Q_{coarse}|}{|Q_{fine}|}
# \frac1{r^p-1}.
# $
# It is meaningful only when the grids are in a defensible asymptotic regime.
#
# ## 15.8 Analytical benchmarks
#
# Typical CFDX analytical targets include Couette:
# $
# u(y)=U\frac yH,
# $
# and plane Poiseuille:
# $
# u(y)=\frac{G}{2\mu}y(H-y)
# $
# for the corresponding pressure-gradient convention.
#
# ## 15.9 Validation uncertainty
#
# A comparison should distinguish numerical uncertainty $U_{num}$, input uncertainty $U_{input}$, experimental uncertainty $U_{exp}$ and model-form discrepancy. A simple combined standard uncertainty may be
# $
# U_c=\sqrt{U_1^2+U_2^2+\cdots}
# $
# only when the stated independence/combination assumptions apply.
#
# ## 15.10 Evidence package
#
# Retain at minimum:
# case/setup revision, mesh and topology, physical constants, numerical schemes, solver tolerances, iteration history, output fields, derived quantities, software revision, hardware/runtime information when relevant, and post-processing scripts.
#
# ## 15.11 Qualification matrix
#
# A qualification statement must identify:
# population, geometry family, physics, numerical method, parameter range, acceptance criterion, evidence artifact, software revision and expiration/review rule.
#
# ## 15.12 CFDX evidence locations
#
# Executable checks live in [tests/unit](https://github.com/stevendaix/cfdx_dev/tree/master/tests/unit/), [tests/numerical](https://github.com/stevendaix/cfdx_dev/tree/master/tests/numerical/) and [tests/validation](https://github.com/stevendaix/cfdx_dev/tree/master/tests/validation/). Current evidence is retained under [docs/validation](https://github.com/stevendaix/cfdx_dev/tree/master/docs/validation/). Stable V&V publication/navigation is under [docs/vv](https://github.com/stevendaix/cfdx_dev/tree/master/docs/vv/).
#
# ## 15.13 Executable order check
# %%
import numpy as np
E=np.array([1e-2,2.5e-3,6.25e-4])
assert np.allclose(np.log(E[:-1]/E[1:])/np.log(2),2.0)


# %% [markdown]
# ## 15.13 From theory to executable V&V
#
# The V&V contract is a traceability chain:
#
# $
# \text{claim}\rightarrow\text{oracle}\rightarrow\text{metric}\rightarrow\text{criterion}\rightarrow\text{artifact}.
# $
#
# For each campaign, retain the software revision, setup, mesh, numerical schemes, solver criteria, raw result and post-processing definition. The evidence must identify whether it is code verification, solution verification, validation or qualification.
#
# ## 15.14 No residual-only validation
#
# A converged residual establishes an algebraic stopping condition, not physical correctness. Physical/reference comparison and discretisation uncertainty require independent evidence.
#
# ## 15.15 Qualification boundary
#
# Qualification is a bounded population statement. A passing case may support a capability claim only within the declared geometry, physics, numerical method and parameter range. Missing population members keep the claim incomplete.
