# %% [markdown]
# 14 — Numerical Analysis
#
# ## 14.1 Consistency and truncation error
#
# Let $L$ be the continuous operator and $L_h$ its discrete approximation:
# $
# \tau_h=L_hu-Lu.
# $
# Consistency means
# $
# \|\tau_h\|\to0\qquad(h\to0).
# $
# A manufactured solution is especially useful because $Lu$ is known exactly and the discrete residual can be measured without conflating it with an unknown physical solution.
#
# ## 14.2 Stability and amplification
#
# For
# $
# \frac{dU}{dt}=AU,
# $
# a one-step scheme has
# $
# U^{n+1}=G(\Delta t A)U^n.
# $
# Stability requires bounded repeated application in the declared norm/regime. For scalar $u'=\lambda u$, define $z=\lambda\Delta t$ and amplification factor $G(z)$. Absolute stability is the set
# $
# \mathcal S=\{z:|G(z)|\le1\}.
# $
#
# ## 14.3 Consistency, stability and convergence
#
# Under the hypotheses of the relevant numerical theorem, consistency plus stability yields convergence. In a CFD implementation the observed error can also contain
# $
# E=E_{space}+E_{time}+E_{iteration}+E_{geometry}+E_{BC}+\cdots.
# $
# Therefore a converged nonlinear iteration is necessary but does not prove discretisation convergence.
#
# ## 14.4 Accuracy and observed order
#
# If
# $
# E_h=Ch^p+o(h^p),
# $
# then for $r=h_c/h_f$,
# $
# p_{obs}=\frac{\ln(E_c/E_f)}{\ln r}.
# $
# For a three-level study, the ratio
# $
# \frac{E_c-E_m}{E_m-E_f}
# $
# can be inspected for asymptotic behaviour before accepting a formal order claim.
#
# ## 14.5 Richardson extrapolation and GCI
#
# $
# Q_h=Q+Ch^p+\cdots,
# $
# gives
# $
# Q_{ext}\approx Q_h+
# \frac{Q_h-Q_H}{r^p-1}.
# $
# A common GCI expression is
# $
# GCI_{fine}=
# F_s\frac{|Q_2-Q_1|}{|Q_1|}
# \frac1{r^p-1}.
# $
# The fine/coarse ordering, refinement ratio and safety factor must be documented. GCI is an uncertainty-style estimate under its assumptions, not a model-validation result.
#
# ## 14.6 Modified-equation view
#
# A discrete operator can often be written formally as
# $
# L_hu=Lu+h^pL_{p+1}u+\cdots.
# $
# The leading extra term explains numerical diffusion, dispersion or anisotropy. This is particularly useful for comparing convection schemes beyond a single benchmark number.
#
# ## 14.7 Fourier dissipation and dispersion
#
# For
# $
# u_t+cu_x=0,
# \qquad
# u=e^{i(kx-\omega t)},
# $
# the exact relation is
# $
# \omega=ck.
# $
# A discrete method produces a numerical $\omega_h(k)$. Its real part controls phase error and its imaginary part controls numerical amplification/damping, subject to the Fourier convention used.
#
# ## 14.8 Monotonicity, boundedness and M-matrix structure
#
# For a scalar diffusion-like operator a desirable coefficient structure is often
# $
# a_P>0,\qquad a_{PN}\le0,
# $
# together with an appropriate diagonal dominance condition. This is not a universal requirement for every high-order or non-diffusive scheme. The property must be stated for the declared PDE, discretisation and sign convention.
#
# ## 14.9 Conditioning
#
# $
# \kappa(A)=\|A\|\|A^{-1}\|.
# $
# For a perturbed right-hand side,
# $
# \frac{|\delta x|}{|x|}
# \lesssim
# \kappa(A)\frac{|\delta b|}{|b|}.
# $
# For non-normal systems, transient behaviour can be much more subtle than this condition-number estimate.
#
# ## 14.10 Verification error budget
#
# A practical decomposition is
# $
# E_{total}\lesssim
# E_{model}+E_{space}+E_{time}+E_{iteration}
# +E_{roundoff}+E_{input}.
# $
# This is bookkeeping, not a norm identity. To measure a spatial order, the other dominant terms must be reduced below the spatial discretisation error.
#
# ## 14.11 CFDX qualification logic
#
# An operator qualification should document:
# 1. exact or manufactured reference;
# 2. controlled mesh family;
# 3. time-step refinement where relevant;
# 4. nonlinear and linear convergence;
# 5. conservation/boundedness invariants;
# 6. observed order and uncertainty;
# 7. limitations and failure regimes.
#
# A green regression test without a quantitative reference is a software regression check, not by itself a numerical qualification.
#
# ## 14.12 CFDX traceability
#
# Numerical operators live under [core/numerics](../../../src/cfdx/core/numerics/) and [core/linalg](../../../src/cfdx/core/linalg/). Representative validation includes [test_mesh_refinement_order.cpp](../../../tests/validation/test_mesh_refinement_order.cpp), [test_mms_scalar_diffusion.cpp](../../../tests/validation/test_mms_scalar_diffusion.cpp) and conservation/boundedness tests.
#
# ## 14.13 Executable observed-order check
# %%
import numpy as np
E=np.array([1e-2,2.5e-3,6.25e-4])
p=np.log(E[:-1]/E[1:])/np.log(2)
assert np.allclose(p,2.0)
