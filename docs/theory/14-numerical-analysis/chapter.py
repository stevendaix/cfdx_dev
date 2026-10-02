# %% [markdown]
# 14 — Numerical Analysis
#
# ## 14.1 Consistency
#
# Let \(L\) be the continuous operator and \(L_h\) its discrete approximation. The truncation error is
# \[
# \tau_h=L_hu-Lu.
# \]
# Consistency means
# \[
# \|\tau_h\|\rightarrow0
# \quad\text{as}\quad h\rightarrow0.
# \]
# It is a property of the operator and exact solution, not of a particular computed result.
#
# ## 14.2 Taylor truncation
#
# For a smooth function:
# \[
# u(x+h)=u(x)+hu'(x)+\frac{h^2}{2}u''(x)+O(h^3).
# \]
# Therefore
# \[
# \frac{u(x+h)-u(x)}h
# =u'(x)+\frac h2u''(x)+O(h^2),
# \]
# a first-order forward difference.
#
# ## 14.3 Stability
#
# For
# \[
# \frac{dU}{dt}=AU,
# \]
# a one-step method produces
# \[
# U^{n+1}=GU^n.
# \]
# Stability over repeated steps requires bounded powers \(G^n\) in the declared regime. For scalar \(u'=\lambda u\), stability is expressed through the amplification factor \(G(z)\), \(z=\lambda\Delta t\).
#
# ## 14.4 Convergence
#
# Consistency and stability together provide convergence under the assumptions of the relevant numerical theorem. In practical CFD, nonlinear iteration error, boundary approximation, geometry error and solver tolerances must also be controlled.
#
# ## 14.5 Accuracy and observed order
#
# If
# \[
# E_h\approx Ch^p,
# \]
# then for refinement ratio \(r=h/H\):
# \[
# p_{obs}=
# \frac{\ln(E_H/E_h)}{\ln(H/h)}.
# \]
# A three-level sequence should be preferred for detecting asymptotic behaviour; a single two-grid slope can be contaminated by pre-asymptotic error.
#
# ## 14.6 Richardson extrapolation
#
# For a scalar quantity
# \[
# Q_h=Q+C h^p+O(h^{p+1}),
# \]
# two levels give
# \[
# Q_{ext}\approx Q_h+
# \frac{Q_h-Q_H}{r^p-1}.
# \]
# This assumes the solution is in the asymptotic regime and \(p\) is known or reliably estimated.
#
# ## 14.7 GCI
#
# A common grid-convergence index is
# \[
# GCI_{fine}=
# F_s\frac{|Q_2-Q_1|}{|Q_1|}
# \frac1{r^p-1},
# \]
# where \(F_s\) is a safety factor and the ordering of fine/coarse grids is declared. GCI is an uncertainty-style metric under its assumptions; it is not proof of model validity.
#
# ## 14.8 Dissipation and dispersion
#
# Numerical dissipation damps resolved amplitudes; numerical dispersion changes phase. For an advection wave \(u_t+cu_x=0\), the exact Fourier mode is
# \[
# u=e^{i(kx-\omega t)},\qquad\omega=ck.
# \]
# A discrete scheme has a numerical dispersion relation \(\omega_h(k)\), from which phase and amplitude errors can be analysed.
#
# ## 14.9 Monotonicity and boundedness
#
# Monotonicity concerns the discrete operator: a perturbation should not create a new local extremum under the declared conditions. A common scalar face bound is
# \[
# \min(\phi_P,\phi_N)\le\phi_f\le\max(\phi_P,\phi_N),
# \]
# but multidimensional monotonicity requires the complete stencil and coefficients.
#
# ## 14.10 Conditioning
#
# \[
# \kappa(A)=\|A\|\,\|A^{-1}\|.
# \]
# Perturbations satisfy a bound of the form
# \[
# \frac{\|\delta x\|}{\|x\|}
# \lesssim
# \kappa(A)\frac{\|\delta b\|}{\|b\|}
# \]
# under the standard assumptions. For nonsymmetric/non-normal systems, conditioning alone does not predict Krylov convergence.
#
# ## 14.11 Error budget
#
# A practical decomposition is
# \[
# E_{total}\lesssim
# E_{model}+E_{space}+E_{time}+E_{iteration}
# +E_{roundoff}+E_{input}.
# \]
# The terms are not generally statistically independent; the expression is a bookkeeping model, not a universal norm identity.
#
# ## 14.12 Verification methodology
#
# A reliable order study requires: exact manufactured/analytical reference; controlled mesh family; controlled time step; converged nonlinear/linear iterations; consistent norms; and documented geometry/BC errors. If one contribution dominates, the measured order belongs to that combined error regime rather than to the intended operator alone.
#
# ## 14.13 CFDX traceability
#
# Numerical operators live under [core/numerics](../../../src/cfdx/core/numerics/) and [core/linalg](../../../src/cfdx/core/linalg/). Representative validation includes [test_mesh_refinement_order.cpp](../../../tests/validation/test_mesh_refinement_order.cpp), [test_mms_scalar_diffusion.cpp](../../../tests/validation/test_mms_scalar_diffusion.cpp) and conservation/boundedness tests.
#
# ## 14.14 Executable observed-order check
# %%
import numpy as np
E=np.array([1e-2,2.5e-3,6.25e-4])
p=np.log(E[:-1]/E[1:])/np.log(2)
assert np.allclose(p,2.0)
