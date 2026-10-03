# %% [markdown]
# 10 — Turbulence
#
# Turbulence modelling is split into exact averaged/filter equations, closure assumptions, wall treatment and numerical transport. CFDX currently exposes laminar, \(k\)-\(\varepsilon\), RNG \(k\)-\(\varepsilon\), realizable \(k\)-\(\varepsilon\), \(k\)-\(\omega\), SST, Spalart–Allmaras, Smagorinsky, WALE, dynamic one-equation, DES, DDES and IDDES model families. “Solver-ready” or “implemented” is not equivalent to qualification.
#
# ## 10.1 Reynolds decomposition
#
# \[
# u_i=\overline u_i+u_i',\qquad \overline{u_i'}=0.
# \]
# Averaging convection gives
# \[
# \overline{u_ju_i}
# =\overline u_j\,\overline u_i+\overline{u_j'u_i'}.
# \]
# Hence the RANS momentum equation contains
# \[
# -\rho\overline{u_i'u_j'}.
# \]
#
# ## 10.2 Boussinesq closure
#
# \[
# k=\frac12\overline{u_i'u_i'},
# \]
# \[
# -\rho\overline{u_i'u_j'}
# =
# 2\mu_tS_{ij}-\frac23\rho k\delta_{ij}.
# \]
# The eddy viscosity \(\mu_t\) is supplied by the turbulence model.
#
# ## 10.3 Generic two-equation transport
#
# \[
# \frac{\partial(\rho k)}{\partial t}
# +\nabla\cdot(\rho\mathbf uk)
# =
# P_k-\rho\beta^*k\omega
# +\nabla\cdot(D_k\nabla k)
# \]
# for a \(k\)-\(\omega\)-type model, while a \(k\)-\(\varepsilon\)-type form is
# \[
# \frac{Dk}{Dt}=P_k-\varepsilon+\nabla\cdot(D_k\nabla k),
# \]
# \[
# \frac{D\varepsilon}{Dt}
# =C_{\varepsilon1}\frac{\varepsilon}{k}P_k
# -C_{\varepsilon2}\frac{\varepsilon^2}{k}
# +\nabla\cdot(D_\varepsilon\nabla\varepsilon).
# \]
#
# ## 10.4 RNG \(k\)-\(\varepsilon\)
#
# RNG modifies the standard \(k\)-\(\varepsilon\) coefficients and adds an additional strain-dependent correction to the dissipation equation. The exact implemented production/correction functions must be taken from [turbulence_solver.h](../../../src/cfdx/physics/turbulence_solver.h); the label “RNG” is not a sufficient equation specification.
#
# ## 10.5 Realizable \(k\)-\(\varepsilon\)
#
# The model uses a variable \(C_\mu\):
# \[
# \mu_t=\rho C_\mu\frac{k^2}{\varepsilon},
# \]
# where \(C_\mu\) depends on local strain/rotation invariants. The realizable closure therefore differs mathematically from a constant-\(C_\mu\) model.
#
# ## 10.6 SST \(k\)-\(\omega\)
#
# A generic SST form is
# \[
# \frac{Dk}{Dt}=P_k-\beta^*k\omega+\nabla\cdot(D_k\nabla k),
# \]
# \[
# \frac{D\omega}{Dt}
# =\gamma\frac{\omega}{k}P_k-\beta\omega^2
# +\nabla\cdot(D_\omega\nabla\omega)+CD_{k\omega}.
# \]
# Blending:
# \[
# \Phi=F_1\Phi_1+(1-F_1)\Phi_2.
# \]
# Eddy viscosity is commonly limited in the form
# \[
# \nu_t=\frac{a_1k}{\max(a_1\omega,SF_2)}.
# \]
# Exact \(F_1,F_2\), cross-diffusion and production limiter formulas are implementation-specific and are covered by the SST tests.
#
# ## 10.7 Spalart–Allmaras
#
# The transported variable is \(\tilde\nu\):
# \[
# \frac{D\tilde\nu}{Dt}
# =
# c_{b1}(1-f_{t2})\tilde S\tilde\nu
# -\left(c_{w1}f_w-\frac{c_{b1}}{\kappa^2}f_{t2}\right)
# \left(\frac{\tilde\nu}{d}\right)^2
# +\frac1\sigma
# \left[
# \nabla\cdot((\nu+\tilde\nu)\nabla\tilde\nu)
# +c_{b2}|\nabla\tilde\nu|^2
# \right].
# \]
# The implementation also evaluates
# \[
# \chi=\frac{\tilde\nu}{\nu},\qquad
# f_{v1}=\frac{\chi^3}{\chi^3+c_{v1}^3}.
# \]
#
# ## 10.8 Wall distance
#
# The exact distance satisfies the Eikonal equation
# \[
# |\nabla d|=1,\qquad d=0\quad\text{on the wall}.
# \]
# CFDX contains exact/search-based, mesh-wave, Poisson, Eikonal, Hamilton–Jacobi, advection–diffusion and hybrid methods. The chosen distance affects SA/SST/hybrid closures and must therefore be treated as a numerical dependency.
#
# ## 10.9 LES
#
# Spatial filtering:
# \[
# \bar u_i=G*u_i,
# \]
# with subgrid stress
# \[
# \tau_{ij}^{sgs}
# =\overline{u_iu_j}-\bar u_i\bar u_j.
# \]
# Smagorinsky:
# \[
# \nu_t=(C_s\Delta)^2|S|.
# \]
# WALE replaces \(|S|\) by a function of the traceless square of the velocity-gradient tensor. Dynamic models determine coefficients from resolved-scale information rather than keeping a fixed \(C_s\).
#
# ## 10.10 DES/DDES/IDDES
#
# Hybrid methods modify the turbulence length scale with grid size and wall distance. A generic DES concept is
# \[
# d_{DES}=\min(d,C_{DES}\Delta).
# \]
# DDES/IDDES add shielding and stress/blending functions. Grid sensitivity is therefore part of their verification/validation contract.
#
# ## 10.11 Admissibility
#
# The solver must preserve declared physical bounds:
# \[
# k\ge k_{min}>0,\qquad
# \varepsilon\ge\varepsilon_{min}>0,\qquad
# \omega\ge\omega_{min}>0,\qquad
# \mu_t\ge0.
# \]
# Floors are admissibility protections; they must not be used to hide a divergent or incorrectly assembled transport equation.
#
# ## 10.12 Verification and qualification
#
# Verify source terms independently, then transport, positivity, wall-distance coupling, model limits and mesh sensitivity. Validation uses canonical turbulent flows and independent reference data. Qualification must state which model, wall treatment, \(y^+\) regime and geometry population are covered.
#
# ## 10.13 CFDX implementation
#
# Model registry/closures: [turbulence_models.h](../../../src/cfdx/physics/turbulence_models.h). Transport/source evaluation: [turbulence_transport.h](../../../src/cfdx/physics/turbulence_transport.h) and [turbulence_solver.h](../../../src/cfdx/physics/turbulence_solver.h). SA: [spalart_allmaras.h](../../../src/cfdx/physics/spalart_allmaras.h). Wall distance: [wall_distance.h](../../../src/cfdx/physics/wall_distance.h). SST: [sst_solver.h](../../../src/cfdx/physics/sst_solver.h).
#
# Tests include [test_turbulence_qualification_equations.cpp](../../../tests/unit/test_turbulence_qualification_equations.cpp), [test_sst_blending.cpp](../../../tests/unit/test_sst_blending.cpp), [test_wall_distance_turbulence.cpp](../../../tests/unit/test_wall_distance_turbulence.cpp) and the Phase-10 qualification/hardening matrix.
#
# ## 10.14 Executable Reynolds decomposition check
# %%
import numpy as np
uprime=np.array([-1.,1.,0.,2.,-2.])
assert np.isclose(np.mean(uprime),0.0)


# %% [markdown]
# ## 10.x Equation-level qualification contract
#
# A turbulence model is not qualified merely because its transport equations compile or a solver can instantiate the model. For each catalogue entry, retain:
#
# 1. exact equations and assumptions;
# 2. model constants and source references;
# 3. dimensional consistency of every term;
# 4. positivity/admissibility constraints;
# 5. wall-distance and wall-treatment dependencies;
# 6. term-by-term regression tests;
# 7. MMS or manufactured checks where meaningful;
# 8. solver-level convergence evidence;
# 9. physical benchmark evidence;
# 10. mesh/resolution and (y^+) sensitivity where applicable.
#
# Issue #473 defines the common 13-model qualification population. Wall-distance verification is a separate dependency: it may block a model path that requires it, but a wall-distance result must never be counted as turbulence-model validation.
#
# The Theory chapter therefore describes the equations and assumptions; the actual maturity decision belongs to the V&V/qualification matrix.
