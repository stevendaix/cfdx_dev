# 00 — Foundations of CFDX

This chapter defines the mathematical language used by every later CFDX chapter. It starts from continuum mechanics and derives the equations that CFDX discretises.

## 1. Continuum mechanics

A continuum replaces a discrete collection of molecules by fields defined at points:
\[
\rho=\rho(\mathbf x,t),\quad p=p(\mathbf x,t),\quad T=T(\mathbf x,t),\quad \mathbf u=\mathbf u(\mathbf x,t).
\]
Here \(\rho\,[kg\,m^{-3}]\) is density, \(p\,[Pa]\) pressure, \(T\,[K]\) temperature and \(\mathbf u\,[m\,s^{-1}]\) velocity. The material derivative is
\[
\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi.
\]
It measures the rate seen by a fluid particle rather than at a fixed spatial location.

**CFDX path:** `src/cfdx/physics/`, `src/cfdx/core/field/`.

## 2. Conservation laws

For a control volume \(V\) with boundary \(\partial V\), outward normal \(\mathbf n\), the generic balance is
\[
\frac{d}{dt}\int_V q\,dV+\oint_{\partial V} \mathbf F_q\cdot\mathbf n\,dA
=\int_V s_q\,dV.
\]
This equation is the bridge from physics to finite volume: storage + net outward flux = volumetric source.

## 3. Continuity

The integral mass balance gives
\[
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.
\]
For constant density,
\[
\nabla\cdot\mathbf u=0.
\]
The incompressible condition is therefore a constraint on the velocity field, not an independent equation replacing momentum.

## 4. Momentum and Navier–Stokes

Cauchy's equation is
\[
\rho\frac{D\mathbf u}{Dt}=\nabla\cdot\boldsymbol\sigma+\rho\mathbf f,
\]
with
\[
\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau.
\]
For a Newtonian fluid,
\[
\boldsymbol\tau=2\mu\mathbf D+\lambda(\nabla\cdot\mathbf u)\mathbf I,
\quad
\mathbf D=\frac12(\nabla\mathbf u+\nabla\mathbf u^T).
\]
For incompressible constant-\(\mu\) flow this becomes
\[
\rho\left(\frac{\partial\mathbf u}{\partial t}+\mathbf u\cdot\nabla\mathbf u\right)
=-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f.
\]
Every term has dimensions \(N\,m^{-3}\): unsteady acceleration, convection, pressure force, viscous diffusion and body force.

**CFDX paths:** `src/cfdx/physics/incompressible.h`, `src/cfdx/physics/finite_volume_transport.h`.

## 5. Energy

With total specific energy \(E=e+|\mathbf u|^2/2\),
\[
\frac{\partial(\rho E)}{\partial t}
+\nabla\cdot[(\rho E+p)\mathbf u]
=\nabla\cdot(\boldsymbol\tau\mathbf u-\mathbf q)+\rho\mathbf f\cdot\mathbf u+\dot q_v,
\]
and Fourier conduction is
\[
\mathbf q=-k\nabla T.
\]
For a constant-property stationary conduction problem,
\[
\nabla\cdot(k\nabla T)+\dot q_v=0.
\]

**CFDX paths:** `src/cfdx/physics/thermal.h`, `src/cfdx/core/numerics/laplacian.h`.

## 6. Boundary and initial conditions

A well-posed transient problem requires initial data
\[
\mathbf u(\mathbf x,0)=\mathbf u_0(\mathbf x)
\]
and boundary conditions. Dirichlet specifies the field, \(\phi=\phi_b\); Neumann specifies its normal derivative,
\[
\nabla\phi\cdot\mathbf n=g_N;
\]
Robin combines both,
\[
a\phi+b\nabla\phi\cdot\mathbf n=c.
\]
CFDX represents these through boundary-condition and boundary-field abstractions rather than silently inferring a physical condition from a variable name.

**CFDX paths:** `src/cfdx/core/boundary/boundary_condition.h`, `src/cfdx/core/boundary/boundary_field.h`.

## 7. Dimensional analysis

The dimensions of the Navier–Stokes terms lead to
\[
Re=\frac{\rho UL}{\mu}=\frac{UL}{\nu},
\quad
Ma=\frac{U}{a},
\quad
Pe=\frac{UL}{\alpha},
\quad
Pr=\frac{\nu}{\alpha}.
\]
The Reynolds number compares inertia with viscosity; \(Pe\) compares advection with thermal diffusion. Nondimensionalisation is also a diagnostic: equations whose terms have inconsistent dimensions cannot be physically correct.

## 8. From PDE to CFDX

The chain is
\[
\text{physics}\rightarrow\text{PDE}\rightarrow\text{integral balance}
\rightarrow\text{mesh/control volumes}\rightarrow\text{face fluxes}
\rightarrow\text{sparse algebraic system}\rightarrow\text{linear/nonlinear solve}.
\]
The finite-volume discretisation is not an approximation applied after physics; it is the mechanism by which the integral conservation laws become computable.

## 9. Error taxonomy

For a numerical solution \(u_h\),
\[
e_h=u_h-u,
\qquad
\|e_h\|_2=\left(\sum_i w_i e_i^2\right)^{1/2}.
\]
Consistency concerns the truncation of the differential operator; stability concerns bounded amplification of perturbations; convergence means \(u_h\to u\) as \(h\to0\). Solver residual reduction is not, by itself, proof of discretisation convergence.

## 10. CFDX traceability

Each Theory equation must end with:
- exact implementation path;
- exact test path;
- benchmark/V&V evidence;
- assumptions and known limits.

**Core implementation:** `src/cfdx/physics/`, `src/cfdx/core/`, `src/cfdx/application/`.
**Core tests:** `tests/unit/`, `tests/validation/` where present.
**Reference:** Versteeg & Malalasekera, Ferziger & Perić, Moukalled et al.; see `docs/references/bibliography.bib`.

## 11. Verification target

The first algebraic checks are dimensional consistency, constant-field identities, conservation closure and manufactured-solution residuals. Physical validation is treated separately in chapter 15.

> **Status vocabulary:** Implemented, Verified, Validated and Qualified are distinct. This chapter documents theory and traceability; it does not promote any CFDX numerical method.
