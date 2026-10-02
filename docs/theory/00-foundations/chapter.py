# %% [markdown]
# 00 — Fluid Mechanics Foundations
#
# This is the fluid-mechanics course underlying CFDX. The numerical chapters should be read only after the continuous equations, constitutive assumptions, dimensions and boundary conditions are understood.
#
# ## 0.1 Continuum fields and assumptions
#
# The primary macroscopic fields are
# \[
# \rho(\mathbf x,t),\quad p(\mathbf x,t),\quad T(\mathbf x,t),\quad
# \mathbf u(\mathbf x,t),\quad Y_k(\mathbf x,t).
# \]
# The continuum approximation is associated with
# \[
# Kn=\frac{\lambda_{mfp}}{L}\ll1.
# \]
# CFDX classical Navier–Stokes models are intended for this continuum regime.
#
# ## 0.2 Material derivative
#
# Eulerian fields are functions of fixed spatial coordinates. Following a material particle gives
# \[
# \boxed{\frac{D\phi}{Dt}=
# \frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi}.
# \]
# This identity is the bridge between integral conservation laws and transport equations.
#
# ## 0.3 Kinematics
#
# \[
# \nabla\mathbf u=\mathbf D+\mathbf W,
# \]
# \[
# \mathbf D=\frac12(\nabla\mathbf u+\nabla\mathbf u^T),\qquad
# \mathbf W=\frac12(\nabla\mathbf u-\nabla\mathbf u^T).
# \]
# The rate-of-deformation tensor \(\mathbf D\) controls viscous strain in a Newtonian fluid. Vorticity is
# \[
# \boldsymbol\omega=\nabla\times\mathbf u.
# \]
#
# ## 0.4 Mass conservation
#
# Integral mass conservation:
# \[
# \frac{d}{dt}\int_V\rho\,dV+
# \oint_{\partial V}\rho\mathbf u\cdot\mathbf n\,dA=0.
# \]
# Local form:
# \[
# \boxed{\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0}.
# \]
# With constant \(\rho\):
# \[
# \boxed{\nabla\cdot\mathbf u=0}.
# \]
#
# ## 0.5 Momentum conservation
#
# Cauchy's equation is
# \[
# \rho\frac{D\mathbf u}{Dt}
# =\nabla\cdot\boldsymbol\sigma+\rho\mathbf f,
# \]
# with
# \[
# \boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau.
# \]
# Therefore
# \[
# \rho\frac{D\mathbf u}{Dt}
# =-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
# \]
#
# ## 0.6 Newtonian constitutive law
#
# \[
# \boldsymbol\tau
# =2\mu\mathbf D+\lambda(\nabla\cdot\mathbf u)\mathbf I.
# \]
# For constant \(\mu\) and incompressibility:
# \[
# \nabla\cdot\boldsymbol\tau=\mu\nabla^2\mathbf u.
# \]
# The incompressible Navier–Stokes equations become
# \[
# \boxed{
# \rho\left(
# \frac{\partial\mathbf u}{\partial t}
# +\mathbf u\cdot\nabla\mathbf u
# \right)
# =-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f
# }.
# \]
#
# ## 0.7 Physical interpretation
#
# \[
# \underbrace{\rho\partial_t\mathbf u}_{local\ inertia}
# +\underbrace{\rho\mathbf u\cdot\nabla\mathbf u}_{convection}
# =
# \underbrace{-\nabla p}_{pressure}
# +\underbrace{\mu\nabla^2\mathbf u}_{viscous\ diffusion}
# +\underbrace{\rho\mathbf f}_{body\ force}.
# \]
# The numerical operator implementing each term must preserve the physical sign convention.
#
# ## 0.8 Reynolds number and scaling
#
# \[
# Re=\frac{\rho UL}{\mu}=\frac{UL}{\nu}.
# \]
# Nondimensionalisation with \(x=Lx^*\), \(u=Uu^*\), \(t=(L/U)t^*\), \(p=\rho U^2p^*\) gives
# \[
# \frac{\partial\mathbf u^*}{\partial t^*}
# +\mathbf u^*\cdot\nabla^*\mathbf u^*
# =-\nabla^*p^*+\frac1{Re}\nabla^{*2}\mathbf u^*+\mathbf f^*.
# \]
#
# ## 0.9 Vorticity equation
#
# For constant-density, constant-viscosity flow:
# \[
# \frac{\partial\boldsymbol\omega}{\partial t}
# +(\mathbf u\cdot\nabla)\boldsymbol\omega
# =
# (\boldsymbol\omega\cdot\nabla)\mathbf u
# +\nu\nabla^2\boldsymbol\omega+\nabla\times\mathbf f.
# \]
# The stretching term vanishes in strictly two-dimensional flow.
#
# ## 0.10 Energy
#
# Total specific energy:
# \[
# E=e+\frac12|\mathbf u|^2.
# \]
# A conservative total-energy form is
# \[
# \frac{\partial(\rho E)}{\partial t}
# +\nabla\cdot[(\rho E+p)\mathbf u]
# =
# \nabla\cdot(\boldsymbol\tau\mathbf u-\mathbf q)
# +\rho\mathbf f\cdot\mathbf u+\dot q_v.
# \]
# Fourier conduction:
# \[
# \boxed{\mathbf q=-k\nabla T}.
# \]
#
# ## 0.11 Thermodynamic closure
#
# For an ideal gas:
# \[
# p=\rho RT,\qquad
# h=h(T),\qquad
# a^2=\left(\frac{\partial p}{\partial\rho}\right)_s.
# \]
# For constant \(c_p\), \(h=c_pT\), \(e=c_vT\), \(R=c_p-c_v\), and
# \[
# \gamma=\frac{c_p}{c_v},\qquad
# a=\sqrt{\gamma RT}.
# \]
# Mach number:
# \[
# Ma=\frac Ua.
# \]
# CFDX also contains incompressible, ideal-gas and Peng–Robinson EOS families; the numerical documentation must identify which closure a case uses.
#
# ## 0.12 Low-Mach versus incompressible versus compressible
#
# Incompressible flow imposes
# \[
# \nabla\cdot\mathbf u=0.
# \]
# Compressible flow solves density/pressure/energy consistently with an EOS and supports acoustic waves. Low-Mach methods retain thermal density effects while controlling the ill-conditioning associated with very small acoustic time scales. These are different mathematical models, not merely solver options.
#
# ## 0.13 Boundary and initial conditions
#
# Dirichlet:
# \[
# \phi=\phi_b.
# \]
# Neumann:
# \[
# \nabla\phi\cdot\mathbf n=g.
# \]
# Robin:
# \[
# a\phi+b\nabla\phi\cdot\mathbf n=c.
# \]
# Initial conditions specify the state at \(t=t_0\). A complete CFD problem is therefore PDE + constitutive closure + domain + BC + IC + material properties.
#
# ## 0.14 Boundary layers
#
# No-slip imposes
# \[
# \mathbf u=\mathbf u_w.
# \]
# At high \(Re\), a thin viscous layer forms near the wall. Wall distance \(y\), friction velocity \(u_\tau\) and
# \[
# y^+=\frac{u_\tau y}{\nu}
# \]
# become important for turbulence/wall-treatment modelling.
#
# ## 0.15 Bernoulli as a limiting case
#
# Steady inviscid flow with conservative body force satisfies
# \[
# \frac p\rho+\frac12|\mathbf u|^2+gz=C
# \]
# along a streamline. It is a special reduction, not a replacement for Navier–Stokes.
#
# ## 0.16 PDE character
#
# Convective terms transport information; viscous terms diffuse it; incompressibility introduces an elliptic constraint. This explains why CFD contains both local flux discretisation and globally coupled linear solves.
#
# ## 0.17 From physics to CFDX
#
# \[
# physical\ problem
# \rightarrow mathematical\ model
# \rightarrow PDE
# \rightarrow integral\ balance
# \rightarrow mesh
# \rightarrow reconstruction
# \rightarrow fluxes
# \rightarrow algebraic\ system
# \rightarrow nonlinear/linear\ solution
# \rightarrow V\&V.
# \]
#
# ## 0.18 Source traceability
#
# Physics lives mainly in [src/cfdx/physics](../../../src/cfdx/physics/), with thermodynamic and transport closures, while field/mesh/numerics provide the discretisation infrastructure. This chapter defines the equations; later chapters define their discrete implementation.
#
# ## 0.19 Executable Reynolds check
# %%
import numpy as np
rho,U,L,mu=1000.0,2.0,0.1,1.0e-3
assert np.isclose(rho*U*L/mu,2.0e5)
