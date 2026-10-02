# %% [markdown]
# 00 — Fluid Mechanics Foundations
#
# This chapter is the fluid-mechanics course of CFDX. The objective is to derive the equations before discussing their numerical approximation.
#
# ## 0.1 What CFD solves
#
# A CFD calculation solves a mathematical model defined by
# \[
# \oxed{geometry+fluid\ model+conservation\ laws+constitutive\ laws+initial/boundary\ conditions}
# \]
# The computed field U_h is an approximation of the continuous solution U. Spatial discretisation introduces h, temporal discretisation introduces Δt, and the algebraic solver introduces an iterative error.
#
# ## 0.2 Continuum hypothesis
#
# Macroscopic fields are
# \[
# \rho(\mathbf x,t),\quad p(\mathbf x,t),\quad T(\mathbf x,t),\quad \mathbf u(\mathbf x,t),\quad Y_k(\mathbf x,t).
# \]
# The Knudsen number
# \[
# Kn=\frac{\lambda_{mfp}}{L}
# \]
# measures the separation between molecular and continuum scales. Classical Navier–Stokes CFD assumes the continuum regime.
#
# ## 0.3 Eulerian and Lagrangian descriptions
#
# Eulerian fields are observed at fixed x; a Lagrangian description follows a material particle. For any scalar:
# \[
# \frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi.
# \]
# The first term is local change; the second is convection through the spatial gradient.
#
# ## 0.4 Kinematics
#
# \[
# \nabla\mathbf u=\mathbf D+\mathbf W,
# \quad
# \mathbf D=\frac12(\nabla\mathbf u+\nabla\mathbf u^T),
# \quad
# \mathbf W=\frac12(\nabla\mathbf u-\nabla\mathbf u^T).
# \]
# D measures deformation and W rigid-body rotation. Vorticity is
# \[
# \boldsymbol\omega=\nabla\times\mathbf u.
# \]
#
# ## 0.5 Mass conservation
#
# \[
# \frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.
# \]
# Equivalently:
# \[
# \frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0.
# \]
# Constant-density incompressible flow therefore satisfies
# \[
# \boxed{\nabla\cdot\mathbf u=0}.
# \]
#
# ## 0.6 Momentum and Cauchy stress
#
# \[
# \rho\frac{D\mathbf u}{Dt}=\nabla\cdot\boldsymbol\sigma+\rho\mathbf f,
# \qquad
# \boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau.
# \]
# Hence
# \[
# \rho\frac{D\mathbf u}{Dt}=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
# \]
#
# ## 0.7 Newtonian constitutive law
#
# \[
# \boldsymbol\tau=2\mu\mathbf D+\lambda(\nabla\cdot\mathbf u)\mathbf I.
# \]
# For constant μ and incompressibility:
# \[
# \nabla\cdot\boldsymbol\tau=\mu\nabla^2\mathbf u.
# \]
# Therefore:
# \[
# \boxed{\rho(\partial_t\mathbf u+\mathbf u\cdot\nabla\mathbf u)
# =-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f}.
# \]
#
# ## 0.8 Physical meaning of Navier–Stokes
#
# \[
# \underbrace{\rho\partial_t\mathbf u}_{local\ inertia}
# +\underbrace{\rho(\mathbf u\cdot\nabla)\mathbf u}_{convective\ inertia}
# =\underbrace{-\nabla p}_{pressure}
# +\underbrace{\mu\nabla^2\mathbf u}_{viscous\ diffusion}
# +\underbrace{\rho\mathbf f}_{body\ force}.
# \]
# Numerical terms must not be confused with these physical contributions.
#
# ## 0.9 Reynolds number
#
# With scales U and L:
# \[
# Re=\frac{\rho UL}{\mu}=\frac{UL}{\nu}.
# \]
# It measures inertia relative to viscosity. Large Re does not, by itself, prove turbulence.
#
# ## 0.10 Pressure as a constraint
#
# In incompressible flow pressure enforces the divergence constraint. Taking the divergence of momentum together with ∇·u=0 generates a pressure equation. This is the origin of pressure–velocity coupling.
#
# ## 0.11 Vorticity
#
# For incompressible constant-viscosity flow:
# \[
# \partial_t\boldsymbol\omega+(\mathbf u\cdot\nabla)\boldsymbol\omega
# =(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega+\nabla\times\mathbf f.
# \]
# The stretching term is a specifically three-dimensional mechanism.
#
# ## 0.12 Energy
#
# \[
# E=e+\frac12|\mathbf u|^2,
# \]
# \[
# \partial_t(\rho E)+\nabla\cdot[(\rho E+p)\mathbf u]
# =\nabla\cdot(\boldsymbol\tau\mathbf u-\mathbf q)+\rho\mathbf f\cdot\mathbf u+\dot q_v,
# ]
# with Fourier law
# \[
# \mathbf q=-k\nabla T.
# \]
#
# ## 0.13 Compressibility and equation of state
#
# An ideal-gas closure is
# \[
# p=\rho RT.
# \]
# The speed of sound is
# \[
# a=\sqrt{\left(\frac{\partial p}{\partial\rho}\right)_s},
# qquad Ma=\frac{U}{a}.
# \]
# Compressible solvers retain acoustic dynamics; low-Mach/incompressible formulations remove or constrain them.
#
# ## 0.14 Boundary and initial conditions
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
# A transient problem additionally requires initial data.
#
# ## 0.15 Boundary layers
#
# No-slip walls impose u=u_w. At large Re, viscous effects can be concentrated in thin wall-normal layers. This explains the importance of wall resolution, wall functions and wall-distance algorithms.
#
# ## 0.16 Bernoulli as a special case
#
# For steady inviscid flow with conservative body force:
# \[
# \frac p\rho+\frac12|\mathbf u|^2+gz=C.
# \]
# Bernoulli is a reduced special case, not the general Navier–Stokes model.
#
# ## 0.17 Nondimensional Navier–Stokes
#
# With x=Lx*, u=Uu*, t=(L/U)t*, p=ρU²p*:
# \[
# \partial_{t^*}\mathbf u^*+\mathbf u^*\cdot\nabla^*\mathbf u^*
# =-\nabla^*p^*+\frac1{Re}\nabla^{*2}\mathbf u^*+\mathbf f^*.
# \]
# This exposes dominant balances and supports similarity analysis.
#
# ## 0.18 PDE character
#
# Convective, diffusive and pressure constraints have different mathematical behaviour. Hyperbolic components propagate information along characteristics; elliptic constraints have global influence; parabolic terms diffuse information. This classification influences boundary conditions, fluxes and solver algorithms.
#
# ## 0.19 From fluid mechanics to CFD
#
# \[
# \boxed{physical\ model\rightarrow PDE\rightarrow integral\ balance
# \rightarrow mesh\rightarrow face\ fluxes\rightarrow algebraic\ system
# \rightarrow iterative\ solution}
# \]
# This chain is the organising principle for all later CFDX chapters.
#
# ## 0.20 CFDX traceability
#
# Field abstractions are under src/cfdx/core/field/. Mesh abstractions are under src/cfdx/core/mesh/. Boundary infrastructure is under src/cfdx/core/boundary/. Physical models are under src/cfdx/physics/.
#
# A source path establishes traceability, not qualification; implementation status must be established from executable evidence.
#
# %%
from __future__ import annotations
import numpy as np
rho,U,L,mu=1000.0,2.0,0.1,1.0e-3
Re=rho*U*L/mu
assert np.isclose(Re,2.0e5)
