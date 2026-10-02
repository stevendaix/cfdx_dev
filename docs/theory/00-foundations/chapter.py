# %% [markdown]
# Foundations

## 0.1 Continuum hypothesis

CFD represents macroscopic quantities by fields
\[
\rho=\rho(\mathbf x,t),\quad \mathbf u=\mathbf u(\mathbf x,t),\quad
p=p(\mathbf x,t),\quad T=T(\mathbf x,t).
\]
The continuum hypothesis replaces molecular-scale fluctuations by smooth fields at the scale resolved by the solver.

## 0.2 Kinematics

The material derivative is
\[
\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi.
\]
The velocity gradient is decomposed into deformation and rotation:
\[
\nabla\mathbf u=\mathbf D+\mathbf W,
\quad
\mathbf D=\frac12(\nabla\mathbf u+\nabla\mathbf u^T),
\quad
\mathbf W=\frac12(\nabla\mathbf u-\nabla\mathbf u^T).
\]

## 0.3 Conservation

For density q, flux F and source s,
\[
\frac{\partial q}{\partial t}+\nabla\cdot\mathbf F=s.
\]
For a fixed control volume,
\[
\frac{d}{dt}\int_Vq\,dV+\oint_{\partial V}\mathbf F\cdot\mathbf n\,dA
=\int_Vs\,dV.
\]
This integral statement is the mathematical foundation of finite volume discretisation.

## 0.4 Constitutive laws

Newtonian stress:
\[
\boldsymbol\sigma=-p\mathbf I+
2\mu\mathbf D+\lambda(\nabla\cdot\mathbf u)\mathbf I.
\]
Fourier heat flux:
\[
\mathbf q=-k\nabla T.
\]
These relations are constitutive assumptions and must not be confused with conservation laws.

## 0.5 Navier--Stokes

Mass:
\[
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.
\]
Momentum:
\[
\rho\frac{D\mathbf u}{Dt}
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
\]
For constant viscosity and incompressible flow:
\[
\nabla\cdot\mathbf u=0,\qquad
\rho\frac{D\mathbf u}{Dt}
=-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f.
\]

## 0.6 Boundary and initial conditions

Dirichlet:
\[
\phi=\phi_b.
\]
Neumann:
\[
\nabla\phi\cdot\mathbf n=g_N.
\]
Robin:
\[
a\phi+b\nabla\phi\cdot\mathbf n=c.
\]
Transient problems additionally require an initial condition.

## 0.7 Dimensionless groups

\[
Re=\frac{\rho UL}{\mu},\qquad
Ma=\frac{U}{a},\qquad
Pe=\frac{UL}{\alpha},\qquad
Pr=\frac{\nu}{\alpha}.
\]
Every term of a governing equation must have compatible dimensions. Nondimensionalisation also exposes dominant balances.

## 0.8 CFDX traceability

Field representation: src/cfdx/core/field/field.h  
Mesh representation: src/cfdx/core/mesh/mesh.h  
Boundary conditions: src/cfdx/core/boundary/  
Incompressible physics: src/cfdx/physics/incompressible.h  
Steady solver: src/cfdx/physics/steady_incompressible_solver.h

The source paths establish traceability, not numerical qualification.

## 0.9 Numerical chain

\[
\text{PDE}\rightarrow\text{integral balance}\rightarrow
\text{discrete operator}\rightarrow A x=b
\rightarrow\text{numerical solution}\rightarrow\text{evidence}.
\]

## 0.10 Evidence vocabulary

Implemented means code exists. Verified means a specified mathematical/software property has executable evidence. Validated means comparison with independent physical/reference evidence supports the declared case. Qualified means the declared population and acceptance gates are complete.

# %%
from __future__ import annotations
import numpy as np

rho,U,L,mu=1000.0,2.0,0.1,1e-3
Re=rho*U*L/mu
assert np.isclose(Re,2e5)
assert np.isfinite(Re)
