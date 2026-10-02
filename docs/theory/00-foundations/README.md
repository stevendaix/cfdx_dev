# 00 — Foundations of CFDX Computational Fluid Dynamics

**Status: IN PROGRESS — first full scientific chapter.**

This chapter is the beginning of the CFDX Theory course. It explains the physics before the implementation and establishes the notation used by later chapters.

## 1. The CFD problem

CFD transforms a continuous physical problem into a finite set of discrete unknowns:

~~~text
Physical problem
  geometry + materials + IC + BC
              |
              v
Conservation laws + constitutive models
              |
              v
Continuous PDE system
              |
              v
Finite-volume integral equations
              |
              v
Face fluxes + source terms + time discretisation
              |
              v
Algebraic system A x = b
              |
              v
Linear / nonlinear / temporal iterations
              |
              v
Fields + residuals + conservation + QoIs
              |
              v
Verification -> validation -> qualification
~~~

Every arrow can introduce modelling, discretisation, implementation, iterative or round-off error. The Theory therefore explains every transformation rather than presenting equations as isolated formulas.

## 2. Continuum fields

The classical Navier–Stokes description uses:

- density rho;
- velocity u;
- pressure p;
- temperature T;
- total specific energy E.

Additional fields appear for turbulence, species, multiphase flow, radiation and other physics.

A CFD variable has a physical dimension, a mesh location, a governing equation, boundary conditions, a discretisation and a verification status.

## 3. Mass conservation

For a general fluid:

$$
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=S_\rho.
$$

For a closed single-component fluid:

$$
\boxed{\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.}
$$

The first term is local accumulation. The second is net outward mass flux. Thus:

$$
\text{accumulation}+\text{net outward flux}=0.
$$

This interpretation becomes the discrete finite-volume conservation equation.

## 4. Momentum conservation

The Cauchy momentum equation is

$$
\frac{\partial(\rho\mathbf u)}{\partial t}
+\nabla\cdot(\rho\mathbf u\otimes\mathbf u)
=
-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
$$

The terms are respectively momentum accumulation, convective transport, pressure force, viscous force and body force.

For a Newtonian fluid:

$$
\boldsymbol\tau =
2\mu\mathbf S+\lambda(\nabla\cdot\mathbf u)\mathbf I,
$$

with

$$
\mathbf S=\frac12(\nabla\mathbf u+\nabla\mathbf u^T).
$$

For incompressible constant-viscosity flow:

$$
\nabla\cdot\mathbf u=0,
$$

and

$$
\boxed{
\rho\left(
\frac{\partial\mathbf u}{\partial t}
+\mathbf u\cdot\nabla\mathbf u
\right)
=
-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f.
}
$$

The Theory must always distinguish conservative and advective forms because their discrete conservation properties are not identical.

## 5. Energy

A general total-energy equation is

$$
\frac{\partial(\rho E)}{\partial t}
+\nabla\cdot[(\rho E+p)\mathbf u]
=
\nabla\cdot(\boldsymbol\tau\cdot\mathbf u)
-\nabla\cdot\mathbf q+S_E.
$$

Fourier conduction is

$$
\mathbf q=-k\nabla T.
$$

The sign convention must be stated explicitly: q is the conductive heat-flux vector.

## 6. Closure

Conservation laws require constitutive relations.

Ideal-gas equation of state:

$$
p=\rho R T.
$$

Newtonian stress:

$$
\boldsymbol\tau=2\mu\mathbf S+\lambda(\nabla\cdot\mathbf u)\mathbf I.
$$

Fourier law:

$$
\mathbf q=-k\nabla T.
$$

The fundamental distinction is:

$$
\boxed{\text{conservation laws}+\text{constitutive models}=\text{closed PDE model}.}
$$

This distinction becomes essential for turbulence, radiation and multiphysics.

## 7. Incompressibility

For constant-density incompressible flow,

$$
\nabla\cdot\mathbf u=0.
$$

Pressure acts as a Lagrange multiplier enforcing this constraint. This is why pressure and velocity cannot be documented independently of the coupling algorithm.

## 8. Boundary conditions

A complete PDE problem contains:

$$
\text{PDE}+\text{domain}+\text{IC}+\text{BC}+\text{closure data}.
$$

Dirichlet:

$$
\phi=\phi_b.
$$

Neumann:

$$
\frac{\partial\phi}{\partial n}=g_b.
$$

Robin:

$$
a\phi+b\frac{\partial\phi}{\partial n}=c.
$$

CFDX separates boundary role, mathematical condition and field constraint. The Theory must first explain why these are different concepts, then explain their implementation.

## 9. Nondimensionalisation

With characteristic length L, velocity U and pressure rho U²:

$$
x=Lx^*,\quad u=Uu^*,\quad t=\frac{L}{U}t^*,\quad p=\rho U^2p^*.
$$

The incompressible momentum equation becomes

$$
\frac{\partial\mathbf u^*}{\partial t^*}
+\mathbf u^*\cdot\nabla^*\mathbf u^*
=
-\nabla^*p^*
+\frac1{Re}\nabla^{*2}\mathbf u^*,
$$

where

$$
\boxed{Re=\frac{\rho U L}{\mu}}.
$$

Later chapters will define Pe, Pr, Ma, Ra, Gr and Kn when they are actually used by a CFDX model or benchmark.

## 10. From PDE to finite volume

For

$$
\frac{\partial U}{\partial t}+\nabla\cdot\mathbf F=S,
$$

integrate over a control volume V:

$$
\frac{d}{dt}\int_V U\,dV+
\int_{\partial V}\mathbf F\cdot\mathbf n\,dS
=
\int_V S\,dV.
$$

Using Gauss' theorem:

$$
\boxed{
\frac{d}{dt}\int_V U\,dV+
\sum_f \mathbf F_f\cdot\mathbf S_f
=
\int_V S\,dV.
}
$$

The face area vector is

$$
\mathbf S_f=\mathbf n_f A_f.
$$

This is the central mathematical equation of the CFDX FVM chain.

## 11. Discrete conservation

For an internal face shared by cells P and N:

$$
\Phi_{P,f}=-\Phi_{N,f}.
$$

Therefore internal fluxes cancel when cell equations are summed.

This is why conservation verification must inspect face-flux antisymmetry independently of the final residual.

## 12. Numerical hierarchy

~~~text
Continuous physics
      |
      +-- modelling assumptions
      |
      v
Conservation/PDE form
      |
      +-- gradients
      +-- interpolation
      +-- convection
      +-- diffusion
      +-- source terms
      +-- boundary discretisation
      +-- time integration
      |
      v
Algebraic system
      |
      +-- linearisation
      +-- pressure/velocity coupling
      +-- linear solver
      +-- preconditioner
      |
      v
Numerical state
      |
      +-- residuals
      +-- conservation
      +-- boundedness
      +-- QoIs
~~~

A benchmark must state which layer it verifies. For example, Couette or Poiseuille diffusion tests do not automatically constitute a complete pressure-velocity Navier–Stokes validation.

## 13. Error taxonomy

CFDX Theory separates:

- modelling error;
- discretisation error;
- iterative error;
- round-off error;
- geometry/input error.

This prevents a numerical convergence result from being presented as physical validation.

## 14. Benchmark standard

Every benchmark page must contain:

1. physical problem;
2. governing equations;
3. assumptions;
4. exact/reference solution;
5. geometry and mesh family;
6. material properties;
7. IC/BC;
8. dimensionless parameters;
9. numerical schemes;
10. expected order;
11. error norms;
12. observed order;
13. conservation checks;
14. independent QoIs;
15. limitations;
16. bibliography;
17. possible improvements.

For two refinement levels:

$$
p_{obs}=\frac{\log(E_1/E_2)}{\log(h_1/h_2)}.
$$

An observed order is meaningless without specifying the mesh family and error metric.

## 15. CFDX implementation traceability

The current repository gives the following physical chain:

~~~text
core/mesh + core/geometry
          |
          v
core/field + core/boundary
          |
          v
core/numerics + core/fvm
          |
          +-- gradient
          +-- interpolation
          +-- divergence
          +-- Laplacian
          +-- convection
          +-- flux
          +-- temporal
          |
          v
physics/*
          |
          v
core/linalg/*
          |
          v
application / runtime / io
~~~

This is the starting point for the file-by-file audit. Each major equation will link to the exact implementation file, tests, benchmark and bibliography.

## 16. Required standard for every equation

Every important equation page must answer:

- Which physical principle produces it?
- Which assumptions are made?
- What does every symbol mean?
- What are the dimensions?
- What is the sign convention?
- Which equivalent forms exist?
- How is it discretised in CFDX?
- Which code implements it?
- What data does that code consume?
- What algebraic system is produced?
- Which benchmark verifies it?
- What is known to fail or degrade?
- What improvements are technically possible?
- Which references justify the formulation?

## 17. References

The Theory bibliography will use textbooks, standards, solver documentation and original benchmark papers. ASME V&V 20 provides the framework for quantitative validation and uncertainty-aware comparison; it is not itself a substitute for CFDX verification evidence.

Core references include Ferziger & Perić, Moukalled–Mangani–Darwish, Versteeg & Malalasekera, Roache, and the original benchmark papers cited by the V&V campaign.
