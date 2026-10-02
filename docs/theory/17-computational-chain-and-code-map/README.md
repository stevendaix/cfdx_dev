# 17 — CFDX Computational Chain and Code-to-Physics Map

**Status: IN PROGRESS — implementation traceability chapter.**

This chapter explains how a physical equation becomes a computed quantity in CFDX.

## 1. Global calculation chain

~~~text
CASE / MESH INPUT
      |
      v
IO + schema validation
      |
      v
Mesh topology
      |
      v
Geometry
      |
      v
Fields + boundary conditions
      |
      v
Numerical-method selection
      |
      v
FVM operators
      |
      +-- gradient
      +-- interpolation
      +-- divergence
      +-- diffusion / Laplacian
      +-- convection
      +-- flux
      +-- source
      +-- time integration
      |
      v
Physical equations
      |
      +-- continuity
      +-- momentum
      +-- energy
      +-- turbulence
      +-- radiation / multiphysics
      |
      v
Linearisation / coupling
      |
      v
Sparse linear system
      |
      v
Krylov / AMG / Schur / preconditioner
      |
      v
Updated fields
      |
      +-- residuals
      +-- conservation
      +-- boundedness
      +-- QoIs
      +-- output
      |
      v
V&V evidence
~~~

## 2. Mesh and geometry

The C++ core contains core/mesh and core/geometry families.

The mathematical distinction is:

$$
\text{topology}\rightarrow\text{geometry}\rightarrow\text{control-volume integrals}.
$$

Important quantities include cell volume, cell centroid, face centroid, face area and oriented face-area vector.

The Theory will give each quantity a definition, construction formula, orientation convention and mesh-quality consequence.

## 3. Fields

A scalar cell field is conceptually

$$
\phi:\mathcal C\rightarrow\mathbb R,
$$

and a vector field is

$$
\mathbf u:\mathcal C\rightarrow\mathbb R^3.
$$

Storage location matters: a cell value, face value and point value are different mathematical objects.

Relevant implementation families are core/field/field.h, core/field/storage.h and their source files.

## 4. Boundary chain

~~~text
mesh patch
    |
    v
boundary role
    |
    v
field condition
    |
    v
mathematical constraint
    |
    v
FVM boundary contribution
~~~

Relevant code families include core/boundary and physics/boundary_constraint_fvm.h.

The Theory will derive the actual discrete contribution of each supported condition.

## 5. Gradient chain

For scalar phi:

$$
\nabla\phi
$$

is a continuous differential quantity.

The discrete chain is:

~~~text
cell values
    |
    +-- neighbour stencil
    +-- face geometry
    +-- boundary information
    |
    v
gradient reconstruction
    |
    +-- Green-Gauss
    +-- least-squares
    +-- weighted least-squares
    +-- boundary policy
    +-- limiter where applicable
    |
    v
cell gradient
~~~

The Theory must distinguish gradient calculation from face-value reconstruction.

Relevant code includes core/numerics/gradient.h, gradient_stencil.h and core/fvm/least_squares_gradient.h.

## 6. Diffusion chain

For

$$
\nabla\cdot(\Gamma\nabla\phi),
$$

the finite-volume contribution is

$$
\sum_f \Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
$$

The numerical chain contains coefficient interpolation, face-normal gradient evaluation, non-orthogonal correction, limiting if present and matrix assembly.

The current implementation is centered around core/numerics/laplacian.h and its geometry/gradient dependencies.

## 7. Convection chain

For

$$
\nabla\cdot(\rho\mathbf u\phi),
$$

define

$$
\dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
$$

Then

$$
\Phi_{\phi,f}=\dot m_f\phi_f.
$$

Thus convection is not one algorithm: it combines mass flux, face-state reconstruction and boundary treatment.

The Theory will document each actual CFDX scheme and its applicable accuracy/boundedness evidence.

## 8. Incompressible momentum

The continuous equation is

$$
\rho\frac{\partial\mathbf u}{\partial t}
+\rho\nabla\cdot(\mathbf u\otimes\mathbf u)
=
-\nabla p+\nabla\cdot(2\mu\mathbf S)+\mathbf S_m.
$$

The computational chain is:

~~~text
U, p, materials
      |
      +-- grad(U)
      +-- convection
      +-- diffusion
      +-- pressure gradient
      +-- source terms
      |
      v
momentum equations
      |
      v
pressure-velocity coupling
      |
      v
velocity / pressure correction
      |
      v
continuity verification
~~~

Relevant current implementation families include physics/incompressible.h, physics/finite_volume_transport.h, physics/steady_incompressible_solver.h and physics/pressure_velocity*.h.

## 9. Pressure-velocity coupling

For incompressible flow:

$$
\nabla\cdot\mathbf u=0.
$$

Pressure enforces this constraint.

The Theory must therefore explain the momentum predictor, pressure correction, pressure Poisson/Schur structure, relaxation, null-space handling and continuity criterion.

Algorithm names such as SIMPLE, SIMPLEC, PISO or PIMPLE are not sufficient by themselves; the equations and update sequence must be documented.

## 10. Linear algebra

After discretisation:

$$
A\mathbf x=\mathbf b.
$$

The chain is:

~~~text
FVM assembly
    |
    v
Sparse A + RHS b
    |
    +-- CG
    +-- BiCGStab
    +-- GMRES
    +-- coupled/block methods
    |
    v
Preconditioner
    |
    +-- AMG
    +-- Schur
    +-- field-split
    |
    v
solution + true residual
~~~

The Theory must distinguish algebraic convergence from physical convergence and conservation.

## 11. Thermal and radiation

Thermal transport includes equations such as

$$
\rho c_p
\left(
\frac{\partial T}{\partial t}
+\mathbf u\cdot\nabla T
\right)
=
\nabla\cdot(k\nabla T)+S_T.
$$

Relevant current families include physics/thermal.h, energy_solver.h, cht_solver.h and radiation-related models/solvers.

Each model will receive its own derivation, assumptions and benchmark matrix.

## 12. Turbulence

RANS begins with

$$
\mathbf u=\overline{\mathbf u}+\mathbf u'.
$$

Averaging introduces additional stresses and therefore a closure problem.

~~~text
Navier-Stokes
      |
      v
Reynolds decomposition
      |
      v
averaged equations
      |
      v
closure problem
      |
      v
turbulence model
      |
      v
additional variables / effective stresses
~~~

Current code contains turbulence, SST, Spalart-Allmaras and wall-distance families. The Theory will explain their equations and distinguish implementation from qualification.

## 13. Time integration

For

$$
\frac{d\mathbf u}{dt}=\mathcal R(\mathbf u,t),
$$

backward Euler is

$$
\frac{\mathbf u^{n+1}-\mathbf u^n}{\Delta t}
=
\mathcal R(\mathbf u^{n+1},t^{n+1}).
$$

Every temporal method must be documented with consistency, order, stability, history and restart semantics.

## 14. Verification chain

~~~text
physical equation
       |
       v
discrete equation
       |
       v
implementation file
       |
       v
unit/operator test
       |
       v
MMS or analytical benchmark
       |
       v
mesh/time refinement
       |
       v
independent diagnostic
       |
       v
V&V evidence
~~~

A test that merely executes a code path is not automatically verification.

## 15. File-by-file audit standard

For each important numerical/physics source file, the Theory audit records:

| Item | Required information |
|---|---|
| File | exact repository path |
| Mathematical role | equation/operator/model |
| Inputs | fields, geometry, material data |
| Outputs | field, flux, matrix or diagnostic |
| Equation | governing expression |
| Assumptions | modelling/discretisation assumptions |
| Sign conventions | orientation and flux convention |
| Tests | exact test files |
| Benchmark | reference problem |
| Status | implemented / verified / qualified |
| Limitation | known negative evidence |
| Improvement | technically justified next step |
| References | bibliography entries |

This is the standard for the continuing source-tree audit.

## 16. Status vocabulary

The Theory deliberately uses three separate labels:

- Implemented: a code path exists and is exercised.
- Verified: mathematical properties have executable evidence.
- Qualified: the defined V&V campaign supports the intended scope.

These labels must never be collapsed into a single generic “supported” statement.
