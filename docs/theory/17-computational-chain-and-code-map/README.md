> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 17 — CFDX Computational Chain and Code-to-Physics Map

**Status: CODE-AUDITED; this is the traceability backbone for the Theory chapters.**

## 1. Global chain

    case / mesh input
          ↓
    import + schema validation
          ↓
    mesh topology
          ↓
    geometry
          ↓
    fields + boundary conditions
          ↓
    numerical-method selection
          ↓
    FVM operators
          ↓
    physical equations
          ↓
    linearisation / coupling
          ↓
    sparse linear systems
          ↓
    Krylov / AMG / Schur / preconditioning
          ↓
    updated fields
          ↓
    residual + conservation + boundedness + QoI
          ↓
    V&V evidence

The chain is a traceability model, not a claim that every branch is production-qualified.

## 2. Mesh → geometry → operators

Mesh topology is stored in core/mesh. Geometry is derived in core/geometry.

The fundamental relationship is

$$
\text{topology}\rightarrow\text{geometry}\rightarrow
\text{control-volume integrals}.
$$

The geometry layer provides the quantities consumed by gradients, interpolation, diffusion and flux assembly.

## 3. Fields

Field metadata distinguishes location, dimension, unit, precision and storage. Current locations include cell, face, point and boundary.

A scalar cell field is

$$
\phi:\mathcal C\rightarrow\mathbb R.
$$

A vector field is

$$
\mathbf u:\mathcal C\rightarrow\mathbb R^3.
$$

A face value is not the same mathematical object as a cell value or cell gradient.

## 4. Boundary chain

    mesh patch
       ↓
    boundary role
       ↓
    field condition
       ↓
    mathematical constraint
       ↓
    FVM boundary contribution

Current implementation families are core/boundary and physics/boundary_constraint_fvm.h. The Theory must trace each supported condition to its discrete matrix/RHS contribution.

## 5. Numerical-method chain

The current numerical family contains separate contracts/dispatch for:
- gradient;
- stencil and boundary policy;
- interpolation;
- convection and assembly;
- divergence;
- flux;
- Laplacian;
- source terms;
- temporal integration;
- conservation;
- matrix-free operators;
- method registry and case selection.

This separation is the basis for N1: changing a gradient must not silently change face interpolation or temporal discretisation.

## 6. Diffusion

For

$$
\nabla\cdot(\Gamma\nabla\phi),
$$

the FV flux is

$$
F_f=\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
$$

The current Laplacian consumes geometry and reconstructed gradients. Therefore N2 and N3 cannot be qualified independently when the gradient is part of the non-orthogonal correction.

## 7. Convection

For

$$
\nabla\cdot(\rho\mathbf u\phi),
$$

$$
\dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f,
\qquad
F_f=\dot m_f\phi_f.
$$

The numerical method consists of mass flux, face reconstruction, limiter/boundedness policy and boundary treatment.

## 8. Physics

Current physics families include incompressible flow, pressure-velocity algorithms, scalar transport, thermal/energy, radiation, turbulence, wall distance, axisymmetric flow, low-Mach/compressible paths, multiphysics and conservation/boundedness controls.

The source tree is deliberately not collapsed into one generic solver equation: each model has its own assumptions and source terms.

## 9. Linear algebra

The algebraic system is

$$
A\mathbf x=\mathbf b.
$$

Current infrastructure contains Krylov solvers, sparse matrices, linear operators, block operators, AMG, MGR, Schur, field split, null-space, matrix-free and mixed-precision paths.

Algebraic convergence and physical convergence are separate evidence layers.

## 10. File-level audit rule

Each important source file must eventually have a row with:

| Item | Evidence |
|---|---|
| exact path | repository |
| equation/operator | Theory derivation |
| inputs/outputs | API and implementation |
| assumptions | code + numerical contract |
| sign convention | geometry/flux contract |
| tests | exact test target/file |
| benchmark | analytical/MMS/reference case |
| maturity | implemented / verified / qualified |
| limitations | negative evidence |
| improvement | technically justified |
| bibliography | reference identifiers |

The audit register in chapter 18 is the index; detailed equation pages remain in the domain chapters.

## 11. Promotion rule

A code path may be implemented without being verified. A verified operator may remain unvalidated. A validated case does not qualify every mesh family or solver backend.

This separation is mandatory throughout CFDX Theory and V&V.


## Scientific explanation standard

For every important equation, algorithm or data-model rule, document the motivation; definitions/units/sign conventions; derivation or formal rationale; discrete/FVM representation where applicable; actual CFDX execution path; example; error/limitation/sensitivity analysis; exact implementation files; executable verification; benchmark/V&V evidence; and stable bibliography/reference identifiers. Tables summarize explanations and do not replace them. Keep **Implemented / Verified / Validated / Qualified** distinct; documentation maturity never promotes numerical maturity.