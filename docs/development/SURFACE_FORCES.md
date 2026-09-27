# Surface forces and moments

CFDX uses a solver-independent surface force evaluator for aerodynamic/hydrodynamic forces.

## Convention

For a body surface, the supplied unit normal points **from the body into the fluid**:

\[
\mathbf{t}=(-p\mathbf{I}+\boldsymbol{\tau})\cdot\mathbf{n}
\]

and the fluid force on the body is

\[
\mathbf{F}=\int_S \mathbf{t}\,dS.
\]

Pressure and viscous contributions are retained separately.

The general C++ API is:

- `cfdx/physics/forces.h`
- `forces::integrate(samples, reference)`

It returns dimensional pressure/viscous/total force and moment, drag/lift/side projections, and pressure/viscous/total force coefficients.

## Reference quantities

For aerodynamic coefficients:

\[
q_\infty=\frac12\rho_\infty U_\infty^2,
\qquad
C_D=\frac{F_D}{q_\infty A_{ref}}.
\]

The reference area is supplied explicitly by the caller; CFDX does not infer it from the mesh. The reference length for moment coefficients is also supplied explicitly, so the force API does not assume a particular body scale or area/length relationship.

## Axisymmetric surfaces

The native axisymmetric evaluator consumes a meridional surface sample and applies the physical revolution measure exactly once:

\[
dS=2\pi r\,ds.
\]

For a no-swirl axisymmetric body:

\[
F_x=\int_\Gamma
\left[-p n_x+
\tau_{xx}n_x+
\tau_{xr}n_r\right]2\pi r\,ds.
\]

The axisymmetric evaluator is independent of the VMFL036 solver and can therefore be reused for bodies of revolution.

## Required diagnostics

Force extraction must reject:

- negative/non-finite surface measures;
- non-unit or non-finite normals;
- non-finite pressure;
- non-finite viscous traction/stress;
- invalid reference density, velocity or area.

A force result is diagnostic data, not a convergence criterion by itself. Validation must also report continuity, momentum residuals, pressure-correction convergence and mesh refinement.

This design follows the standard pressure/viscous surface-force decomposition used by established finite-volume CFD post-processing tools. citeturn0search0turn0search7
