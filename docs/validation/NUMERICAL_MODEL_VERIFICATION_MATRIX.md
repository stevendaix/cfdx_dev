# Numerical model verification matrix

This matrix complements the Fluent VMFL campaign. It targets numerical methods and newly added CFDX models for which a small deterministic oracle is more useful than a large application case.

## Acceptance levels

- **A — algebraic/closed-form verification:** exact formula or invariant.
- **B — numerical solver verification:** known linear-system or ODE solution.
- **C — discretisation/model integration:** coupled numerical component with a physical invariant.

| ID | Model / numerical method | Reference | Level |
|---|---|---|---|
| N001 | Conjugate Gradient | exact SPD 3x3 solve | B |
| N002 | BiCGStab | exact non-symmetric 3x3 solve | B |
| N003 | GMRES | exact non-symmetric 3x3 solve | B |
| N004 | ILU(0) | finite/stable preconditioner application + residual reduction | B |
| N005 | Block-Jacobi | exact local block solve residual | B |
| N006 | Mixed-precision dot/norm | exact arithmetic reference | A |
| N007 | Fused reductions | exact dot/norm/max | A |
| N008 | TVD limiters | local boundedness invariant | A |
| N009 | Adaptive CFL | analytical controller equation | A |
| N010 | Pseudo-transient CFL | growth/saturation equation | A |
| N011 | Low-storage RK2/RK3 | exact exponential decay | B |
| N012 | Local time stepping | CFL/flux definition | A |
| N013 | Rusanov flux | identical-state physical flux consistency | A |
| N014 | Low-Mach preconditioning/EOS | ideal-gas density and scaling | A |
| N015 | Spalart–Allmaras closures | published closure equations | A |
| N016 | Boussinesq model | linear density/buoyancy law | A |
| N017 | Sutherland/Prandtl/Schmidt | transport correlations | A |
| N018 | Source linearisation | exact reconstruction of $S=S_u+S_p\phi$ | A |
| N019 | Blackbody/Rosseland/DOM | Stefan–Boltzmann and quadrature invariants | A |
| N020 | Radiation regime selector | optical-thickness thresholds | A |
| N021 | M1 pressure/turbulence closures | implemented closure equations | A |
| N022 | Thermal/CHT closures | resistance-network solution | A |

## Required extension for production-level confidence

These tests are deliberately small and deterministic. They do **not** replace multi-cell solver verification.

For every production numerical scheme we also require:

1. manufactured solutions;
2. at least three geometrically similar meshes;
3. $L_1$, $L_2$ and $L_\infty$ errors;
4. observed order of convergence;
5. temporal refinement for transient schemes;
6. conservation error independently from algebraic residuals;
7. serial/parallel equivalence when MPI is relevant;
8. negative/degenerate input tests;
9. regression of the reference result in the PDF report.

This follows established finite-volume verification practice: MMS plus systematic grid refinement is specifically used to determine observed order and expose formulation errors.

## Next solver-level waves

The highest-value next waves are:

- scalar Poisson MMS on orthogonal and non-orthogonal grids;
- convection-diffusion MMS for upwind/linear/TVD schemes;
- transient diffusion MMS for Euler/CN/BDF2;
- Navier–Stokes MMS for SIMPLE pressure/velocity coupling;
- turbulence MMS for k-epsilon, SST and Spalart–Allmaras transport;
- energy/CHT MMS;
- radiation equilibrium and participating-media transport;
- compressible Euler manufactured solution and shock-tube references;
- MPI partition invariance;
- matrix-free versus assembled operator equivalence.

The Ansys VMFL suite should remain the application-level verification layer. Ansys explicitly describes its VMFL cases as verification tests rather than complete model validation, so the two layers serve different purposes.
