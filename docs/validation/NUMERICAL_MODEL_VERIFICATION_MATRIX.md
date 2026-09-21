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

## Extended solver/discretisation verification

| ID | Domain | Verification | Level | Current gate |
|---|---|---|---|---|
| N023 | convection-diffusion | manufactured sinusoidal solution + 16/32/64/128 refinement | C | implemented |
| N024 | convection | upwind boundedness and refinement trend | C | implemented through N023 transport operator |
| N025 | transient integration | RK2 temporal refinement on u'=-u | C | implemented |
| N026 | transient integration | RK3 temporal refinement on u'=-u | C | algebraic; extend campaign |
| N027 | transient diffusion | PDE MMS with temporal refinement | C | pending solver driver |
| N028 | incompressible NS | Ghia Re=100/400 grid convergence | C | implemented |
| N029 | incompressible NS | manufactured divergence-free velocity/pressure solution | C | pending forcing interface |
| N030 | pressure-velocity coupling | SIMPLE zero-state + continuity/physical residual gates | C | implemented |
| N031 | turbulence | k-epsilon transport positivity/source response | C | implemented |
| N032 | turbulence | SST transport positivity/source response | C | implemented |
| N033 | turbulence | SA closure algebraic verification | A | implemented |
| N034 | turbulence | SA manufactured transport equation | C | pending SA transport solver |
| N035 | energy | 1-D conduction refinement | C | implemented |
| N036 | energy | transient energy manufactured solution | C | pending PDE campaign driver |
| N037 | CHT | two-region interface temperature/flux continuity | C | implemented |
| N038 | CHT | resistance-series analytical solution vs solver | C | pending coupled quantitative case |
| N039 | radiation | blackbody equilibrium | A/B | implemented |
| N040 | radiation | participating-media transport equilibrium | C | implemented |
| N041 | radiation | manufactured intensity/source solution | C | pending source-manufacture interface |
| N042 | compressible | Euler flux identical-state/free-stream preservation | A/B | implemented at flux level |
| N043 | compressible | 1-D shock tube against exact Riemann solution | C | pending compressible solver |
| N044 | compressible | isentropic vortex/free-stream preservation | C | pending compressible solver |
| N045 | low-Mach | asymptotic pressure/density consistency | A/B | implemented |
| N046 | transport properties | Sutherland/Prandtl/Schmidt regression | A | implemented |
| N047 | parallel | 1/2/4/8-rank numerical equivalence | C | pending distributed solver case |
| N048 | linear algebra | assembled vs matrix-free solution equivalence | C | pending common operator driver |
| N049 | mesh quality | orthogonal/skew/non-orthogonal convergence comparison | C | pending mesh campaign |
| N050 | reproducibility | deterministic repeated-run QoI and residual comparison | C | pending campaign harness |
| N051 | negative verification | invalid BC/NaN/negative property rejection | A/B | partly covered by unit tests |
| N052 | uncertainty | mesh/time uncertainty estimate from refinement | C | pending report aggregation |
| N053 | experimental/DNS | independent reference datasets with stated uncertainty | D | case-by-case |
| N054 | Fluent VMFL | VMFL reference/oracle matrix | E | matrix implemented; solver coverage varies |

### Required interpretation

An `implemented` gate means the executable gate exists; it does not mean that CFDX has passed the case. A PASS is only emitted after the executable, quantitative error criterion, convergence criterion, and reference comparison succeed. `pending` marks capabilities for which the current solver API does not yet expose the forcing or coupled driver required to perform a legitimate verification; these are not converted into artificial PASS results.

For quantitative campaigns the report records the mesh/time sequence, QoI, reference value, absolute error, relative error, residual/conservation measures, and observed order whenever the data permit it.
