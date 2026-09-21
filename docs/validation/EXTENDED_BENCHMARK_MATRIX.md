# CFDX M1-M4 Extended Benchmark Matrix

This document expands the analytical campaign with a broad set of canonical CFD verification cases. The cases are deliberately split into **physics validation**, **numerical verification**, and **model/reference checks** so that a passing algebraic oracle is not confused with a solver-to-reference validation.

## Benchmark families

| ID | Family | Case | Level | Quantity |
|---|---|---|---|---|
| M1-01 | Incompressible | Couette | A | linear velocity profile |
| M1-02 | Incompressible | Plane Poiseuille | A | profile and flow rate |
| M1-03 | Incompressible | 2-D Taylor-Green vortex | A | exponential velocity/energy decay |
| M1-04 | Incompressible | Closed-domain continuity | A | integrated mass balance |
| M1-05 | Incompressible | Lid-driven cavity / Ghia | B | centreline velocity profiles |
| M2-01 | Turbulence | k-epsilon eddy viscosity | A | nu_t |
| M2-02 | Turbulence | SST limiting cases | A | nu_t |
| M2-03 | Turbulence | Smagorinsky scaling | A | delta^2 scaling |
| M2-04 | Turbulence | DES length-scale limiter | A | limiting length |
| M2-05 | Turbulence | wall-function viscous/log limits | A | U+ |
| M2-06 | Turbulence | laminar flat plate | B | Cf correlation |
| M2-07 | Turbulence | turbulent flat plate | B | Cf correlation |
| M2-08 | Turbulence | Moser-Kim-Mansour channel | B | y+/U+ reference normalization |
| M3-01 | Thermal | 1-D diffusivity | A | alpha |
| M3-02 | Thermal | two-layer slab | A | heat flux |
| M3-03 | CHT | interface conductance | A | q'' and flux continuity |
| M3-04 | Thermal | fully-developed heated duct | A | bulk-temperature rise |
| M3-05 | Thermal | transient diffusion oracle | A | Fourier-number decay |
| M3-06 | Thermal | source linearisation | A | exact linearisation at T0 |
| M4-01 | Radiation | Stefan-Boltzmann | A | Eb |
| M4-02 | Radiation | gray-surface emission | A | epsilon Eb |
| M4-03 | Radiation | black parallel surfaces | A | sigma(T1^4-T2^4) |
| M4-04 | Radiation | gray-surface exchange | A | resistance network |
| M4-05 | Radiation | view-factor closure/reciprocity | A | enclosure identities |
| M4-06 | Radiation | P1 equilibrium | A | q_rad = 0 |
| M4-07 | Radiation | DOM quadrature | A | sum(weights)=4 pi |
| M4-08 | Radiation | optically thin limit | A | zero absorption source |
| NUM-01 | Numerics | constant-field gradient | A | grad(phi)=0 |
| NUM-02 | Numerics | constant-field Laplacian | A | laplacian(phi)=0 |
| NUM-03 | Linear algebra | CG manufactured system | A | exact x |
| NUM-04 | Linear algebra | BiCGStab manufactured system | A | exact x |
| NUM-05 | Linear algebra | GMRES manufactured system | A | exact x |
| NUM-06 | Temporal | explicit Euler scalar ODE | A | one-step exact RHS |
| NUM-07 | Temporal | Crank-Nicolson scalar ODE | A | one-step amplification factor |
| NUM-08 | Thermodynamics | ideal-gas round trip | A | rho/p/T consistency |
| NUM-09 | Transport | Sutherland/Prandtl | A | reference values |
| NUM-10 | Momentum | zero-state fixed point | A | all momentum terms zero |
| NUM-11 | Turbulence | zero-strain production | A | Pk=0 |

## Reference benchmarks to add next

The following are intentionally listed as **future solver-level validation**, not as unit-test oracles:

1. **Ghia lid-driven cavity** at Re=100, 400, 1000, 3200, 5000, 7500 and 10000. Compare both centreline profiles, vortex-centre locations and extrema. Ghia et al. is the canonical reference family.
2. **Taylor-Green vortex 2-D and 3-D**. Compare instantaneous velocity, kinetic energy and dissipation. The 2-D case has an analytical solution; the 3-D case is a published high-fidelity benchmark.
3. **Blasius zero-pressure-gradient flat plate**. Compare boundary-layer thickness, displacement/momentum thickness and local skin friction against the laminar similarity solution.
4. **Turbulent channel** at Re_tau around 180, 395 and 590. Compare U+ profiles, wall shear, bulk velocity and Reynolds stresses against DNS.
5. **Backward-facing step**. Compare reattachment length and pressure recovery against experimental/reference data.
6. **Natural-convection heated cavity**. Compare Nu_avg, local Nu and centreline velocity/temperature against De Vahl Davis / Betts-Bokhari style references.
7. **Graetz thermal entrance**. Compare local Nusselt number and bulk temperature for laminar internal flow.
8. **Cylinder cross-flow** at low Reynolds number. Compare drag, lift and Strouhal number for Re=20/40/100.
9. **NACA 0012 / RAE 2822** once compressible/external-aerodynamics infrastructure is mature enough. Compare Cp, lift, drag and shock position where applicable.
10. **Compressible verification**: Sod shock tube, isentropic nozzle, oblique shock and manufactured Euler/Navier-Stokes solutions once the compressible solver is enabled.
11. **Mesh-manufactured solutions**: scalar Poisson, diffusion, advection-diffusion and momentum MMS on at least four grids, reporting L1/L2/Linf and observed order.
12. **Parallel reproducibility**: identical benchmark on 1, 2, 4 and 8 partitions, comparing field hashes, conservation residuals and QoIs.

## Acceptance policy

- A test must compute a CFDX quantity before it can be called a CFDX verification case.
- Analytical identities are Level A verification.
- Published DNS/experimental comparisons are Level B validation.
- Coupled turbulence/energy/radiation/CHT cases are Level C integration verification.
- A reference-data integrity test must never be reported as a solver-vs-reference PASS.
- Mesh refinement must be used for spatial-order claims.
- Temporal refinement must be used for time-order claims.
- Conservation error must be reported independently from algebraic solver residuals.
- No benchmark is considered scientifically validated merely because its setup compiles.

## Sources and benchmark rationale

The benchmark selection follows established CFD V&V practice: analytical solutions and manufactured solutions for verification, canonical numerical solutions/DNS for model validation, and experimental reference cases for application-level validation.

The Taylor-Green family is widely used for quantitative verification and validation of incompressible/high-fidelity CFD. NASA's CFL3D and Turbulence Modeling Resource provide established flat-plate, backward-facing-step, channel and airfoil cases. OpenFOAM's V&V catalogue likewise separates laminar, turbulent and heat-transfer reference cases.

The current PR therefore contains 38 executable analytical/component benchmarks, while reserving geometry-heavy solver benchmarks for dedicated drivers rather than creating misleading placeholder PASS tests. The existing Level-B test separately carries the Moser-Kim-Mansour DNS reference subset.
