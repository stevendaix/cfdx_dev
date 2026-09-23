# Phase 9 — M1 incompressible acceptance

## Scope

This document is the acceptance contract for the M1 steady incompressible solver. It deliberately separates:

1. implementation of the pressure-velocity algorithm;
2. analytical/operator verification;
3. solver-level physical validation;
4. mesh convergence;
5. conservation;
6. regression/CI evidence.

The primary published reference for the lid-driven cavity is Ghia, Ghia & Shin, *Journal of Computational Physics* 48 (1982), 387–411, DOI 10.1016/0021-9991(82)90058-4.

## Acceptance matrix

| ID | Requirement | Evidence in this PR |
|---|---|---|
| 9.1 | Momentum | Integrated pressure-driven Navier–Stokes solve; final equation residual is checked independently |
| 9.2 | Continuity | Cell-wise L1/Linf and normalized continuity gates |
| 9.3 | SIMPLE | Analytical pressure-driven Poiseuille execution |
| 9.4 | SIMPLEC | Analytical pressure-driven Poiseuille execution |
| 9.5 | PISO | Multiple pressure correctors exercised on the analytical case |
| 9.6 | PIMPLE | PIMPLE path exercised with explicit outer/pressure-corrector controls |
| 9.7 | Rhie–Chow | Existing primitive and coupled-solver regressions plus integrated acceptance |
| 9.8 | Pressure gauge | Pure-Neumann non-zero initial pressure is reset to the declared reference |
| 9.9 | Convection schemes | Upwind and second-order-upwind paths execute through the coupled solver |
| 9.10 | Bounded deferred correction | Bounded and unbounded paths are executed and compared against the same analytical solution |
| 9.11 | Residual vs physical convergence | Acceptance requires linear residuals, final momentum-equation residual, continuity and field-change metrics |
| 9.12 | Mass conservation | Continuity L1/Linf/normalized gates are mandatory |
| 9.13 | Cavity/Ghia | Existing executable Ghia Re=100/400 campaign remains a required CI gate |
| 9.14 | Mesh convergence | Existing 32/64/128 Ghia Re=100 refinement gate remains a required CI gate |
| 9.15 | Canonical benchmarks | Pressure-driven Poiseuille is solved by the coupled NS solver against its exact quadratic profile |
| 9.16 | Independent reference comparison | Ghia published-reference comparison is the independent Level-B acceptance; OpenFOAM is treated as an additional implementation/reference comparison, not as the mathematical oracle |
| 9.17 | Regression suite | Phase-9 acceptance, steady-solver and Ghia executables are all part of the validation CI |

## Analytical pressure-driven Poiseuille case

The unit channel uses:

- L = 1;
- H = 1;
- rho = 1;
- nu = 0.1;
- p(inlet) = 1;
- p(outlet) = 0;
- no-slip at y=0 and y=1;
- zero-gradient velocity at inlet/outlet;
- zero-gradient pressure on walls.

The exact fully-developed solution is:

    u(y) = (-dp/dx)/(2 nu) * y * (1-y)

and therefore, for -dp/dx = 1 and nu = 0.1:

    u(y) = 5 y (1-y)

with u_max = 1.25.

The test evaluates the CFDX cell-centred field at cell centres and requires a bounded L2+Linf profile error, a converged coupled solve, and finite physical diagnostics.

## Pressure-velocity algorithms

The same physical problem is run through:

- SIMPLE;
- SIMPLEC;
- PISO with two pressure correctors;
- PIMPLE with two outer correctors and two pressure correctors.

The acceptance does not rank the algorithms. It verifies that each produces the same physical solution within the analytical tolerance and satisfies the same conservation/physical-residual gates.

## Pressure gauge

A separate pure-Neumann run starts with an arbitrary non-zero pressure field and no fixed-pressure boundary. The declared reference cell/value must remove the pressure null space without leaving the arbitrary initial offset.

## Convection and boundedness

The coupled solver executes both:

- first-order upwind;
- second-order upwind;

and both bounded and unbounded convection paths.

The acceptance is deliberately based on the analytical pressure-driven solution and physical convergence, not merely on successful construction of a scheme enum.

## Cavity / Ghia

test_ghia_cavity remains the independent reference gate. It executes:

- Re=100 on 32x32, 64x64 and 128x128 grids;
- Re=400 on 64x64;
- both centreline velocity profiles;
- RMS and maximum absolute errors;
- continuity;
- final momentum-equation residual;
- observed refinement order.

The reference is the published Ghia et al. dataset. Ghia et al. used high-resolution uniform grids for the driven-cavity benchmark, including a 129x129 solution at Re=1000.

## Scientific rule

A low linear residual is never sufficient for a Phase-9 PASS.

A run must also satisfy:

- finite solution fields;
- physical momentum-equation residual;
- continuity L1/Linf/normalized limits;
- field-change convergence;
- analytical/reference QoI error limits;
- mesh/refinement evidence where applicable.

## CI gate

.github/workflows/cfdx-validation.yml builds and executes:

- test_phase9_acceptance;
- test_steady_incompressible_solver;
- test_ghia_cavity;

in addition to the existing M1–M4 validation executables.

No Phase-9 roadmap item is promoted to green solely because a target exists or a code path compiles. The green state is tied to executable quantitative evidence.
