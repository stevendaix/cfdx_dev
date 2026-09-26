# VMFL036 — Literature-aligned sphere case at Re=100

## Purpose

This is the clean CFDX reference configuration for the laminar flow past a sphere benchmark.

- sphere diameter: D = 1 m
- density: rho = 1 kg/m^3
- inlet velocity: U = 1 m/s
- dynamic viscosity: mu = 0.01 Pa.s
- Reynolds number: Re = rho U D / mu = 100
- reference projected area: A = pi D^2 / 4
- reference drag coefficient: Cd = 1.0895
- flow model: steady, incompressible, laminar
- geometry: 2-D axisymmetric sphere section
- outer circular domain radius: 50 D

The Cd=1.0895 datum is the literature anchor associated with Tabata & Itakura at Re=100 and cited by Ansys for VMFL036. The published paper covers sphere drag coefficients over Re=10–200.

## Important distinction from the Ansys VMFL036 page

The current Ansys VMFL036 page lists rho=1 kg/m^3, U=1 m/s, D=1 m and mu=0.02 Pa.s. Those displayed quantities imply Re=50, while the page reports Cd=1.0895, the Re=100 literature datum.

CFDX must not hide this discrepancy by widening a validation tolerance.

This file therefore defines the literature-aligned Re=100 case. It does not claim byte-for-byte identity with the proprietary Fluent project file VMFL036_FLUENT.cas.h5. Exact Fluent-project traceability remains tracked in Issue #440.

## Geometry

- sphere radius: R=0.5 m
- sphere surface: no-slip wall
- axis: axisymmetric boundary
- outer circular boundary: radius=50 m
- use a body-fitted mesh around the sphere
- resolve the separated wake; do not replace the axisymmetric case with an arbitrary Cartesian-box surrogate

## Quantity of interest

The primary quantity is total streamwise drag coefficient:

Cd = F_D / (0.5 rho U^2 A), with A = pi D^2 / 4.

Wall force must be obtained from the traction integral F = integral_S (-p n + tau.n) dS and decomposed into pressure drag, viscous drag and total drag.

## Acceptance status

This is a case definition/reference specification, not a solver PASS.

A solver-level PASS requires actual CFDX execution, convergence and conservation evidence, independent wall-force integration, coarse/medium/fine mesh refinement, quantitative comparison against Cd=1.0895, CTest/CI execution, and reproducible logs/report.

No numerical Cd tolerance is frozen by this reference-definition change.

## Calculation sheet

All dimensional inputs are SI.

| Quantity | Expression | Value |
|---|---|---:|
| Diameter | D | 1.000000 m |
| Radius | D/2 | 0.500000 m |
| Density | rho | 1.000000 kg/m^3 |
| Velocity | U | 1.000000 m/s |
| Dynamic viscosity | mu | 0.010000 Pa.s |
| Kinematic viscosity | nu=mu/rho | 0.010000 m^2/s |
| Reynolds number | rho U D / mu | **100.000000** |
| Dynamic pressure | 0.5 rho U^2 | 0.500000 Pa |
| Projected area | pi D^2 / 4 | 0.7853981634 m^2 |
| Reference Cd | literature datum | 1.089500 |
| Reference drag force | Cd (0.5 rho U^2 A) | **0.4278456495 N** |
| Outer radius | 50 D | 50.000000 m |

The reference force above is an oracle derived from the literature Cd, not a CFDX result. CFDX must independently integrate the sphere wall traction and obtain pressure and viscous contributions before comparison.

### Drag decomposition required from CFDX

For the sphere wall S:

F = integral_S (-p n + tau.n) dS

and, for the streamwise direction e_x:

F_D = e_x . F

Cd_pressure = F_D,pressure / (0.5 rho U^2 A)

Cd_viscous = F_D,viscous / (0.5 rho U^2 A)

Cd_total = Cd_pressure + Cd_viscous

The validation result must retain all three values. A total Cd matching the reference while one component is unavailable or silently omitted is not sufficient evidence.

### Re=100 consistency check

Re = (1.0 kg/m^3)(1.0 m/s)(1.0 m) / (0.010 Pa.s) = 100.

This is the intentional literature-aligned configuration. The documented Ansys value mu=0.02 Pa.s would instead produce Re=50 and therefore must not be mixed into this CFDX oracle.

## Physical CFD execution

The test `test_vmfl036_reference_case` now:
1. generates the body-fitted sphere mesh with `scripts/generate_vmfl036_mesh.py`;
2. imports the Gmsh mesh through the CFDX mesh I/O path;
3. solves steady incompressible Navier–Stokes at Re=100 with SIMPLE;
4. records convergence/continuity/momentum diagnostics;
5. integrates the sphere wall traction independently from the solver result;
6. reports pressure, viscous and total drag and the resulting (C_D).

The outer cylinder has radius (10D) and extends from (x=-10D) to (x=+10D). This is a finite-domain CFD calculation, not an analytical reconstruction of the reference force.

The CI comparison is intentionally diagnostic at this stage. The number printed as `Cd_reference=1.0895` is still the literature oracle; `Cd_total` is the CFDX result.

## Build / CI evidence

At the current PR head, both required GitHub Actions workflows have started and reached the compilation stage:

- CFDX CI, run #1903 — Configure: PASS; Build: running.
- CFDX MPI and parallel HDF5 integration, run #3531 — Configure: PASS; Build: running.

No compile PASS is claimed until the Build steps complete successfully. The CTest and validation-report stages are downstream of the build and are therefore not yet evidence for this case.

The CI workflow is deliberately configured to diagnose compilation failures before running CTest and to collect diagnostics before the final gate. This record should be updated with the final run conclusions and, if relevant, the exact failing job/log before the PR is considered ready.