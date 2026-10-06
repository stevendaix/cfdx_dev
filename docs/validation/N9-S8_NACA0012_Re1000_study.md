# N9-S8 — Detailed validation study: NACA 0012 at Re=1000

## Purpose

This case is the first external-flow qualification study for N9. The objective is to determine whether CFDX pressure--velocity algorithms produce converged, conservative and reproducible solutions for a canonical low-Reynolds-number airfoil without changing the physical problem between algorithms. The case is intentionally **laminar and incompressible**; it is not a turbulence-model validation case.

## Frozen physical problem

- Geometry: NACA 0012
- Chord: c = 1
- Freestream velocity: U_inf = 1
- Density: rho = 1
- Reynolds number: Re_c = 1000
- Kinematic viscosity: nu = 1e-3
- Angle of attack: alpha = 0 deg
- Steady first population
- Incompressible Navier--Stokes
- No-slip airfoil wall
- Uniform freestream at the external boundary

Re_c = U_inf c / nu = 1000.

Aerodynamic coefficients use dynamic-pressure normalization: Cd = D/(0.5 rho U_inf^2 c), Cl = L/(0.5 rho U_inf^2 c).

## Independent numerical oracle

The primary reference is Kurtulus (2015), International Journal of Micro Air Vehicles, 7(3), 301--326, DOI 10.1260/1756-8293.7.3.301. The study explicitly treats NACA 0012 at Re=1000 as laminar incompressible flow. At alpha=0 deg it reports approximately:

- Cd_ref = 0.12
- Cl_ref = 0 by symmetry

Kurtulus also reports the onset of characteristic unsteady vortex shedding at about alpha=8 deg. Thus alpha=0 deg is a deliberately conservative steady first population.

The Cd value is a **numerical literature reference**, not an exact analytical solution. CFDX must not relax tolerances merely to match it.

## Independent cross-checks

Kunz & Kroo (2000) provide low-Re laminar two-dimensional incompressible Navier--Stokes calculations for Re approximately 1000--12000 using C-grid topologies. Kouser et al. (2021) independently study NACA 0012 at Re=1000 using DNS and compare two- and three-dimensional aerodynamic behaviour. These references support the benchmark family and are retained as independent cross-checks.

## Mesh population

The qualification mesh family is the independent NASA/TMR NACA 0012 C-grid:

| Level | Grid | Purpose |
|---|---:|---|
| Coarse | 225 x 65 | debugging/refinement |
| Medium | 449 x 129 | primary result |
| Fine | 897 x 257 | refinement evidence |

The meshes are external validation assets, not CFDX-generated meshes. Materialization is controlled and offline; CI must never silently download them.

## Six-method N9 matrix

Every N9 pressure--velocity family must use exactly the same mesh, physics, boundary conditions, spatial discretization, stopping criteria and force-integration definition. No method-specific mesh, relaxed tolerance, reduced conservation requirement or hidden fallback is permitted.

For every run retain nonlinear history, pressure and momentum iterations, independently recomputed true residuals, continuity/mass imbalance, force balance, Cd/Cl/Cm, pressure and viscous drag, Cp(x/c), runtime and iteration count.

## Verification gates

1. **Residual convergence:** satisfy the frozen N9 convergence contract; retain true residual evidence.
2. **Conservation:** satisfy the declared continuity/mass-balance and force-balance gates.
3. **Symmetry:** at alpha=0 deg, Cl must converge toward zero without imposing Cl=0 in the solver.
4. **Mesh convergence:** coarse/medium/fine solutions must stabilize Cd, Cl, Cm, pressure drag, viscous drag and Cp(x/c).

Observed order is reported only where the refinement and error definition make an order estimate mathematically meaningful.

## Literature comparison

For each method and mesh report the literature discrepancy

E_Cd = |Cd_CFDX - 0.12| / 0.12.

This is a reported validation metric, not an arbitrary hard pass/fail threshold. For lift, report E_Cl = |Cl_CFDX| because symmetry supplies the independent reference value.

Where reference distributions are available, compare Cp(x/c) and Cf(x/c) using RMS and/or L-infinity error after documenting interpolation and surface-coordinate conventions. Do not reduce the validation to one integrated coefficient.

## Physical interpretation

At Re=1000 the flow is strongly viscous. An inviscid thin-airfoil result is therefore not an appropriate drag oracle. Likewise, NASA/TMR NACA 0012 results based on Spalart--Allmaras at high Reynolds number must not be mixed into this laminar Re=1000 validation case.

The oracle hierarchy is: (1) exact symmetry Cl -> 0; (2) convergence and conservation; (3) mesh refinement; (4) independent low-Re literature comparison, principally Cd ~= 0.12; (5) detailed Cp/Cf comparison where reference data are available.

## Qualification decision

N9-S8 remains **NOT QUALIFIED** until the controlled Total Validation campaign executes all three meshes with the complete six-method N9 matrix and retains the evidence. A close Cd match alone is never sufficient.

## References

Authoritative BibTeX entries are maintained in docs/references/bibliography.bib:

- kurtulus2015naca0012re1000
- kunzkroo2000ultralowre
- kouser2021naca0012re1000
- swansonlanger2016naca0012
- nasa_tmr_naca0012
