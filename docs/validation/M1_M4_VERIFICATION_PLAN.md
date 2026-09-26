# CFDX M1-M4 verification plan

## Level A — exact analytical verification

Level A isolates discretization and implementation errors against closed-form solutions.

- M1 Couette
- M1 Poiseuille
- M3 1-D conduction
- M3 CHT resistance in series
- M4 two-surface black/gray radiation

Quantities include L2, Linf, conservation/oracle checks and observed spatial order.

### Important M1 scope note

The current Couette and Poiseuille Level-A executables solve the corresponding 1-D finite-volume diffusion equation with analytical velocity references. They therefore verify the scalar diffusion/discretisation path used by the momentum assembly, but they are **not yet a full pressure-velocity-coupled Navier-Stokes validation**. The latter requires an actual incompressible solver driver with velocity boundary conditions, pressure correction, continuity monitoring and grid convergence.

The Poiseuille gate is deliberately diagnostic rather than fail-fast. `test_poiseuille_diagnostics` evaluates geometry/topology invariants, solver convergence, finite values, independent true residual, conservation, matrix structure, profile symmetry, velocity/flow-rate QoIs, analytical errors, and refinement order on N=8/16/32/64 before making a single final PASS/FAIL decision. Coarse-grid accuracy is reported as a diagnostic on N=8 and N=16; absolute L2/Linf, maximum-velocity, cell-centered flow-rate, and wall-flux accuracy gates apply on N=32 and N=64. The cell-centered flow-rate uses the midpoint quadrature `Q_h = 1/6 + h^2/3` for this stencil, so its coarse-grid error is expected and is not treated as a solver failure. An independent wall-flux QoI is also reported to separate field/quadrature effects from boundary-flux behavior.

## Extended component verification

The extended benchmark matrix (`test_benchmark_matrix`) contains 38 executable checks. These are **component-level verification/reference-oracle checks**, not solver-level CFD validation. They exercise implemented CFDX functions for turbulence closures, radiation limits, thermophysical models, linear algebra, temporal operators and numerical invariants. They are run by the campaign driver but must not be interpreted as 38 independent physical CFD cases.

A solver-level benchmark requires CFDX to advance a spatially discretised field and compare the resulting field or QoI with an analytical, DNS or experimental reference.

## Level B — independent references

Level B removes the implementation from the reference definition.

- M2 turbulent-channel DNS reference (Moser-Kim-Mansour, Re_tau≈180)
- M2 flat-plate ZPG correlations
- M2 turbulence closure identities
- M4 limiting-case radiation

The current Level-B executable verifies reference-data integrity and independent reference formulas. It does **not** compare a CFDX field to the DNS/engineering reference and therefore is not a solver-validation PASS. That comparison remains a Level-B solver benchmark to be implemented once the required multi-cell channel/flat-plate drivers exist.

## Level C — integration

Level C checks interactions among implemented M1-M4 components.

- turbulence transport + energy
- radiation + energy
- two-region CHT interface coupling

Acceptance is based on convergence, positivity and interface energy-flux balance.

## Campaign command

    python scripts/run_validation.py --build-dir build

The script executes the Level-A solver benchmarks, the extended component matrix, and Levels B/C. It fails the campaign if any registered executable is missing or returns non-zero. GitHub Actions builds exactly these validation executables and runs them with CTest.

## Status semantics

PASS means the executable was run and its explicit gates passed. It does not mean that a higher-level physical validation has been completed when the required production case is not yet implemented.


## V&V methodology basis

The campaign follows the usual separation between code verification, solution verification and validation. ASME V&V 20 emphasizes quantitative comparison at specified validation variables, while MMS provides a systematic way to generate exact solutions and assess observed order under grid refinement.

References:
- ASME V&V 20-2009 (R2021): https://www.asme.org/codes-standards/find-codes-standards/standard-for-verification-and-validation-in-computational-fluid-dynamics-and-heat-transfer
- Roache, P. J., "Code Verification by the Method of Manufactured Solutions", Journal of Fluids Engineering, 124(1), 4-10, 2002: https://doi.org/10.1115/1.1436090


## Validation taxonomy — scalar analytical vs coupled Navier–Stokes

The validation matrix now distinguishes two different verification classes; a green result in one class must not be reported as evidence for the other.

### A. Scalar analytical / discretisation verification

These tests solve an independent scalar finite-volume diffusion problem and compare it with a closed-form manufactured/analytical solution:

- `test_analytical_benchmarks`: Couette diffusion, Poiseuille diffusion, 1-D conduction, CHT resistance and radiation oracles.
- `test_poiseuille_diagnostics`: Poiseuille scalar FV campaign with geometry, matrix, conservation, QoI and refinement diagnostics.
- `test_mesh_refinement_order`: independent Poiseuille refinement/order check.
- `test_mms_scalar_diffusion`: manufactured scalar-diffusion verification.

These tests establish the correctness and convergence of the scalar FV/discretisation path. They **do not** by themselves validate pressure-velocity coupling, incompressibility, SIMPLE/PISO/PIMPLE/PIMPLE-like algorithms, or a coupled Navier–Stokes solution.

### B. Coupled incompressible Navier–Stokes verification

`test_phase9_acceptance` is classified separately as `coupled-ns`. It exercises the steady pressure-velocity system on the Couette problem with the exposed coupling algorithms and checks the velocity profile, pressure uniformity, transverse velocity, continuity and convergence metrics.

The Couette gate is deliberately based on the discrete cell-centred analytical profile, including `L2 < 1e-6`, `Linf < 1e-6`, and the exact discrete `Umax = 0.96875` for the N=16 cell-centred direction. Algorithm-invariance is additionally constrained to `max_abs_dU < 1e-5`.

The coupled classification is important for reporting: scalar Poiseuille/Couette diffusion results and coupled Couette results are complementary evidence, not interchangeable validation claims. The pressure-driven coupled Poiseuille campaign is tracked separately in PR #420 and is intentionally not duplicated here.
