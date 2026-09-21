# CFDX M1-M4 verification plan

## Level A — exact analytical verification

Level A isolates discretization and implementation errors against closed-form solutions.

- M1 Couette
- M1 Poiseuille
- M3 1-D conduction
- M3 CHT resistance in series
- M4 two-surface black/gray radiation

Quantities include L2, Linf, conservation/oracle checks and observed spatial order.

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
