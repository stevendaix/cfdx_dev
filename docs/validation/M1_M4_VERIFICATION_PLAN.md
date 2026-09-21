# CFDX M1-M4 verification plan

## Level A — exact analytical verification

Level A isolates discretization and implementation errors against closed-form solutions.

- M1 Couette
- M1 Poiseuille
- M3 1-D conduction
- M3 CHT resistance in series
- M4 two-surface black/gray radiation

Quantities include L2, Linf, conservation/oracle checks and observed spatial order.

## Level B — independent references

Level B removes the implementation from the reference definition.

- M2 turbulent-channel DNS reference (Moser-Kim-Mansour, Re_tau=178.12)
- M2 flat-plate ZPG correlations
- M2 turbulence closure identities
- M4 limiting-case radiation

The current Level-B executable verifies the independent reference data and closure/reference algebra. A resolved channel or flat-plate comparison is intentionally not represented as PASS until a multi-cell CFDX driver produces the required statistics.

## Level C — integration

Level C checks interactions among implemented M1-M4 components.

- turbulence transport + energy
- radiation + energy
- two-region CHT interface coupling

Acceptance is based on convergence, positivity and interface energy-flux balance.

## Campaign command

    python scripts/run_validation.py --build-dir build

The script executes all three levels and fails the campaign if any registered executable fails or is missing.

## Status semantics

PASS means the executable was run and its explicit gates passed. It does not mean that a higher-level physical validation has been completed when the required production case is not yet implemented.
