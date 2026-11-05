# CFDX GitHub CI architecture

The repository deliberately separates four CI responsibilities.

## 1. CFDX CI

`.github/workflows/cfdx-ci.yml`

Mandatory software build/regression layer:
- Release and Debug/UBSan builds;
- fast CTest regression;
- Python regression;
- diagnostics and artifacts before the final gate.

It does not contain N-level qualification logic.

## 2. Numerical maturity

`.github/workflows/cfdx-numerical.yml`

One workflow exposes independent jobs for N1-N16 from issue #461. The jobs are intentionally independent: a failure in N7 does not skip N8.

Each N-level job:
1. configures the same Release numerical build;
2. runs only its declared numerical test set;
3. prints diagnostics before the final gate;
4. uploads evidence.

N16 has no executable acceptance gate because advanced methods are explicitly deferred in the current maturity contract; its job validates the machine-readable audit instead of pretending that deferred work is qualified.

The test mapping is kept in this workflow so GitHub remains the orchestration layer while CTest remains the source of executable test definitions.

## 3. Total validation

`.github/workflows/cfdx-validation.yml`

Long/full qualification campaigns only:
- long CTest;
- wall-distance campaigns;
- validation report;
- thermal/radiation regression;
- application E2E;
- retained evidence.

It is intentionally not a second N1 gate.

## 4. Maturity records

`.github/workflows/cfdx-maturity.yml`

Validates the machine-readable #461 maturity and requirements audits on every PR/push.

## 5. Independent thermohydraulic reference-contract CIs

These are optional, case-scoped benchmark checks, not repository-wide required checks:

- `.github/workflows/th01-channel-dns.yml` — Moser–Kim–Mansour turbulent channel.
- `.github/workflows/th03-backward-facing-step.yml` — Driver–Seegmiller straight-wall backward-facing step.
- `.github/workflows/th07-mixed-convection.yml` — Wagner–Wetzel mixed-convection DNS.
- `.github/workflows/th08-t-junction.yml` — OECD/NEA–Vattenfall T-junction.

The `.github/workflows/th-physical-prerequisites.yml` workflow executes Ghia,
thermal, buoyancy and de Vahl Davis production preflight plus the comparator
mutation tests. It does not run TH-01/03/07/08 turbulent CFDX simulations and
must not be counted as their physical validation.

Each reference-contract workflow runs on changes to its own reference contract, its own workflow, or the shared reference-contract checker; it can also be dispatched manually. They validate reference metadata and derived dimensionless values only. Their artifacts explicitly state that CFDX physical validation remains BLOCKED and include a fail-closed solver-evidence status report. They must not be configured as mandatory branch-protection checks until solver-level campaigns exist.

## Required-check policy

Required checks should point only to workflows that always create a deterministic check:
- CFDX CI;
- the N-level checks that the project chooses to make blocking;
- CFDX numerical maturity audit.

Do not use path-filtered workflows for required checks. A required workflow must exist on every relevant PR; otherwise GitHub can report a missing/skipped check even when no numerical test failed.

## Removed duplication

The old standalone N1 workflow, the overlapping MPI/HDF5 build workflow, and the path-filtered numerical maturity workflow are replaced by the four-layer structure above.
