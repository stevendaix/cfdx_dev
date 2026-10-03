# CFDX User Guide

**Status: FOUNDATION ONLY — this domain is intentionally not finished.**

This directory currently provides only the navigation and minimum conceptual contract needed to prepare the future User Guide. It must not be interpreted as a complete user manual.

## Current scope

The present skeleton covers:

- getting started;
- basic concepts;
- the intended CFDX workflow;
- a reference/documentation boundary.

Detailed operational documentation is still to be written and verified against the actual application/API.

## Planned user-guide scope

The completed guide will progressively cover:

1. installation and supported environments;
2. case creation;
3. mesh import and inspection;
4. physics and material setup;
5. boundary and initial conditions;
6. numerical-method selection;
7. solver execution and convergence monitoring;
8. restart/checkpoint workflows;
9. post-processing and visualisation;
10. probes and time series;
11. parallel execution;
12. GPU execution;
13. CLI/Python/API usage;
14. troubleshooting;
15. complete verified examples.

No command, API signature or GUI workflow should be documented here until it has been checked against the corresponding implementation.

## Documentation boundary

- **Theory** explains the mathematics.
- **User** explains how to operate CFDX.
- **Developer** explains implementation contracts.
- **V&V** defines evidence and acceptance.
- **docs/validation/** stores executed evidence.

Existing user-facing documents remain migration sources until each capability is reconciled with the actual implementation and evidence.

## Maturity rule

A user-guide page may describe a capability as available only when its implementation status is known. Validation/qualification claims belong to V&V and validation evidence, not to this guide.
