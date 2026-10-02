# CFDX Validation and V&V Evidence

**Status: evidence-store framework — quantitative results are populated only from executed campaigns.**

This tree stores executed evidence. It is deliberately separate from `docs/vv/`, which defines the methodology and acceptance rules.

## Rule

A report may contain a PASS, FAIL, NOT QUALIFIED or NOT RUN result only when that state is backed by an identifiable execution artifact. Never invent or copy CFDX results from a reference source.

## Campaign classes

- code verification
- solution verification
- MMS
- spatial and temporal convergence
- conservation and boundedness
- linear algebra
- pressure–velocity coupling
- gradients, diffusion and convection
- turbulence
- thermal and radiation
- multiphysics
- benchmark validation
- qualification

## Required campaign record

Revision, case/setup, mesh, physics, numerics, solver criteria, oracle, metric, acceptance criterion, raw result, derived result, environment, status and limitations.

The benchmark directory is the validation phase. A benchmark result does not replace lower-level verification evidence.
