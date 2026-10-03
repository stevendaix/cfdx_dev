# CFDX Verification & Validation

**Status: PARTIAL — governance is established; evidence migration and campaign execution are ongoing.**

The V&V Guide is the operational evidence manual. Theory explains mathematics; Developer explains implementation; V&V defines how claims are demonstrated and retained.

```{toctree}
:maxdepth: 2

00-vv-governance/chapter
01-requirements-and-claims/chapter
02-code-verification/chapter
03-solution-verification/chapter
04-mms/chapter
05-spatial-convergence/chapter
06-temporal-convergence/chapter
07-conservation/chapter
08-boundedness/chapter
09-linear-solver-verification/chapter
10-pressure-velocity-verification/chapter
11-gradient-verification/chapter
12-diffusion-verification/chapter
13-convection-verification/chapter
14-turbulence-verification/chapter
15-thermal-verification/chapter
16-radiation-verification/chapter
17-multiphysics-verification/chapter
18-benchmark-validation/chapter
19-qualification-matrix/chapter
20-evidence-and-reproducibility/chapter
```

## Evidence vocabulary

- **Code verification:** implementation matches the specified mathematical algorithm.
- **Solution verification:** discretisation and iteration errors are controlled for the declared calculation.
- **Validation:** physical/reference comparison with stated uncertainty.
- **Qualification:** an explicitly bounded capability population satisfies declared acceptance criteria.

A green CI job is not automatically validation or qualification.

## Minimum evidence record

Every quantitative claim should identify, where applicable: software revision, case/setup revision, mesh and geometry, physics/model choices, numerical schemes, solver/stopping criteria, metric and oracle, acceptance criterion, execution environment and retained artifacts.

`docs/validation/` is the executed evidence store. Never substitute a reference result for a CFDX result.
