# CFDX Core Skills

Phase B provides concrete reusable, vendor-neutral skills used by the CFDX agents.

## Core skills

- `repository-navigation` — establish authoritative repository context.
- `code-search` — layered textual/structural search.
- `cpp-development` — C++20 implementation and audit discipline.
- `python-development` — Python tooling and binding discipline.
- `cmake-build` — configure/build/diagnose workflow.
- `testing` — test selection, execution and evidence classification.
- `ci-analysis` — GitHub Actions diagnosis.
- `code-review` — independent change/evidence review.

## Skill composition

Core skills are deliberately small and composable. Agents may combine them with later CFD/numerical skills. Skills do not own repository state and do not replace the canonical source tree.

## Evidence

Every skill must distinguish implementation, testing, verification, validation and qualification. Numerical claims require the applicable numerical/V&V evidence rather than inference from software behavior alone.

## Next

The next skill layer is CFD/numerical: FVM, gradients/reconstruction, mesh, boundary conditions, interpolation/fluxes, pressure-velocity coupling, temporal discretisation, linear/nonlinear solvers, turbulence, wall distance, multiphysics, conservation/boundedness and validation.
