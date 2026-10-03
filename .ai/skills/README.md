# CFDX Skills

The `.ai/skills/` tree contains small, vendor-neutral, composable skills used by CFDX agents.

## Core engineering skills

- `repository-navigation` — establish authoritative repository context.
- `code-search` — layered textual/structural search.
- `cpp-development` — C++20 implementation and audit discipline.
- `python-development` — Python tooling and binding discipline.
- `cmake-build` — configure/build/diagnose workflow.
- `testing` — test selection, execution and evidence classification.
- `ci-analysis` — GitHub Actions diagnosis.
- `code-review` — independent change/evidence review.

## Foundational CFD/numerical skills

- `fvm-fundamentals` — conservative finite-volume formulation and discrete balance.
- `mesh-topology-quality` — topology, geometry, metrics and mesh-quality reasoning.
- `gradients-reconstruction` — Green–Gauss, least-squares, weighted LS, conditioning and reconstruction separation.
- `boundary-conditions` — mathematical BC semantics and discrete boundary contributions.
- `face-interpolation-fluxes` — face values, numerical fluxes, conservation and boundedness.
- `numerical-evidence` — implementation → testing → verification → validation → qualification evidence ladder.

## Solver and time skills

- `pressure-velocity-coupling` — pressure correction, continuity, flux correction, nullspaces and coupled formulations.
- `temporal-discretisation` — physical/pseudo-time schemes, stability, restart and temporal-order verification.
- `linear-nonlinear-solvers` — matrix properties, Krylov/preconditioning, nonlinear convergence and stopping criteria.
- `amg-schur-preconditioning` — block Schur complements, AMG hierarchy, Galerkin consistency and energy/error evidence.

## Composition and evidence

Skills are deliberately small and composable. Agents may combine them for a task instead of relying on a monolithic prompt. Skills do not own repository state and do not replace the canonical CFDX source tree.

Every numerical skill must distinguish implementation, testing, verification, validation and qualification. Numerical claims require explicit evidence: method/configuration, mesh/refinement, metric/reference, tolerance, commit/PR and limitations.

## Next

Further CFD/numerical skills will cover turbulence, wall distance, multiphysics, conservation/boundedness and V&V workflows.
