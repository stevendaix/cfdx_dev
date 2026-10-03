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

## Domain numerical skills

- `turbulence-modelling` — turbulence equations, closure terms, wall treatment and model evidence.
- `wall-distance` — distance-field accuracy, monotonicity, positivity and reconstruction evidence.
- `multiphysics-coupling` — coupled-field equations, exchange conservation and coupling verification.
- `conservation-boundedness` — discrete conservation, flux antisymmetry, positivity and boundedness.

## Next

Further skills will cover V&V workflows, code intelligence, MCP operation, GUI/TUI and runtime/development agents.


## V&V workflow skills

- `vv-workflow` — reproducible verification/validation workflow
- `mms-verification` — manufactured-solution verification
- `refinement-convergence` — mesh/time refinement and observed-order analysis
- `benchmark-reference` — benchmark and independent-reference comparison
- `qualification-matrix` — auditable maturity/qualification mapping
