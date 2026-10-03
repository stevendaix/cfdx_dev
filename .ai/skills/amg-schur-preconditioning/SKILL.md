# AMG, Schur complements and preconditioning

## Purpose

Guide CFDX block-preconditioner and multigrid development.

## Required reasoning

1. State the block matrix and ordering.
2. Define the Schur complement approximation explicitly.
3. Distinguish exact Schur oracle checks from production approximations.
4. For AMG, identify hierarchy construction, coarsening, interpolation, restriction, smoother and cycle.
5. Check Galerkin consistency where applicable.
6. Evaluate convergence using both residual and an appropriate error/energy measure.
7. Report setup cost, iteration count and robustness across mesh refinement.

## Verification

Use small exact block systems, manufactured operators, spectral/energy diagnostics, and mesh-refinement studies. Compare approximate Schur/preconditioners against direct or exact-oracle references when possible.

## Anti-patterns

- Do not claim AMG quality from a single residual plot.
- Do not replace a failed convergence gate by increasing iteration limits indefinitely.
- Do not silently fall back from GPU to CPU or from the requested preconditioner to another method.
