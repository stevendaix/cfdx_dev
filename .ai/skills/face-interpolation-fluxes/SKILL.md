# Face interpolation and fluxes

## Purpose

Audit the separation between reconstruction, face interpolation, numerical fluxes, and assembled FVM contributions.

## Required reasoning

1. Identify cell-to-face interpolation independently of the gradient implementation.
2. State whether interpolation is central, upwind, limited, bounded, or otherwise specialised.
3. Trace face values into mass, momentum, energy, diffusion, and other fluxes.
4. Check orientation and conservation across internal faces.
5. Check boundedness/positivity requirements for quantities where the scheme requires them.
6. Keep flux reconstruction independent from the gradient implementation unless the method explicitly couples them.

## Verification

Use constant preservation, linear-field tests, symmetry, flux antisymmetry, conservation, and boundedness/positivity checks. For higher-order schemes, use refinement studies and known reference solutions.

## Anti-patterns

- Do not equate a small residual with flux correctness.
- Do not bury limiter behaviour inside unrelated gradient code.
- Do not use a single benchmark to qualify every interpolation/flux path.
