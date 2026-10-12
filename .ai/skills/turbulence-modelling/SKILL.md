# Skill: Turbulence Modelling

## Purpose
Guide CFDX agents through turbulence-model implementation, review and verification without conflating model equations with solver convergence.

## Scope
- RANS: Spalart–Allmaras and k–omega SST.
- LES/DES-family extensions where explicitly implemented.
- Turbulent transport, production/destruction, wall treatment and limiter interactions.
- Model constants, dimensional consistency and boundary conditions.

## Method
1. State the continuous model and closure assumptions.
2. Map every model term to the discrete CFDX operators.
3. Check dimensions, signs, source linearisation and boundary contributions.
4. Verify limiting/degenerate states and positivity-sensitive quantities.
5. Separate turbulence-model verification from mean-flow solver convergence.
6. Compare against analytical/benchmark/reference-code evidence only with configuration and model details documented.

## Required evidence
- Unit/dimensional checks.
- Constant/zero-turbulence limiting cases where mathematically applicable.
- Manufactured or controlled transport tests.
- Canonical benchmark evidence for the implemented model.
- Mesh/time refinement where order claims are made.

## Review rules
Do not claim a turbulence model is validated merely because a flow converges. Do not alter model constants or clipping solely to obtain benchmark agreement without documenting the mathematical and physical rationale.
