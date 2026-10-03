# Temporal discretisation

## Purpose

Guide CFDX transient-solver implementation and verification.

## Required reasoning

1. Identify the time integration formula and its formal order.
2. Separate spatial and temporal truncation errors.
3. State explicit/implicit treatment and stability restrictions.
4. Check initialisation/startup treatment for multistep methods.
5. Check nonlinear iteration and pseudo-time coupling separately from physical time accuracy.
6. Define restart/checkpoint semantics and ensure time/state consistency.

## Verification

Use manufactured transient solutions, spectral or modal problems with known decay rates, periodic advection, and Richardson refinement. Demonstrate temporal order only after spatial error is controlled or analytically removed.

## Anti-patterns

- Do not infer temporal order from one time step.
- Do not mix spatial and temporal error in an order claim.
- Do not alter timestep or iteration stopping criteria between refinement levels without documenting the effect.
