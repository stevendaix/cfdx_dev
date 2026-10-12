# Boundary conditions

## Purpose

Reason about CFDX boundary-condition semantics and their discrete contribution.

## Required reasoning

1. Identify the physical condition and its mathematical form.
2. Map it to the CFDX boundary/BC/value-provider architecture.
3. Identify which face fluxes, coefficients, or source terms are modified.
4. Check consistency with the interior discretisation.
5. Check conservation implications and signs at boundaries.
6. Verify transient and steady semantics separately when relevant.
7. Ensure no legacy automatic mapping or hidden fallback changes the requested condition.

## Verification

Use canonical Dirichlet, Neumann, mixed/Robin, symmetry, wall, inlet/outlet, and periodic cases as applicable. Prefer analytical or manufactured solutions with explicit expected boundary traces and fluxes.

## Anti-patterns

- Do not infer a boundary condition from a variable name alone.
- Do not silently replace unsupported BCs with another type.
- Do not qualify a BC merely because the solver runs.
