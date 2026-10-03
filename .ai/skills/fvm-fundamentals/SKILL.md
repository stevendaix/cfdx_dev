# FVM fundamentals

## Purpose

Provide a repository-grounded skill for finite-volume discretisation and audit of CFDX operators.

## Required reasoning

1. Identify the governing integral balance and control-volume definition.
2. Separate volume sources from face fluxes.
3. State the discrete unknown location and geometric measures.
4. Verify that interior-face contributions are antisymmetric when a conservative formulation requires it.
5. Check dimensional consistency and sign conventions.
6. Distinguish exact geometric identities from approximations.
7. Trace the discrete operator into the actual CFDX implementation before making claims.

## Conservation contract

For a conservative FVM operator, an internal face shared by cells P and N must contribute equal and opposite fluxes to the two cell balances, modulo deliberate source/sink terms and floating-point effects. Boundary fluxes must be accounted for explicitly.

A residual reaching a small value is not evidence of conservation by itself. Inspect assembled coefficients, fluxes, and global/cell balance diagnostics.

## Verification

Use manufactured solutions, constant-field tests, polynomial fields where appropriate, mesh refinement, and independent integral checks. Report observed order and the norm/quantity used.

## Anti-patterns

- Do not infer conservation from residual convergence alone.
- Do not silently change signs or normal conventions to make a test pass.
- Do not weaken tolerances or remove failing cases to manufacture qualification.
- Do not call a partial operator fully qualified.

## Evidence

Cite the relevant CFDX source, test, benchmark, and authoritative numerical reference when making a numerical claim.
