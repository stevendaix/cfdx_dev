# Skill: Pressure–Velocity Coupling

## Purpose
Provide a repository-grounded method for designing, implementing, reviewing, and verifying pressure–velocity coupling for incompressible or low-Mach finite-volume solvers.

## Scope
- Segregated methods: SIMPLE, SIMPLEC, PISO, PIMPLE and fractional-step variants.
- Fully coupled/block formulations.
- Pressure correction, continuity enforcement, pressure nullspaces and reference pressure.
- Face mass-flux correction and discrete conservation.
- Coupling diagnostics and block-preconditioning interfaces.

## Method
1. State the continuous equations and assumptions.
2. Identify the exact discrete momentum and continuity operators used by CFDX.
3. Derive the pressure equation/correction from the actual discrete fluxes; do not rely on a residual-only argument.
4. Verify that corrected face mass fluxes satisfy continuity and remain antisymmetric across internal faces.
5. Make pressure nullspace handling explicit for closed Neumann systems.
6. Distinguish algorithmic convergence from discretisation error.
7. For block solvers, document the matrix block structure and the Schur-complement approximation.

## Required evidence
- Constant-state preservation.
- Discrete continuity/flux balance.
- Pressure-nullspace/reference-pressure test.
- Manufactured or analytical pressure/velocity verification where applicable.
- Refinement evidence for the spatial scheme.
- Iteration/stopping criteria documented independently from accuracy claims.

## Review rules
- Never infer conservation from a small residual alone.
- Do not silently substitute a different coupling algorithm or CPU fallback.
- Do not hide instability by loosening tolerances or iteration limits.
