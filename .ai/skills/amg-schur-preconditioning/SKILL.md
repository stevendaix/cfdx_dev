# Skill: AMG and Schur Preconditioning

## Purpose
Guide algebraic multigrid and block-Schur preconditioning work for CFDX.

## Scope
- Block systems and velocity/pressure Schur complements.
- Exact Schur or small-system oracles.
- Approximate Schur operators.
- AMG hierarchy, coarsening, interpolation/restriction, smoothing and cycles.
- Galerkin consistency and energy/error measures.

## Method
1. Write the block matrix explicitly and define every block.
2. Establish a trusted exact or direct oracle on controlled systems.
3. Define the Schur complement mathematically before approximating it.
4. For AMG, document strength/coarsening, transfer operators, coarse operator construction and smoother.
5. Verify Galerkin consistency where applicable.
6. Measure error/energy contraction in addition to residual reduction.
7. Repeat tests under mesh refinement and representative coefficient variation.

## Required evidence
- Exact-Schur comparison on small controlled systems.
- AMG hierarchy/operator consistency checks.
- Energy/error contraction evidence.
- Refinement robustness.
- Diagnostics exposing breakdown, rank deficiency or poor conditioning.

## Review rules
Do not call a V-cycle robust because it reduces one residual. Do not silently replace an unavailable accelerator/preconditioner with another implementation. Any fallback must be explicit and tested.
