# Skill: Conservation and Boundedness

## Purpose
Guide verification of discrete conservation, flux antisymmetry, positivity and boundedness.

## Method
1. Identify the conserved quantity and its discrete balance.
2. Trace every internal-face contribution and verify antisymmetry.
3. Check boundary fluxes independently.
4. Distinguish local cell balance from global balance.
5. For bounded variables, identify the mathematical mechanism enforcing the bound.
6. Test positivity/boundedness independently of residual convergence.

## Required evidence
- Internal-face flux antisymmetry.
- Global and per-cell balance.
- Boundary flux accounting.
- Constant-state preservation.
- Positivity/boundedness tests on representative cases.
- Diagnostics that identify the first violating cell/face and contribution.

## Review rules
Residual convergence does not prove conservation or boundedness. Never suppress a violating case or relax bounds merely to obtain a green test.
