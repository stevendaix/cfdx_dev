# Skill: Linear and Nonlinear Solvers

## Purpose
Provide a rigorous workflow for linear-system, nonlinear-iteration, convergence, and stopping-criterion work.

## Scope
- Krylov methods and stationary iterations.
- Preconditioners and scaling.
- Nonlinear fixed-point/Newton-like iterations.
- Residual norms, error estimates, stagnation and termination.
- Matrix properties: symmetry, definiteness, conditioning and nullspaces.

## Method
1. Characterise the actual assembled operator before choosing a solver.
2. State matrix ordering, block structure and expected algebraic properties.
3. Separate absolute/relative residual criteria from solution-error estimates.
4. Check scaling and nullspaces explicitly.
5. Distinguish linear solver error, nonlinear iteration error and discretisation error.
6. Record iteration counts and convergence histories for benchmark cases.

## Required evidence
- Known small-system or manufactured oracle where possible.
- Solver convergence history.
- Robustness across mesh refinement.
- Tests for singular/nullspace systems when mathematically required.
- Failure diagnostics for divergence/stagnation.

## Review rules
Residual reduction is evidence about the algebraic solve, not automatically about physical accuracy. Never relax tolerances solely to make a benchmark pass.
