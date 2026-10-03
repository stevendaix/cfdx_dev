# Linear and nonlinear solvers

## Purpose

Guide implementation and audit of CFDX linear/nonlinear solver stacks.

## Required reasoning

1. Identify the matrix/operator, unknown ordering and symmetry properties.
2. Identify Krylov method, smoother/preconditioner and stopping criterion.
3. Distinguish algebraic residual reduction from error reduction.
4. Check scaling, conditioning and null spaces.
5. For nonlinear solves, distinguish nonlinear residual convergence from discretisation error.
6. Record absolute and relative tolerances and iteration limits explicitly.

## Verification

Use manufactured linear systems, known spectra where available, exact/direct references, iteration histories, residual/error relations, and mesh-size scaling. For nonlinear systems, use controlled problems with known solutions.

## Anti-patterns

- Do not use residual reduction alone as an accuracy claim.
- Do not inflate tolerances to make a solver appear robust.
- Do not silently change preconditioners or solver classes.
