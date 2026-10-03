# CFDX concepts

## Case, state and output

CFDX separates **case/setup**, **state/checkpoint**, and **output**. This prevents a restart checkpoint from becoming an implicit configuration format.

## Numerical chain

**physical problem → mathematical model → PDE → integral balance → mesh/control volumes → discrete operators → algebraic system → nonlinear iterations → linear solves → convergence diagnostics → V&V**.

Each transition has contracts that can be verified independently.

## Maturity vocabulary

- **Implemented:** a code path exists.
- **Verified:** a specified mathematical/software property has executable evidence.
- **Validated:** independent physical/reference evidence supports the declared case.
- **Qualified:** the declared capability population and acceptance gates are complete.

These labels are not interchangeable.

## Boundary conditions

Boundary conditions contribute algebraically to the governing equations. Their physical type, sign convention, units and numerical treatment must be documented; an API spelling does not by itself establish mathematical correctness.

## Parallel and accelerator execution

Parallel execution must preserve numerical contracts across partition interfaces. GPU execution is an explicit execution path; failure to use the requested accelerator must not be silently presented as successful GPU execution.
