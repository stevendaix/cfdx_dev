# CFDX Theory

The Theory documentation is a structured CFD and numerical-method course for CFDX.

Its purpose is to explain the mathematics, numerical analysis and physical assumptions behind the methods implemented by CFDX. It is not a list of API calls and it is not a copy of the source code.

## Course structure

```text
00-foundations
01-conservation-laws
02-finite-volume-method
03-meshes
04-gradients-reconstruction
05-fluxes
06-time-integration
07-pressure-velocity-coupling
08-linear-algebra
09-amg-mgr-schur
10-turbulence
11-heat-transfer
12-radiation
13-multiphysics
14-numerical-analysis
15-verification-validation
```

Each chapter may contain several lessons. The directory number defines conceptual order, not implementation maturity.

## Standard chapter contract

A substantive theory chapter should normally contain:

1. Motivation and physical problem.
2. Mathematical formulation.
3. Assumptions and scope.
4. Discretisation.
5. Conservation properties.
6. Consistency, stability and accuracy.
7. Error sources and failure modes.
8. Numerical examples.
9. CFDX implementation consequences.
10. Verification methodology.
11. References.

When a chapter contains an executable numerical experiment, the experiment must identify the mathematical problem, reference solution, mesh/discretisation, numerical configuration, measured quantities, expected behaviour and reproducibility command.

## Executable scientific content

Python is the preferred execution language for theory demonstrations. The long-term model is:

```text
theory source
    -> Python executable experiment
    -> machine-readable result
    -> generated figure/table
    -> web rendering
    -> Typst publication
```

Committed notebooks are not the default. Prefer Python source files that can be executed in CI and rendered into the documentation.

## Rules

- Do not claim a numerical property from a plot alone.
- Do not claim an order of accuracy without a defined refinement study.
- Do not hide failed experiments.
- Do not hard-code a result merely to reproduce a figure.
- Do not duplicate implementation contracts from Developer documentation.
- Link every non-trivial theoretical assertion to the authoritative reference.
