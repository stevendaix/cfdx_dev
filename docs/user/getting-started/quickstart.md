# Quickstart

## Workflow

The normal CFDX workflow is:

1. define the physical problem;
2. select the mathematical model and closures;
3. prepare or import the mesh;
4. define fields and boundary conditions;
5. choose spatial, temporal, coupling and linear-solver methods;
6. write the case/setup representation;
7. run code and solution verification;
8. execute the relevant benchmark or validation case;
9. inspect residuals **and** physical invariants;
10. retain the evidence needed for reproducibility.

## First calculation

Use a small canonical case such as Couette flow or Poiseuille flow when learning the numerical workflow. These cases have analytical reference solutions and isolate different parts of the discretisation.

Do not treat a successful process exit as proof of physical correctness. Compare the declared metric with its oracle and acceptance criterion.

## Failure triage

When a calculation fails, localise the first failed invariant:

**mesh → geometry → boundary condition → operator/flux → algebraic system → nonlinear coupling → linear solve → physical metric**.

Fix the earliest failed contract before tuning later solver stages.
