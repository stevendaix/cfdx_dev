# CFDX user workflow

## 1. Setup

Define geometry, mesh, physical properties, models, boundary and initial conditions. Record units and reference values.

## 2. Numerics

Choose spatial discretisation, gradient/reconstruction, convection and diffusion treatment, time integration, pressure–velocity coupling, linear solver and preconditioner. Numerical choices should be explicit in the case.

## 3. Run

Execute with a declared stopping policy. Monitor residuals **and** physical invariants: conservation defects, boundedness/positivity, pressure compatibility, nonlinear convergence and quantities of interest.

A residual threshold alone is insufficient evidence of a correct CFD solution.

## 4. Restart

A restart must reproduce the intended state without silently changing the governing setup or numerical interpretation. Retain checkpoint revision and provenance.

## 5. Post-process

Inspect fields, integrated quantities, probes and time histories. For validation, compute the documented metric against its declared reference oracle.

## 6. Report

Retain software revision, case, mesh, physics, numerics, solver criteria, metric, criterion, environment and artifacts.

## 7. Promote maturity

Use the hierarchy:

**code verification → solution verification → validation → qualification**.

Do not promote a capability because a single demonstration case converges.
