# VMFL036 — Literature-aligned sphere case at Re=100

## Purpose

This is the clean CFDX reference configuration for the laminar flow past a sphere benchmark.

- sphere diameter: D = 1 m
- density: rho = 1 kg/m^3
- inlet velocity: U = 1 m/s
- dynamic viscosity: mu = 0.01 Pa.s
- Reynolds number: Re = rho U D / mu = 100
- reference projected area: A = pi D^2 / 4
- reference drag coefficient: Cd = 1.0895
- flow model: steady, incompressible, laminar
- geometry: 2-D axisymmetric sphere section
- outer circular domain radius: 50 D

The Cd=1.0895 datum is the literature anchor associated with Tabata & Itakura at Re=100 and cited by Ansys for VMFL036. The published paper covers sphere drag coefficients over Re=10–200.

## Important distinction from the Ansys VMFL036 page

The current Ansys VMFL036 page lists rho=1 kg/m^3, U=1 m/s, D=1 m and mu=0.02 Pa.s. Those displayed quantities imply Re=50, while the page reports Cd=1.0895, the Re=100 literature datum.

CFDX must not hide this discrepancy by widening a validation tolerance.

This file therefore defines the literature-aligned Re=100 case. It does not claim byte-for-byte identity with the proprietary Fluent project file VMFL036_FLUENT.cas.h5. Exact Fluent-project traceability remains tracked in Issue #440.

## Geometry

- sphere radius: R=0.5 m
- sphere surface: no-slip wall
- axis: axisymmetric boundary
- outer circular boundary: radius=50 m
- use a body-fitted mesh around the sphere
- resolve the separated wake; do not replace the axisymmetric case with an arbitrary Cartesian-box surrogate

## Quantity of interest

The primary quantity is total streamwise drag coefficient:

Cd = F_D / (0.5 rho U^2 A), with A = pi D^2 / 4.

Wall force must be obtained from the traction integral F = integral_S (-p n + tau.n) dS and decomposed into pressure drag, viscous drag and total drag.

## Acceptance status

This is a case definition/reference specification, not a solver PASS.

A solver-level PASS requires actual CFDX execution, convergence and conservation evidence, independent wall-force integration, coarse/medium/fine mesh refinement, quantitative comparison against Cd=1.0895, CTest/CI execution, and reproducible logs/report.

No numerical Cd tolerance is frozen by this reference-definition change.
