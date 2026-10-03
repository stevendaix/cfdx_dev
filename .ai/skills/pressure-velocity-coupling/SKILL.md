# Pressure–velocity coupling

## Purpose

Guide implementation and audit of pressure–velocity coupling algorithms in CFDX.

## Required reasoning

1. Identify the segregated or coupled formulation and the discrete continuity/momentum operators it couples.
2. State the pressure equation or block system actually assembled.
3. Trace Rhie–Chow or equivalent face-flux treatment when collocated storage is used.
4. Distinguish SIMPLE/SIMPLEC, PISO, PIMPLE, fractional-step and fully coupled formulations.
5. Check pressure null-space handling and reference-pressure semantics.
6. Check conservation of corrected face mass fluxes.
7. For block preconditioners, preserve the mathematical block structure and document approximations.

## Verification

Use continuity closure, constant-pressure/velocity preservation, manufactured solutions, canonical laminar cases, iteration/convergence behaviour and mesh refinement. For coupled solvers, compare against an independent block or direct reference where practical.

## Anti-patterns

- Do not call a pressure correction converged because momentum residuals alone are small.
- Do not hide pressure null-space problems with arbitrary offsets.
- Do not silently switch coupling algorithms or fallback to a different solver.
- Do not qualify a coupling algorithm without checking mass conservation and solution accuracy.
