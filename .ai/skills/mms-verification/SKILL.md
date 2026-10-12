---
name: cfdx-mms-verification
description: Use manufactured solutions to isolate spatial, temporal and operator verification errors.
version: 1
---

# Purpose

Provide a rigorous MMS workflow for verifying individual PDE operators and coupled solver behavior.

# Procedure

- State the exact manufactured field and PDE/source construction.
- Verify boundary and initial conditions match the manufactured solution.
- Separate spatial and temporal error by holding one contribution sufficiently small or using an appropriate refinement design.
- Measure norms with fixed definitions and units.
- Use at least three meaningful refinement levels when estimating an observed order.
- Inspect monotonicity and asymptotic behavior rather than relying on a fitted slope alone.
- Attribute observed order loss to the relevant operator, boundary treatment, reconstruction, mesh or time integration where evidence permits.

# Evidence

Report mesh/time levels, errors, observed orders, solver tolerances and whether the solution is in the asymptotic regime.

# Failure handling

Do not manufacture a passing order by discarding unfavorable refinement levels without documenting and justifying the exclusion.
