# Ansys Fluent reference dossier

## Evidence scope

Fluent evidence is based on current public Ansys documentation. Proprietary implementation details not exposed by the documentation are not inferred.

## Relevant documented capabilities

Fluent documents Green–Gauss cell-based, Green–Gauss node-based and least-squares cell-based gradient reconstruction. Its pressure-based solver documentation also exposes SIMPLE, SIMPLEC, PISO, Fractional Step and Coupled pressure-velocity families. First/second-order spatial options and selected higher-order options such as QUICK/MUSCL are documented.

Fluent documentation also explicitly discusses boundary treatment of least-squares gradients, demonstrating why N14 must include boundary reconstruction as part of the formulation comparison.

## Sources

- https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/flu_ug/flu_ug_sec_solve_check.html
- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_sec_overset.html
- https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_chp_solve.html
- https://innovationspace.ansys.com/knowledge/forums/topic/i-computed-the-species-mass-fraction-in-a-setup-where-the-convective-velocity-at-the-inlets-outlets-is-negligible-that-is-the-species-transport-over-inlets-and-outlets-is-diffusion-dominated-i-acti/

## N14 interpretation

The Fluent name “Least Squares Cell Based” does not establish equality with CFDX WLS. Boundary correction mode, weighting, stencil and face reconstruction must be compared explicitly.
