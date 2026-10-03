# OpenFOAM reference dossier

## Evidence scope

OpenFOAM is open source, so the strongest N14 evidence can combine official documentation and inspectable source. The comparison must still identify the exact release/version because numerical options evolve.

## Relevant documented capabilities

- Surface-normal gradient schemes include orthogonal, corrected and limited-corrected formulations.
- Polynomial-fit gradient/surface-normal schemes include `quadraticFit` and `linearFit`.
- Least-squares and extended-stencil approaches are documented for distorted meshes.
- SIMPLE/PIMPLE solver-control families are part of the documented OpenFOAM workflow.

## N14 interpretation

OpenFOAM is particularly useful for comparing FVM operator contracts, stencil construction and non-orthogonal correction. Matching a scheme name is insufficient: the exact stencil, weighting, face reconstruction and correction coefficient must be compared.

## Sources

- https://openfoam.org/release/2-2-1/
- https://openfoam.org/release/2-3-0/numerics/
- https://openfoam.org/release/2-0-0/run-time-control-code-compilation/

## Limitation

The cited pages are official historical documentation used to establish mathematical capability. They are not claims that current CFDX must reproduce a specific OpenFOAM version byte-for-byte.
