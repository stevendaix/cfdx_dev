# SU2 reference dossier

## Evidence scope

SU2 is an open-source finite-volume CFD code. N14 uses official SU2 technical documentation/publications and, where necessary, source-level inspection.

## Relevant documented capabilities

The SU2 code-structure material documents convective schemes including Roe, JST, AUSM, HLLC and related flux families, plus gradient/viscous-term infrastructure. SU2 technical literature documents Green–Gauss and weighted least-squares gradients and higher-order/limited spatial reconstruction.

## N14 interpretation

SU2 is a useful reference for contrasting cell/node-centred finite-volume reconstruction, edge-based flux evaluation and compressible upwind flux families with CFDX's cell-centred FVM contracts.

A SU2 Roe/JST/AUSM capability is not an equivalent of CFDX upwind/MUSCL unless the governing flux, state reconstruction, wave-speed treatment and limiter are shown to match.

## Sources

- https://su2code.github.io/documents/Introduction_to_the_SU2_Code_Structure.pdf
- https://su2code.github.io/documents/SU2_AIAA_SciTech2014.pdf

## Limitation

The N14 matrix intentionally does not claim equivalence for compressible fluxes that CFDX does not currently implement as the same mathematical operator.
