# CFDX Theory — Numerical Methods Course

CFDX Theory is the scientific course and mathematical reference for the code. It is separate from Developer implementation contracts and V&V status.

Every mature chapter follows:

**physical problem → governing equations → assumptions → finite-volume discretisation → numerical properties → CFDX implementation → verification → benchmark → limitations → improvement paths → references**

## Course map

| Chapter | Topic | Status |
|---|---|---|
| 00 | Foundations | IN PROGRESS |
| 01 | Conservation laws | IN PROGRESS |
| 02 | Finite-volume method | IN PROGRESS |
| 03 | Meshes | REPOSITORY-GROUNDED |
| 04 | Gradients and reconstruction | PILOT COMPLETE |
| 05 | Fluxes | REPOSITORY-GROUNDED |
| 06 | Time integration | REPOSITORY-GROUNDED |
| 07 | Pressure–velocity coupling | REPOSITORY-GROUNDED |
| 08 | Linear algebra | REPOSITORY-GROUNDED |
| 09 | AMG / MGR / Schur | REPOSITORY-GROUNDED |
| 10 | Turbulence | REPOSITORY-GROUNDED |
| 11 | Heat transfer | REPOSITORY-GROUNDED |
| 12 | Radiation | REPOSITORY-GROUNDED |
| 13 | Multiphysics | REPOSITORY-GROUNDED |
| 14 | Numerical analysis | REPOSITORY-GROUNDED |
| 15 | Verification and validation | REPOSITORY-GROUNDED |
| 16 | Data model and file formats | CODE-AUDITED |
| 17 | Computational chain and code map | CODE-AUDITED |
| 18 | Source-tree physics audit | AUDIT REGISTER |

These labels describe documentation maturity, not numerical qualification.

## Mandatory equation page

Every important equation states its physical origin, assumptions, symbols and units, sign convention, continuous equation, control-volume integral, discrete approximation, implementation path, tests, benchmark or oracle, limitations, improvement path and references.

## Status vocabulary

- Implemented: a code path exists.
- Verified: a defined mathematical property has executable evidence.
- Qualified: the declared V&V population supports the intended scope.

Theory never promotes a capability merely because an implementation exists.
