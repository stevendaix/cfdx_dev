# CFDX Theory — Numerical Methods Course

CFDX Theory is intended to read as a coherent CFD and numerical-analysis course, not as an API manual.

Each chapter follows the same scientific progression:

**physical problem → mathematical formulation → discretisation → numerical properties → implementation consequences → verification → executable experiment → references**

## Course map

| Chapter | Topic | Status |
|---|---|---|
| 00 | Foundations | TOPO |
| 01 | Conservation laws | TOPO |
| 02 | Finite-volume method | TOPO |
| 03 | Meshes | TOPO |
| 04 | Gradients and reconstruction | PILOT COMPLETE |
| 05 | Fluxes | TOPO |
| 06 | Time integration | TOPO |
| 07 | Pressure–velocity coupling | TOPO |
| 08 | Linear algebra | TOPO |
| 09 | AMG / MGR / Schur | TOPO |
| 10 | Turbulence | TOPO |
| 11 | Heat transfer | TOPO |
| 12 | Radiation | TOPO |
| 13 | Multiphysics | TOPO |
| 14 | Numerical analysis | TOPO |
| 15 | Verification and validation | TOPO |

**TOPO means that the subject and intended contents have been identified, not that the CFDX capability is implemented.**

## Complete chapter template

Every mature chapter should contain:

1. Motivation
2. Physical problem
3. Mathematical formulation
4. Discretisation
5. Numerical properties
6. Accuracy
7. Stability
8. Conservation
9. Limitations
10. CFDX implementation consequences
11. Verification methodology
12. Executable numerical experiment
13. References

## Pilot

[Gradients and Reconstruction](04-gradients-reconstruction/README) is the reference format. It contains a complete 17-section course chapter and an independent executable Python study.

The pilot is the template to use before expanding the remaining chapters.

## Source-of-truth rule

Theory explains mathematics and methods. It does not own implementation status, qualification status or generated numerical results. Those remain in Developer/V&V authoritative sources.
