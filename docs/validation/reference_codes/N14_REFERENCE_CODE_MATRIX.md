# N14 — Reference-code comparison matrix

## Scope

This matrix compares representative CFDX numerical capabilities against publicly documented capabilities of five reference codes. It is deliberately method-level rather than feature-name-level.

### Classification

| Classification | Meaning |
|---|---|
| SAME MATHEMATICAL METHOD | The public formulation is equivalent at the mathematical-contract level, subject to explicitly listed conventions. |
| SIMILAR METHOD | Same numerical family/objective, but a material formulation, stencil, limiter, coupling or implementation difference remains. |
| DIFFERENT METHOD / SAME PURPOSE | Different mathematical algorithm addressing the same numerical task. |
| UNAVAILABLE IN CFDX | The reference capability is documented, but CFDX does not currently provide that capability as a production method. |

### CFDX maturity is not implied by comparison

Each row carries independent CFDX statuses:
`implementation_status`, `verification_status`, `validation_status`, `qualification_status`.

## Comparison matrix

| Capability | CFDX formulation | OpenFOAM | SU2 | Fluent | CFX | STAR-CCM+ |
|---|---|---|---|---|---|---|
| Cell-centred least-squares gradient | weighted/least-squares reconstruction with explicit rank/conditioning policy | Similar: least-squares gradient families; stencil choices differ | Similar: weighted LS gradients documented for node-based FV | Similar: LS cell-based gradient | Different purpose: CFX documents nodal/control-volume gradients | Similar: FV gradient/reconstruction is documented, exact proprietary contract not public |
| Green–Gauss gradient | point/cell Green–Gauss variants | Same family; exact stencil/options differ | Same family available in SU2 literature/code | Same family: cell/node variants | Similar mathematical family but CFX uses integration-point/nodal formulation | Similar family; exact public implementation details limited |
| Non-orthogonal diffusion | orthogonal term plus explicit correction policy | Similar; corrected/limited/non-orthogonal schemes documented | Similar FV viscous reconstruction, formulation differs | Similar finite-volume diffusion treatment | Different implementation details; CFX uses its own integration-point formulation | Similar FV diffusion, exact implementation details not public |
| First-order upwind convection | upwind face state | Same mathematical method | Same purpose; SU2 flux families are generally node/edge based | Same mathematical method | Same mathematical method | Same mathematical method |
| Second-order upwind | reconstructed upstream face state | Similar/available family; exact reconstruction differs by scheme | Similar higher-order reconstruction/limiter family | Similar second-order discretisation | Similar advection correction/blend formulation | Same named family; exact reconstruction differs |
| MUSCL | slope/gradient-based high-resolution reconstruction | Similar family | Similar; limiter-based higher-order reconstruction | Similar; MUSCL documented | Different method/same purpose unless exact option is configured | Similar; third-order MUSCL/CD hybrid is documented |
| TVD limiters | MinMod, Van Leer, Van Albada, Superbee, MC | Similar limiter families | Similar limiter families; specific defaults differ | Similar limiter families | Similar bounded high-resolution philosophy; not assumed identical | Similar high-resolution/limiter families |
| QUICK / quadratic upwind | quadratic upwind face reconstruction | Similar higher-order FV family; exact implementation differs | Different method/same purpose in the documented SU2 baseline catalogue | Same family is documented | Different method/same purpose in baseline CFX advection options | Different method/same purpose unless exact second-order-upwind variant is shown equivalent |
| SIMPLE | segregated pressure-velocity algorithm | Same named algorithm, but coefficient/pressure correction contract must still be compared | Different method/same purpose in typical SU2 coupled formulation | Same named family | Different method/same purpose for CFX's coupled/segregated controls | Same named algorithm is documented |
| SIMPLEC | CFDX segregated SIMPLEC formulation | Same named family; algebraic coefficient conventions must be checked | Different method/same purpose | Same named family | Different method/same purpose | Same named family is documented |
| PISO | pressure-correction predictor/corrector | Same named family | Different method/same purpose in documented baseline | Same named family | Different method/same purpose | Same named family |
| PIMPLE | CFDX PIMPLE implementation | Same named family in OpenFOAM | Different method/same purpose | Different method/same purpose | Different method/same purpose | Different method/same purpose |
| Coupled U-p | monolithic coupled pressure/velocity path | Similar coupled solver objective; exact block/preconditioner differs | Similar coupled FVM objective | Same purpose; Fluent coupled solver is not assumed algebraically identical | Similar coupled pressure/velocity formulation | Similar coupled flow solver |
| AMG | native AMG / smoothed aggregation paths | Similar AMG family; hierarchy/coarsening differ | Multigrid available; exact implementation differs | Algebraic multigrid available | Algebraic multigrid available | Algebraic multigrid available |
| Schur preconditioning | block-local Schur plus explicit exact-oracle infrastructure | Different method/same purpose unless algebraically matched | Different method/same purpose | Different method/same purpose | Different method/same purpose | Different method/same purpose |
| BDF2 | constant/variable-step BDF2 | Similar temporal family | Different baseline naming/options; RK/Euler documented | Similar second-order implicit family | Same mathematical family: second-order backward Euler | Similar second-order implicit family |
| Adaptive CFL | auditable adaptive controller | Similar control objective; exact controller differs | Similar CFL/time-step controls | Similar adaptive/pseudo-time controls | Similar timescale/Courant controls | Similar automatic CFL control |

## Important interpretation examples

### Least-squares is not automatically identical

Fluent explicitly documents least-squares cell-based gradients and distinguishes them from node-based Green–Gauss. SU2 documents weighted least-squares gradients at grid nodes. OpenFOAM exposes multiple least-squares/polynomial stencil families. These facts establish a common numerical family, not one common operator.

### High-resolution convection is not automatically TVD-equivalent

CFX's High Resolution scheme uses a nonlinear boundedness construction based on a local stencil. That does not make it identical to CFDX MUSCL+MC, Van Leer, Superbee, etc. A row is therefore `similar_method` unless the actual limiter and reconstruction equations match.

### Schur methods require algebraic evidence

A CFDX Schur approximation is not promoted to equivalence because another code has a “coupled” or “Schur” option. The discrete block matrices, sign convention, inner velocity approximation, null-space treatment and pressure operator must match before `same_mathematical_method` is allowed.

## Reference sources

The detailed source dossiers contain the authoritative URLs and the evidence scope for each code. Public documentation is versioned where possible; proprietary implementation details are explicitly marked unknown.
