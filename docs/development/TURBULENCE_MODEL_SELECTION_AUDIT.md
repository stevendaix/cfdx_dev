# Turbulence model audit and selection architecture

## Scope

This audit compares the current CFDX turbulence stack with the model-selection structure exposed by mature CFD solvers. CFDX must expose only models whose implementation status is explicit, while keeping the case format stable as additional turbulence and transition models are implemented.

## Current CFDX state

### Solver-ready transport models
- Standard k-epsilon
- RNG k-epsilon
- Realizable k-epsilon
- k-omega
- k-omega SST
- Spalart-Allmaras

These have transport-equation implementations in the current C++ turbulence layer and regression/validation-level tests. Presence here does not mean solver-level physical validation is complete; #118 remains the qualification authority.

### Kernel-only / incomplete scale-resolving models
- Smagorinsky
- WALE
- DES
- DDES
- IDDES

The current code explicitly documents these as algebraic closures or length-scale helpers. They must not appear as production solver choices merely because the enum or kernel exists.

### Planned / not yet complete
- Dynamic one-equation LES (DYNAMIC_KEQN)
- SA negative and correction variants
- transition transport
- Reynolds-stress transport
- complete wall-function families
- complete transported LES and hybrid RANS-LES implementations

## Comparison with mature CFD solvers

| Capability | CFDX current state | SU2 | OpenFOAM | Fluent | STAR-CCM+ |
|---|---|---|---|---|---|
| SA RANS | transport baseline | SA + negative/Edwards/corrections | SA family | SA | SA |
| k-epsilon family | standard/RNG/realizable transport | not the central RANS family | broad RAS family | standard/RNG/realizable | standard/realizable + variants |
| k-omega/SST | SST + generic k-omega transport | SST variants | broad SST/RAS family | standard/BSL/SST/GEKO | standard/SST |
| LES | algebraic kernels only | available in broader solver stack | Smagorinsky/WALE/dynamic families | LES/SRS families | LES + dynamic/WALE |
| DES/DDES/IDDES | kernels only | broader turbulence/scale-resolving capabilities | DES/DDES/IDDES families | DES/SRS families | SST/SA/elliptic-blending DES |
| Transition | not implemented | Langtry-Menter family | separate models/extensions | transition SST/k-kl-omega | Gamma/Gamma-ReTheta |
| RSM | not implemented | limited compared with industrial RANS stacks | Reynolds-stress models | RSM | several RSM closures |
| Wall treatment | y+ diagnostic only | SA/SST wall functions | model-dependent wall treatment | several wall treatments | high/low/all-y+ |
| User-facing selection | basic 5-entry dropdown | explicit config options | model dictionaries | hierarchical model panel | hierarchical physics tree |

The external documentation supports these broad capability differences. SU2 documents several SA and SST variants plus Langtry-Menter transition options. OpenFOAM groups turbulence into RAS, LES and DES and exposes multiple RAS/LES/DES families. Fluent exposes SA, three k-epsilon variants, several k-omega variants, transition, RSM, SAS, DES and LES with separate near-wall treatment controls. STAR-CCM+ exposes RANS, RSM, LES, DES and transition families with explicit coupling between transition and the base turbulence model.

## Architectural conclusion

The important lesson is not to reproduce another solver catalogue. Mature solvers separate:
1. turbulence family;
2. base closure;
3. optional corrections;
4. transition model;
5. wall treatment;
6. scale-resolving/hybrid strategy;
7. applicability and dimensional constraints.

CFDX should therefore avoid a flat model=SST list as the long-term case schema.

### Target user workflow

The GUI/TUI should present:
1. Flow regime — Laminar / RANS / LES / Hybrid RANS-LES.
2. Base turbulence model — filtered by family and implementation status.
3. Wall treatment — resolved / wall function / future all-y+ families.
4. Corrections — rotation/curvature, production limiter, compressibility, roughness and model-specific options.
5. Transition — none; future Gamma-ReTheta / Gamma / other model-specific choices.
6. Applicability diagnostics — dimension, steady/transient, y+ target, mesh requirements, implementation status.

The selector should explain recommendations rather than silently changing the case. The chosen model remains explicit in the case and in validation provenance.

## Proposed long-term case representation

physics.turbulence
- family: RANS
- model: SST
- wall_treatment: resolved
- transition: none
- corrections: {...}
- initialization: {...}

This is preferable to encoding every combination as a new enum. The current flat model key remains supported during migration.

## Immediate PR scope

This PR implements the first architectural layer only:
- complete machine-readable turbulence catalogue;
- explicit implementation status;
- solver-ready filtering;
- model descriptions and caveats;
- transparent rule-based recommendations;
- Python schema/tests;
- audit documentation.

It does not claim that LES, DES, DDES, IDDES, transition, RSM or wall functions are implemented. Those remain separate implementation/V&V workstreams.

## Validation policy

A model can move from KERNEL_ONLY or PLANNED to SOLVER_READY only after:
1. complete governing-equation implementation;
2. boundary/initial-condition contract;
3. positivity and boundedness diagnostics where applicable;
4. model-specific unit/regression tests;
5. solver-level verification;
6. quantitative reference/benchmark validation;
7. mesh/y+ sensitivity evidence;
8. CI evidence on the exact final revision;
9. #118 qualification updates where applicable.

No catalogue status is promoted merely because a class, enum or algebraic helper exists.