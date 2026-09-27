# Turbulence correction capability matrix

This document records the distinction between correction kernels that exist in CFDX and model-specific production capabilities that are actually selectable.

The distinction is required by Issue #461 N17: a helper function implementing a mathematical factor is not evidence that the correction is correctly coupled to a turbulence transport model.

## Current production capability

| Base model | Rotation/curvature | Compressibility | Roughness | Production limiter | Kato-Launder | QCR |
|---|---|---|---|---|---|---|
| k-epsilon | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| RNG k-epsilon | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| Realizable k-epsilon | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| k-omega | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| SST | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| Spalart-Allmaras | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| Smagorinsky | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| WALE | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| DES/DDES/IDDES | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |

This is intentionally conservative. The generic correction kernels remain available for unit-level testing, but the capability registry does not advertise them as production model options.

## Promotion gate

A correction can be changed from unavailable to supported only when all of the following exist:

1. model-specific mathematical formulation;
2. explicit coupling point in the production transport/closure equations;
3. coefficient/configuration contract;
4. equation-level regression against the documented formulation;
5. positivity/boundedness checks where applicable;
6. at least one solver-level verification case;
7. validation evidence appropriate to the correction;
8. capability-registry metadata and compatibility tests.

A correction must not become selectable merely because its generic helper produces a finite value.

## Relationship to wall distance

Wall distance is an independent dependency. Roughness, SA variants and hybrid shielding may use wall distance, but a correct wall-distance kernel does not by itself validate the associated turbulence correction. The wall-distance method and the turbulence-model coupling therefore retain separate V&V evidence.

## Next V&V wave

The next turbulence work should qualify the six SOLVER_READY RANS models at equation level:

- standard k-epsilon;
- RNG k-epsilon;
- Realizable k-epsilon;
- k-omega;
- k-omega SST;
- Spalart-Allmaras.

For each model, retain separate evidence for transport equations, source terms, limiting, wall treatment, mesh sensitivity and canonical-flow validation.
