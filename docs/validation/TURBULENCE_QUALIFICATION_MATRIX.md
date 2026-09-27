# Turbulence qualification matrix

This document is the model-level qualification map for Issue #461 N17. It deliberately separates **implementation**, **verification**, and **validation**.

The canonical machine-readable registry is:
`docs/validation/TURBULENCE_QUALIFICATION_MATRIX.json`

## Catalogue

| Family | Model | Implementation status | Wall distance | Current qualification boundary |
|---|---|---|---|---|
| Laminar | Laminar | SOLVER_READY | No | Full solver/physical validation can proceed independently |
| RANS | standard k-epsilon | SOLVER_READY | No | Equation/term/invariant verification now; physical campaign separate |
| RANS | RNG k-epsilon | SOLVER_READY | No | Equation/term/invariant verification now; physical campaign separate |
| RANS | Realizable k-epsilon | SOLVER_READY | No | Equation/term/invariant verification now; physical campaign separate |
| RANS | k-omega | SOLVER_READY | No | Equation/term/invariant verification now; physical campaign separate |
| RANS | k-omega SST | SOLVER_READY | Yes | Equation verification now; wall-dependent qualification remains separate |
| RANS | Spalart-Allmaras | SOLVER_READY | Yes | Equation verification now; wall-dependent qualification remains separate |
| LES | Smagorinsky | KERNEL_ONLY | No | Kernel equations/invariants only |
| LES | WALE | KERNEL_ONLY | No | Kernel equations/invariants only |
| LES | Dynamic one-equation | PLANNED | No | Contract/schema only until transport/filter implementation exists |
| Hybrid | DES | KERNEL_ONLY | Yes | Length-scale kernel only |
| Hybrid | DDES | KERNEL_ONLY | Yes | Shielding/length-scale kernel only |
| Hybrid | IDDES | KERNEL_ONLY | Yes | Shielding/blending kernel only |

## Common procedure

Each model is qualified through the same traceability chain:

1. equations and assumptions;
2. constants and literature source;
3. implementation audit;
4. term-by-term tests;
5. dimensional/invariant/positivity checks;
6. analytical or MMS verification where meaningful;
7. solver verification;
8. physical validation;
9. mesh/resolution and y+ sensitivity where applicable;
10. independent conservation and convergence diagnostics;
11. documentation;
12. capability-matrix promotion only after the applicable evidence exists.

A `SOLVER_READY` status means that a complete solver path exists. It is **not** a validation result.

## Wall-distance separation

Wall distance is treated as an explicit dependency for SST, SA, DES, DDES and IDDES. Its verification is tracked independently under the wall-distance workstream.

Therefore:

- wall-distance PASS does not imply turbulence-model PASS;
- turbulence equation verification does not imply wall-treatment validation;
- models without a wall-distance dependency must not be artificially blocked by #467;
- y+ sensitivity is a model/wall-treatment campaign, not a wall-distance implementation test.

## Promotion gates

A model is not promoted to production validation status from registry presence or kernel tests alone. Physical promotion requires a frozen reference, real CFDX execution, convergence, independent conservation checks, quantitative QoIs, appropriate refinement evidence and retained CI evidence.

No tolerance relaxation, disabled test, oracle substitution or fixed-iteration PASS is permitted.
