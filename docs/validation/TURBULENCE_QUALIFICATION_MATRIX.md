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

## Equation-level audit status

The registry is intentionally not a claim that the current kernels reproduce the named literature models exactly. The first code audit found several items that must remain explicit qualification gaps:

- **standard k-epsilon:** the eddy-viscosity identity is present as nu_t = C_mu k^2/epsilon. Transport/source-term equivalence still requires term-level evidence.
- **RNG k-epsilon:** the RNG correction helper and reference constants are present. The complete production/destruction formulation and its coupling to the transport solver still require an independent equation audit.
- **Realizable k-epsilon:** a variable C_mu kernel exists, but the strain/rotation/invariant construction used by the generic eddy-viscosity path is not yet sufficient to claim full realizable-model equivalence. This is a qualification gap, not a PASS.
- **k-omega:** the baseline transport and closure paths now use the Wilcox-style nu_t = k/omega identity; the SST strain-rate limiter is kept only in the SST path. The remaining transport terms still require independent term-level verification.
- **SST:** blending, production limiting and cross-diffusion are implemented, but the complete F1/F2 dependency chain must be verified with independent reference values, including density/gradient conventions and wall-distance sensitivity.
- **Spalart-Allmaras:** the canonical fv1, fv2, ft2, cw1 and sixth-root fw structure is present. The implemented transport coupling and optional corrections still require term-by-term verification.
- **Smagorinsky:** the algebraic closure uses Delta^2 |S| with the configured Cs and remains kernel-level only.
- **WALE:** the closure now requires the independent tensor invariants S2 and Sd2 instead of fabricating Sd2 from the scalar strain rate. It remains kernel-level only until the full velocity-gradient path is wired.
- **Dynamic one-equation LES:** the generic closure now uses the one-equation SGS relation nu_t = Ck Delta sqrt(k_sgs), with the existing second scalar interpreted as k_sgs. The capability remains PLANNED because the SGS transport, test filtering and dynamic coefficient evaluation are not implemented.
- **DES/DDES/IDDES:** the closure dispatch now separates the DES length scale, DDES shielding-based length scale, and IDDES shielding/stress blend. They remain KERNEL_ONLY because base-RANS coupling, complete shielding inputs and transient hybrid-solver integration are still missing.

### Audit rule

A test that merely checks nu_t > 0, finite values, or successful convergence is a hardening/regression test. It is not evidence that the governing model equations are correct. Qualification tests must compare each independent term or closed-form identity against a frozen reference value or manufactured solution, with the reference convention recorded.

This audit is deliberately kept in #474 so implementation can progress without waiting for wall-distance, while preventing premature capability promotion.
