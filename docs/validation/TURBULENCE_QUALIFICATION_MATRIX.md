# Turbulence qualification matrix

This document is the model-level qualification map for Issue #461 N17. It deliberately separates **implementation**, **verification**, and **validation**.

The canonical machine-readable registry is:
`docs/validation/TURBULENCE_QUALIFICATION_MATRIX.json`

The equation-level closure contract for all 13 models is:
`docs/validation/TURBULENCE_CLOSURE_EQUATIONS.md`

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
| LES | Dynamic one-equation | KERNEL_ONLY | No | Closure kernel and dynamic-coefficient algebra verified; SGS transport/test-filter/complete dynamic path remains incomplete |
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

- **standard k-epsilon:** the eddy-viscosity identity and transport/source-term reference terms are now covered by the equation-level test.
- **RNG k-epsilon:** the RNG C1* correction, production, destruction and diffusion coefficients now have an explicit independent reference test.
- **Realizable k-epsilon:** variable C_mu, invariant inputs, production, epsilon production/destruction and diffusion coefficients now have an explicit equation-level test; physical/mesh qualification remains open.
- **k-omega:** the baseline nu_t = k/omega identity, production, destruction, cross-diffusion sign convention and diffusion coefficients now have an explicit reference test.
- **SST:** blending, production limiting, cross-diffusion and diffusion coefficients have explicit reference terms; the F1/F2 dependency chain and wall-distance sensitivity still require dedicated verification.
- **Spalart-Allmaras:** fv1/fv2/ft2, modified vorticity, production, signed destruction, gradient source and effective diffusion are now covered by an explicit reference test.
- **Smagorinsky:** the algebraic closure uses Delta^2 |S| with the configured Cs and remains kernel-level only.
- **WALE:** the closure now requires the independent tensor invariants S2 and Sd2 instead of fabricating Sd2 from the scalar strain rate. It remains kernel-level only until the full velocity-gradient path is wired.
- **Dynamic one-equation LES:** the kernel relation nu_t = Ck Delta sqrt(k_sgs) and the dynamic coefficient algebra are now tested independently. The capability is KERNEL_ONLY: the closure kernel is implemented and verified, but SGS transport, test filtering and complete dynamic coefficient evaluation are not implemented.
- **DES/DDES/IDDES:** the closure dispatch and independent length-scale references now distinguish DES, DDES shielding and IDDES stress/blend behaviour. They remain KERNEL_ONLY because base-RANS coupling, complete shielding inputs and transient hybrid-solver integration are still missing.

### Audit rule

A test that merely checks nu_t > 0, finite values, or successful convergence is a hardening/regression test. It is not evidence that the governing model equations are correct. Qualification tests must compare each independent term or closed-form identity against a frozen reference value or manufactured solution, with the reference convention recorded.

## Missing work before the next promotion

- Dynamic one-equation LES: implement SGS transport, resolved/test filtering, dynamic coefficient evaluation from filtered fields, consistent coefficient clipping/backscatter policy, LES time integration and quantitative V&V.
- Smagorinsky/WALE: wire the closure kernels into the LES solver path with the actual velocity-gradient/filter-width pipeline, boundary/wall treatment, transient verification and quantitative V&V.
- DES: couple the length scale to a complete RANS base model and hybrid transient solver path, then add 3-D verification/V&V.
- DDES: implement the physical shielding input (including the model-specific r_d chain), couple it to the RANS base model and hybrid transient path, then add 3-D verification/V&V.
- IDDES: implement the complete model-specific shielding, stress/wake functions and blending, couple them to the RANS base model and hybrid transient path, then add 3-D verification/V&V.
- SST/SA: complete the independent wall-distance and wall-treatment qualification, then run model-specific quantitative campaigns.
- All RANS models: complete dedicated quantitative benchmark campaigns; the current SOLVER_READY status means the executable solver path exists, not that physical validation is complete.

This audit is deliberately kept in #474 so implementation can progress without waiting for wall-distance, while preventing premature capability promotion.
