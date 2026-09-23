# Phase 12 / M4 Radiation — Deep Audit and Hardening

Date: 2026-09-23

## Audit scope

The audit covers the M4 radiation implementation currently on master: blackbody and gray-surface primitives; two-surface exchange; view-factor validation; P1 model utilities; DOM quadrature and participating-media transport; radiation/energy coupling; analytical, reference and coupled validation; and separation between implementation evidence and solver-level V&V.

CFDX keeps implementation, verification, validation, conservation and CI evidence as separate claims.

## Findings and changes in this PR

### 1. Surface exchange formula needed an explicit area convention

The original two-surface helper implicitly assumed equal surface areas. For a general two-surface enclosure, the resistance term depends on A1/A2.

This PR adds an area-aware overload:

q''_1 = sigma (T1^4-T2^4) / [(1-e1)/e1 + 1/F12 + (1-e2)/e2 * A1/A2]

The legacy overload remains available and explicitly maps to A1=A2=1.

### 2. View-factor validation was incomplete

The existing matrix check verified bounds and row closure, but not reciprocity. The new area-aware validator additionally checks A_i F_ij = A_j F_ji.

This is a necessary enclosure invariant. It is still a validator, not a geometric view-factor calculator.

### 3. Diffuse-gray wall emission/reflection was missing as an explicit primitive

The wall intensity relation I_w = e I_b + (1-e) G/(4 pi) is now exposed as a reusable primitive.

The current DOM solver still accepts explicit scalar intensity boundary values. A future production wall-BC adapter must update outgoing intensities from temperature, emissivity and incoming irradiation rather than hard-code a single intensity for all directions.

### 4. DOM quadrature validation was too weak

A weight sum of 4 pi is necessary but not sufficient. The validator now checks normalized directions, non-negative weights, zeroth moment sum(w)=4 pi, first moment sum(w s)=0, and second moments sum(w s_i s_j)=4 pi/3 delta_ij.

### 5. P1 was previously only a source/utility layer

The code contained P1 source and coefficient helpers but no dedicated scalar P1 solve path. This PR adds a constant-property gray isotropic P1 solve:

div(D grad G) - kappa_a G + 4 kappa_a sigma T^4 = 0

with D = 1/[3(kappa_a+kappa_s)] for isotropic scattering.

The path returns solver diagnostics and rejects non-physical negative irradiation.

### 6. DOM convergence used dimensional absolute tolerances

Intensity and radiation-source changes have different units, so comparing both directly against one absolute tolerance is not dimensionally meaningful. The outer DOM gate now uses normalized changes relative to the current field magnitude, while the inner linear solve retains its own residual tolerance.

### 7. Physical positivity is now checked

Negative intensity/irradiation is not accepted as a valid converged radiation state. Small round-off-level negative values are clipped to zero; larger negative values fail explicitly.

## Remaining gaps — intentionally not hidden

### A. General geometric view-factor computation

There is no general face-to-face visibility/ray-tracing/hemicube view-factor calculator in CFDX yet. The current code validates supplied view factors.

Required future acceptance: geometry-driven computation; occlusion/visibility handling; area reciprocity; enclosure closure; convergence/accuracy study; persistence/reuse after mesh changes.

### B. Production diffuse-gray DOM wall boundary conditions

The DOM solver currently receives ScalarBoundaryConditions containing prescribed intensity values. It does not yet implement a direction-aware diffuse-gray reflective wall operator.

Required: incoming/outgoing direction classification; hemispherical integration; emissivity-dependent reflection; temperature coupling; wall radiative heat-flux reporting.

### C. DOM angular refinement

The six-direction set is an invariant test, not an angular-convergence campaign. Required: at least 2–3 angular quadratures; a non-isothermal participating-medium case; QoIs versus angular order; conservation and energy-balance gates.

### D. Spatial refinement for radiation transport

There is no dedicated multi-cell DOM/P1 mesh-convergence campaign with observed spatial order. Required: at least three/four meshes; L1/L2/Linf or physically relevant QoIs; energy-balance error; observed order.

### E. Spatially varying optical properties

The participating solver currently takes scalar absorption/scattering controls for the complete domain. Production radiation needs cell/temperature/species-dependent properties, with explicit units and finite/positive validation.

### F. Spectral/non-gray radiation

M4 is currently gray. There is no wavelength-band model, Planck/Rosseland mean infrastructure or band-wise DOM/P1 coupling. This is a separate capability and should not be implied by the current M4 status.

### G. Full Rosseland energy integration

The Rosseland conductivity helper exists, but a dedicated nonlinear Rosseland energy solve and boundary treatment are not yet a closed acceptance path. Rosseland is a diffusion approximation for optically thick media and should be selected and validated only within its stated validity regime.

### H. Independent solver-level radiation reference

The current Level-B radiation check is an independent algebraic limiting-case oracle, not a full CFDX-vs-reference radiation solution. The next Level-B case should be a genuine multi-cell enclosure/participating-medium problem with an independently generated reference.

### I. Coupled radiation-energy validation is currently equilibrium-heavy

The existing Level-C case proves coupling infrastructure and an energy-balance gate, but the isothermal equilibrium setup does not exercise a non-trivial radiative heat transfer path. A stronger case should use non-uniform temperature or distinct radiating walls and report integrated radiative source, wall heat flux, total energy balance, temperature change, and mesh/angular sensitivity.

## Acceptance interpretation

After this PR:

- Implemented: blackbody/gray primitives, area-aware two-surface exchange, view-factor invariants, DOM moment validation, participating-media DOM core, a constant-property P1 solve, radiation-energy coupling.
- Verified at component/analytical level: Stefan-Boltzmann, gray exchange, view-factor identities, DOM moments, P1 equilibrium and source closure.
- Not yet validated as a complete radiation solver: geometric S2S, production diffuse-gray DOM walls, mesh/angular convergence, non-uniform participating-media reference case, full Rosseland path, spectral/non-gray models.

Therefore Phase 12 remains IMPLEMENTED / VALIDATION IN PROGRESS, not fully green.

## External technical references used for the audit

The implementation choices were cross-checked against established radiation-model formulations from Ansys Fluent, OpenFOAM and NASA technical references. These references distinguish Rosseland, P1, S2S and DO models and document their respective equations, boundary conditions and validity assumptions.
