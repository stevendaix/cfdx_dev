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

## Audit closure

All actionable implementation gaps identified in the original audit are now represented in the PR:

- geometry-driven deterministic ray-traced S2S visibility;
- direction-aware diffuse-gray DOM wall operator;
- spatially varying optical properties;
- band-wise non-gray DOM;
- nonlinear Rosseland energy integration;
- radiation balance diagnostics;
- non-isothermal coupled radiation/energy verification;
- controlled angular sampling/refinement through the DOM direction API.

Two scope boundaries remain explicit rather than hidden:

1. the S2S geometry kernel is a deterministic Monte-Carlo estimator, so its numerical accuracy is controlled by sample count rather than an exact closed-form solution;
2. spectral radiation is band-wise gray transport, not continuous line-by-line spectroscopy.

The implementation therefore covers the Phase-12 M4 model family without claiming capabilities outside the phase.

## Acceptance interpretation

After this PR:

- Implemented: blackbody/gray primitives, area-aware two-surface exchange, view-factor invariants, DOM moment validation, participating-media DOM core, a constant-property P1 solve, radiation-energy coupling.
- Verified at component/analytical level: Stefan-Boltzmann, gray exchange, view-factor identities, DOM moments, P1 equilibrium and source closure.
- Not yet validated as a complete radiation solver: geometric S2S, production diffuse-gray DOM walls, mesh/angular convergence, non-uniform participating-media reference case, full Rosseland path, spectral/non-gray models.

Therefore Phase 12 remains IMPLEMENTED / VALIDATION IN PROGRESS, not fully green.

## External technical references used for the audit

The implementation choices were cross-checked against established radiation-model formulations from Ansys Fluent, OpenFOAM and NASA technical references. These references distinguish Rosseland, P1, S2S and DO models and document their respective equations, boundary conditions and validity assumptions.
