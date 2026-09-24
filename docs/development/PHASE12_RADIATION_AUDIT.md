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

The four critical reservations from the external review are now closed in code and regression coverage.

### Rosseland opacity convention

The one-argument `rosseland_conductivity(T, transport_opacity)` API now explicitly defines its second argument as the transport/Rosseland opacity \(\kappa_R\), not absorption alone. For isotropic gray scattering:

\[
\kappa_R = \kappa_a + \kappa_s
\]

and the dedicated scattering-aware solve accepts an asymmetry factor and uses:

\[
\kappa_R = \kappa_a + \kappa_s(1-g)
\]

Thus the diffusion coefficient is not silently based on \(\kappa_a\) when scattering is present.

### Direction-aware diffuse-gray DOM walls

A dedicated DOM wall operator now evaluates the outgoing discrete irradiation on every boundary face:

\[
G_{out} = \sum_{m:\,\mathbf{s}_m\cdot\mathbf{n}>0}
w_m I_m(\mathbf{s}_m\cdot\mathbf{n})
\]

and applies the diffuse-gray incoming condition:

\[
I_{in}=\epsilon I_b(T_w)+(1-\epsilon)G_{out}/\pi
\]

The wall reflection is Picard-lagged by one transport iteration. Incoming ordinates receive fixed face values; outgoing ordinates retain zero-gradient/extrapolated transport. The generic finite-volume face-value adapter now supports this direction-specific override using non-finite sentinels for faces where no incoming condition is applicable.

This is now exercised by a non-black gray-wall equilibrium regression using a half-range-consistent 14-direction quadrature.

### S2S reciprocity

The area-aware view-factor validator checks both enclosure closure and:

\[
A_iF_{ij}=A_jF_{ji}
\]

A regression explicitly accepts a reciprocal unequal-area matrix and rejects a deliberately non-reciprocal matrix. The deterministic ray-traced kernel remains an estimator; it is not presented as an exact closed-form view-factor calculator.

### Angular sensitivity / ray effects

The implementation continues to expose the direction set through the DOM API, so angular refinement can be performed independently of the transport solver. The validation suite now uses a half-range-consistent quadrature for the diffuse-wall oracle, avoiding a false equilibrium caused by a full-sphere-only moment set. Ray-effect convergence remains a V&V campaign item rather than being claimed from a single low-order case.

## Acceptance interpretation

Implemented and regression-covered:

- area-aware surface exchange;
- enclosure and area-weighted view-factor reciprocity;
- deterministic S2S visibility estimation;
- direction-aware diffuse-gray DOM walls;
- DOM full-sphere moment validation;
- normalized DOM convergence and positivity guards;
- constant-property P1;
- Rosseland diffusion with explicit transport-opacity convention and scattering-aware wrapper;
- spatially varying optical properties;
- band-wise gray non-gray infrastructure;
- coupled radiation/energy verification;
- component and equilibrium regression oracles.

Explicit scope boundaries remain:

1. S2S is a deterministic sampling estimator; sample count controls its numerical error.
2. Spectral radiation is band-wise gray DOM, not line-by-line spectroscopy.
3. A broad S2/S4/S6/S8 mesh-and-angle convergence campaign is still a V&V expansion item; the current tests verify the mathematical operators and boundary treatment without claiming a universal angular-convergence result.

## External technical references used for the audit

The implementation choices were cross-checked against established radiation-model formulations from Ansys Fluent, OpenFOAM and NASA technical references. These references distinguish Rosseland, P1, S2S and DO models and document their respective equations, boundary conditions and validity assumptions.
