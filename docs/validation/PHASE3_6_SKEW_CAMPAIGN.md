# Phase 3.6 — Non-orthogonal / skewed operator validation

This campaign closes the evidence gap identified in Issue #34 for Phase 3.6 at the **operator-regression level**.

## Coverage

The executable regression test_nonorthogonal_skew_campaign sweeps controlled face skew levels:

- 0.00
- 0.05
- 0.10
- 0.20
- 0.30
- 0.40
- 1.00
- 2.00
- 4.00
- 8.00

For every level it evaluates the orthogonal, corrected and `LIMITED(0.5)`
Laplacian schemes and records:

- orthogonal operator response;
- corrected operator response;
- limited operator response;
- magnitude of the non-orthogonal correction;
- magnitude of the limited correction;
- global internal-face conservation errors for corrected and limited fluxes.

## Acceptance gates

1. All reported values are finite.
2. The corrected and limited operators are conservative across the two-cell internal face to 1e-12 in volume-integrated form.
3. Zero skew produces zero corrected-vs-orthogonal difference to 1e-12.
4. Every non-zero skew level produces a non-zero correction, proving that the corrected path does not silently collapse to the orthogonal scheme.
5. The existing analytical/skew regression in test_laplacian remains part of the suite.
6. The limited correction never exceeds either the full correction or the
   orthogonal contribution. At severe skew 4 and 8, it is actively clipped.

## Scope boundary

`LIMITED` uses an OpenFOAM-style coefficient in `[0,1]`: zero disables the
correction, one applies the full corrected scheme, and the default `0.5` bounds
the non-orthogonal part by the orthogonal face contribution. This operator
campaign does not replace later full-mesh MMS/refinement campaigns: those
remain separate evidence where required.

The purpose of this PR is to turn the current single-skew regression into a reproducible quantitative skew sweep with explicit acceptance criteria.


## N10 reuse policy

N10 does not duplicate the existing spatial-accuracy campaigns.

The repository already contains executable quantitative V&V for the difficult-mesh dimensions relevant to N10:

- `test_gradient_verification`: smooth-field refinement on orthogonal, affine-skewed and high-aspect-ratio meshes, with L1/L2/Linf metrics and observed order;
- `test_polyhedral_gradient_campaign`: tetrahedral/polyhedral gradient accuracy, weighted least-squares and conditioning/rank-deficiency evidence;
- `test_nonorthogonal_skew_campaign`: controlled skew sweep;
- `test_nonorthogonal_laplacian_campaign`: non-orthogonal diffusion accuracy/regression;
- `test_polyhedral_laplacian_campaign`: polyhedral diffusion consistency limits.

These tests are now also labelled `n10;validation;qualification` so the N10 qualification subset reuses the same production numerical paths and does not maintain a second, redundant smooth-field campaign.

The dedicated N10 test `test_n10_difficult_mesh_robustness` is therefore restricted to the cross-cutting difficult-mesh contracts that are not already owned by N2/N3: mesh-quality diagnostics, solver robustness versus the controlled quality ladder, polyhedral geometry acceptance, near-degenerate valid geometry, and deterministic invalid-geometry rejection.
