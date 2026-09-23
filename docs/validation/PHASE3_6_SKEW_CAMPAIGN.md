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

For every level it evaluates both the orthogonal and corrected Laplacian schemes and records:

- orthogonal operator response;
- corrected operator response;
- magnitude of the non-orthogonal correction;
- global internal-face conservation error.

## Acceptance gates

1. All reported values are finite.
2. The corrected operator is conservative across the two-cell internal face to 1e-12 in volume-integrated form.
3. Zero skew produces zero corrected-vs-orthogonal difference to 1e-12.
4. Every non-zero skew level produces a non-zero correction, proving that the corrected path does not silently collapse to the orthogonal scheme.
5. The existing analytical/skew regression in test_laplacian remains part of the suite.

## Scope boundary

This PR does not claim that the LIMITED non-orthogonal scheme is implemented; the existing explicit unsupported-path regression remains unchanged. It also does not replace later full-mesh MMS/refinement campaigns: those belong to the broader verification matrix and remain separate evidence where required.

The purpose of this PR is to turn the current single-skew regression into a reproducible quantitative skew sweep with explicit acceptance criteria.
