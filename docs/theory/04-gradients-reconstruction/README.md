# Gradients and Reconstruction

This chapter develops spatial gradient reconstruction from continuous mathematics through finite-volume discretisation, polyhedral geometry, reconstruction algorithms and verification.

## Learning path

1. [Introduction](01-introduction.md)
2. [The continuous gradient](02-continuous-gradient.md)
3. [Green–Gauss](03-green-gauss.md)
4. [Least squares](04-least-squares.md)
5. [Weighted least squares](05-weighted-least-squares.md)
6. [Vertex reconstruction](06-vertex-reconstruction.md)
7. [Boundary reconstruction](07-boundary-reconstruction.md)
8. [Face-value reconstruction](08-face-reconstruction.md)
9. [Conditioning and rank deficiency](09-conditioning.md)
10. [Gradient limiting](10-gradient-limiting.md)
11. [Accuracy and convergence](11-accuracy-and-convergence.md)
12. [Comparison of methods](12-comparison-of-methods.md)
13. [CFDX implementation contract](13-cfdx-implementation.md)
14. [Verification](14-verification.md)
15. [Worked example](15-worked-example.md)
16. [Numerical pitfalls](16-numerical-pitfalls.md)
17. [Summary](17-summary.md)

## Executable science

The pilot now contains an independent Python reference study with:

- a deterministic skewed quadrilateral mesh;
- Green–Gauss, least-squares and weighted least-squares reconstruction;
- linear and quadratic manufactured fields;
- volume-weighted relative L2 errors;
- observed-order calculation over four refinement levels;
- conditioning diagnostics for LS/WLS;
- a deterministic self-check.

Run the study from `experiments/` with `python self_check.py` and
`python run_gradient_study.py`. The resulting numerical values are generated
at execution time; no convergence result is hard-coded into the documentation.

## Scientific status

The reference experiment is a mathematical cross-check and teaching instrument.
It is deliberately independent of the CFDX implementation and is therefore not
a qualification oracle. Formal N2 qualification remains governed by the V&V
requirements and executable CFDX evidence.

## Traceability

This chapter is the Theory-side explanation for the N2 gradient/reconstruction work.
Authoritative implementation and V&V records remain in their respective repository domains.
