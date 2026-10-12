# N10 — Difficult-mesh numerical robustness

N10 closes the first executable difficult-mesh robustness layer of Issue #461.

The objective is not to declare every difficult mesh accurate. The objective is to make degradation measurable, deterministic and attributable to mesh quality.

## Quality ladder

The executable campaign `test_n10_difficult_mesh_robustness` contains:

| Level | Fixture | Main evidence |
|---|---|---|
| A | orthogonal | zero skew/non-orthogonality, aspect ratio = 1 |
| B | mildly skewed | non-zero skew/non-orthogonality and convergent linear solve |
| C | strongly skewed | increased geometric distortion and convergent linear solve |
| D | high aspect ratio | aspect ratio >= 20 and convergent linear solve |
| E | tetrahedral/polyhedral | topology/geometry validation and linear-gradient consistency |
| F | near-degenerate but valid | positive volume and finite geometry diagnostics |
| G | invalid | deterministic mesh rejection |

The campaign deliberately distinguishes valid difficult meshes from invalid meshes. A valid but poor-quality mesh must remain executable and observable; an invalid mesh must be rejected rather than repaired silently.

## Mesh-quality metrics

`MeshQualityReport` now exposes:

- maximum face skewness;
- maximum internal-face non-orthogonality;
- maximum cell bounding-box aspect ratio;
- minimum internal-face owner/neighbour volume ratio;
- minimum and maximum cell volume;
- surface-vector closure error.

For internal faces, skewness and non-orthogonality are evaluated from the owner-to-neighbour centre line. Boundary faces use the owner-to-face direction. This avoids classifying an internal face from the wrong geometric reference.

The aspect-ratio metric is a coordinate-aligned bounding-box metric. It is a diagnostic indicator, not a substitute for a full cell-shape-quality metric.

## Solver robustness

For the structured quality ladder, the test assembles the actual CFDX mesh topology and face geometry into a two-point finite-volume diffusion matrix:

$$a_{PN} = -\frac{|S_f|}{|C_N-C_P|}$$

for internal faces, with the corresponding diagonal contribution. Dirichlet boundary contributions use the owner-to-face distance.

The test constructs a manufactured algebraic right-hand side from a known solution and solves it with the native CG implementation. Acceptance uses an independently recomputed true residual:

$$\frac{\|b-Ax\|_2}{\|b\|_2}<10^{-8}.$$

This is a solver-robustness test, not a claim that the two-point diffusion discretisation is accurate on every difficult mesh.

## Polyhedral scope

The tetrahedral fixture verifies:

- valid positive cell volumes;
- finite skewness/non-orthogonality diagnostics;
- finite aspect-ratio diagnostics;
- exact reproduction of a linear field by the existing least-squares gradient on the tested mesh family.

It does not promote N10 to a universal polyhedral accuracy qualification. Smooth-field order on the complete quality ladder remains a subsequent V&V item.

The polyhedral Laplacian campaign enforces operator consistency at every campaign refinement level (n = 4/6/8/16): with a linear-consistent least-squares cell gradient, the corrected non-orthogonal operator must be linear-exact on interior tetrahedra at each level. For the smooth manufactured field under this consistent configuration, the interior residuals and observed orders are measured and reported without threshold: the corrected operator is not quadratic-exact and its gradient-based correction makes it non-symmetric, so its pointwise interior residual is O(1) and no convergence or accuracy order is claimed for it. A strict-decrease gate on these residuals was tried and measured to fail on the tested family; the campaign keeps the honest diagnostic instead. The default two-point Gauss gradient keeps its documented no-order status, and its measured non-decreasing error remains documented rather than hidden.

## Solution-level polyhedral Poisson convergence

The polyhedral Laplacian campaign above measures the operator interior
residual without a threshold. `test_polyhedral_poisson_convergence` closes
the solution-level criterion that those diagnostics cannot claim. On
Kuhn-tetrahedralised grids (n = 4/8/16 interior cubes) it solves the
manufactured Poisson problem Delta phi = -3 pi^2 phi with
phi = sin(pi x) sin(pi y) sin(pi z) and exact Dirichlet values on the
outer ghost shell.

The interior block of the discrete operator is never re-implemented: it is
assembled from black-box applications of the production `compute_laplacian`
operator (corrected scheme with the least-squares gradient) on colour-coded
indicator fields. A distance-5 cube colouring keeps same-colour interior
cells outside the two-hop dependency support of the corrected operator. A
deterministic pseudo-random field then gates the assembly itself: the
assembled matrix must reproduce `compute_laplacian(psi)` on every interior
cell to 1e-9 relative, so an incomplete dependency model or a wrong
colouring fails closed.

The system is solved with the production BiCGStab and accepted only with
CONVERGED status and an independently recomputed true relative residual
below 1e-8. Acceptance is a strict decrease of the volume-weighted L2
solution error at every refinement level n = 4 -> 8 -> 16. Observed orders
are reported as diagnostics only; no accuracy order is claimed for the
corrected polyhedral operator. If the strict-decrease gate fails, the test
fails: that is an honest negative result which blocks the N10 campaign,
not a tolerance to relax.


## Invalid geometry

A deliberately collapsed face is passed to the same production mesh validator. The expected result is rejection with an explicit error. No repair or fallback is permitted.

## Relationship with existing campaigns

N10 deliberately reuses the repository's existing numerical V&V rather than
creating a duplicate smooth-field accuracy suite.

The following production tests are part of the N10 qualification subset through
the `n10;validation;qualification` CTest label:

- `test_gradient_verification`: orthogonal, affine-skewed and high-aspect-ratio
  smooth-field refinement with L1/L2/Linf metrics and observed order;
- `test_polyhedral_gradient_campaign`: tetrahedral/polyhedral gradient accuracy,
  weighted least-squares and conditioning/rank-deficiency checks;
- `test_nonorthogonal_skew_campaign`: controlled skew sweep;
- `test_nonorthogonal_laplacian_campaign`: non-orthogonal diffusion V&V;
- `test_polyhedral_laplacian_campaign`: polyhedral diffusion consistency.
- `test_polyhedral_poisson_convergence`: solution-level polyhedral Poisson convergence under refinement (strict L2 decrease gate, observed orders reported as diagnostics only).

These cases already provide the spatial-accuracy evidence that a second
N10-specific smooth-field campaign would duplicate.

The dedicated N10 test therefore remains focused on the cross-cutting mesh
robustness contracts:

- quality metrics and controlled quality ladder;
- solver convergence/true residual versus mesh quality;
- polyhedral geometry acceptance;
- valid near-degenerate geometry;
- deterministic rejection of invalid geometry.

Production imported fixtures from #423 remain a separate integration item.

## Current N10 status

N10 remains **PARTIAL** after this PR.

Closed by this PR:

- controlled orthogonal/skew/non-orthogonal ladder;
- high-aspect-ratio fixture;
- valid near-degenerate fixture;
- polyhedral geometry fixture;
- solution-level polyhedral Poisson convergence gate on the Kuhn tetrahedral family;
- solver convergence versus quality for the two-point diffusion matrix;
- invalid/degenerate rejection;
- owner-neighbour-aware internal quality metrics.

Remaining:

1. run the full smooth-field/polyhedral accuracy matrix over the N10 ladder; the tetrahedral solution-level Poisson criterion is now executable and merge-gated through `test_polyhedral_poisson_convergence`, while the broader ladder (skew, high-aspect-ratio and imported-fixture families) remains open;
2. connect the campaign to the imported production fixture manifest from #423;
3. retain quantitative evidence from the exact CI HEAD.