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

The polyhedral Laplacian campaign additionally gates a genuine convergence property where one is actually required. With a linear-consistent least-squares cell gradient, the corrected non-orthogonal operator is linear-exact on interior tetrahedra, and the interior L2 error of the smooth manufactured field must strictly decrease under refinement; the observed order is reported but deliberately not thresholded, because the campaign gates convergence rather than a claimed polyhedral accuracy order. The default two-point Gauss gradient is excluded from this gate: it is not linear-consistent on tetrahedra, and its measured non-decreasing error is documented rather than hidden.

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
- solver convergence versus quality for the two-point diffusion matrix;
- invalid/degenerate rejection;
- owner-neighbour-aware internal quality metrics.

Remaining:

1. run the full smooth-field/polyhedral accuracy matrix over the N10 ladder;
2. connect the campaign to the imported production fixture manifest from #423;
3. retain quantitative evidence from the exact CI HEAD.