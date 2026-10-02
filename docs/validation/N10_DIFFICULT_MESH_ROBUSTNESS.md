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

## Invalid geometry

A deliberately collapsed face is passed to the same production mesh validator. The expected result is rejection with an explicit error. No repair or fallback is permitted.

## Relationship with existing campaigns

The existing non-orthogonal/skew Laplacian campaigns remain valid and are not replaced. N10 adds a cross-cutting mesh-quality ladder and solver robustness gate on top of those operator-specific tests.

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