# 04 — Gradients and reconstruction

**Status: IMPLEMENTED — repository-grounded theory and code map; N2 qualification remains a V&V question.**

This chapter is the scientific reference for cell-centred gradient operators in CFDX. It deliberately separates the mathematical definition of an operator from the evidence required to qualify it on the intended mesh families.

## 1. Why gradients matter

For a scalar field \\(\\phi\\), the gradient is

\\[
\\nabla\\phi =
\\begin{bmatrix}
\\partial\\phi/\\partial x\\\\
\\partial\\phi/\\partial y\\\\
\\partial\\phi/\\partial z
\\end{bmatrix}.
\\]

In a finite-volume solver it is not only a post-processing quantity. It enters diffusion, non-orthogonal corrections, face reconstruction, wall models, turbulence closures and higher-order convection schemes.

The central distinction is therefore:

> **A gradient implementation can be mathematically valid and still be unqualified on a given mesh family.**

OpenFOAM exposes the same conceptual separation between Gauss and least-squares gradients, interpolation, surface-normal gradients, divergence and Laplacian schemes. Its documentation also treats gradient limiting as a separate numerical operation. citeturn0search2turn0search1

## 2. Cell-centred finite-volume gradient

For a control volume \\(V_P\\), Gauss' theorem gives

\\[
\\int_{V_P} \\nabla\\phi\\,dV
=
\\oint_{\\partial V_P} \\phi\\,\\mathbf n\\,dS.
\\]

Approximating the volume-average gradient by face values gives

\\[
(\\nabla\\phi)_P
\\approx
\\frac{1}{V_P}
\\sum_{f\\in\\partial P}
\\phi_f\\mathbf S_{f,P},
\\]

where

\\[
\\mathbf S_{f,P}=\\mathbf n_{f,P}A_f.
\\]

The sign is cell-local: CFDX uses the stored owner orientation and reverses \\(\\mathbf S_f\\) when the current cell is the neighbour.

This is the mathematical reason that owner/neighbour orientation and face-area-vector consistency are part of the gradient contract.

## 3. CFDX cell-based Gauss gradient

The implementation in `src/cfdx/core/numerics/gradient.h` computes a face value from owner and neighbour cell values:

\\[
\\phi_f=(1-w)\\phi_P+w\\phi_N,
\\]

with

\\[
w=
\\frac{|\\mathbf C_f-\\mathbf C_P|}
{|\\mathbf C_N-\\mathbf C_P|}.
\\]

For a boundary face, the current Module-0 operator uses

\\[
\\phi_f=\\phi_P,
\\]

which corresponds to a zero-gradient treatment in this low-level operator. This must not be confused with the complete physical boundary-condition machinery.

The resulting operator is

\\[
(\\nabla\\phi)_P
=
\\frac{1}{V_P}
\\sum_f \\phi_f\\mathbf S_{f,P}.
\\]

### Numerical interpretation

The method is attractive because it is conservative in its surface-integral construction and inexpensive. Its accuracy, however, depends on the face-value reconstruction and on mesh geometry. On affine/appropriate meshes, a linear field can be reproduced exactly; skewness and non-orthogonality can introduce error.

## 4. Point-based / vertex Green–Gauss

CFDX also contains a vertex-based Green–Gauss path.

The implementation first constructs a vertex value by inverse-distance weighting:

\\[
\\phi_v =
\\frac{\\sum_{c\\in\\mathcal C(v)}
w_{cv}\\phi_c}
{\\sum_{c\\in\\mathcal C(v)}w_{cv}},
\\qquad
w_{cv}=\\frac{1}{|\\mathbf P_v-\\mathbf C_c|}.
\\]

It then averages vertex values to faces and applies the same Green–Gauss volume formula.

This is a **secondary reconstruction path**, not automatically equivalent to cell-based Gauss. The interpolation layer has its own consistency and accuracy requirements.

## 5. Least-squares gradient

For cell \\(P\\), CFDX solves

\\[
\\min_{\\mathbf g}
\\sum_{N\\in\\mathcal N(P)}
w_N
\\left[
\\phi_N-\\phi_P-
\\mathbf g\\cdot(\\mathbf x_N-\\mathbf x_P)
\\right]^2.
\\]

The normal equations are

\\[
\\mathbf A\\mathbf g=\\mathbf b,
\\]

with

\\[
\\mathbf A=
\\sum_N w_N\\mathbf d_N\\mathbf d_N^T,
\\qquad
\\mathbf b=
\\sum_N w_N\\mathbf d_N(\\phi_N-\\phi_P),
\\]

and

\\[
\\mathbf d_N=\\mathbf x_N-\\mathbf x_P.
\\]

CFDX implements the corresponding symmetric 3×3 system with partial pivoting.

The default low-level least-squares implementation uses inverse-distance-squared weighting:

\\[
w_N=|\\mathbf d_N|^{-2}.
\\]

OpenFOAM likewise provides least-squares gradients as a distinct second-order gradient scheme for suitable unstructured stencils. citeturn0search7turn0search9

## 6. Why rank matters

The matrix \\(\\mathbf A\\) must contain enough independent geometric directions to determine the requested gradient.

CFDX records:

- number of samples;
- requested dimension;
- numerical rank;
- minimum and maximum pivots;
- a pivot-based condition estimate;
- geometric scale;
- full-rank status;
- finite-data status.

A rank-deficient stencil cannot uniquely determine all components of a 3-D gradient.

For example, points lying on a plane cannot determine the normal derivative without additional information. A robust CFD code must detect this rather than silently interpreting an underdetermined system as a fully determined 3-D gradient.

## 7. Boundary stencil policy

CFDX explicitly defines:

~~~text
EXCLUDE_BOUNDARY
ZERO_GRADIENT_GHOST
REJECT_BOUNDARY_STENCIL
~~~

The policy is a numerical choice, not merely an implementation detail.

For `EXCLUDE_BOUNDARY`, boundary faces are not treated as neighbouring cells in the current least-squares stencil.

For `REJECT_BOUNDARY_STENCIL`, encountering a boundary face in the rejected configuration produces an explicit error.

A future fully coupled boundary reconstruction must inject the actual mathematical boundary condition rather than assuming zero gradient.

## 8. Extended and quadratic fits

A one-ring linear least-squares fit uses the model

\\[
\\phi(\\mathbf x)
\\approx
\\phi_P+\\mathbf g\\cdot\\mathbf d.
\\]

It is linear-exact, but curvature contributes to the truncation error on general unstructured meshes.

CFDX therefore contains extended-stencil and quadratic-basis paths. The quadratic model has the form

\\[
\\phi_N-\\phi_P
\\approx
g_xd_x+g_yd_y+g_zd_z
+m_{xx}d_x^2+m_{yy}d_y^2+m_{zz}d_z^2
+m_{xy}d_xd_y+m_{xz}d_xd_z+m_{yz}d_yd_z.
\\]

The fitted coefficient vector has nine components. The first three are the gradient.

This distinction is critical for N2: **linear exactness is not evidence of second-order convergence of the gradient on arbitrary polyhedral meshes.**

## 9. Weighted least squares

CFDX supports the weighting choices:

~~~text
UNIFORM
INVERSE_DISTANCE
INVERSE_DISTANCE_SQUARED
~~~

Changing \\(w_N\\) changes the metric of the reconstruction. It therefore changes sensitivity to near and far neighbours, conditioning and truncation error.

The choice must be evaluated experimentally rather than selected because it is conventional.

## 10. Scheme-selection API

The central dispatch distinguishes:

~~~text
GAUSS_TWO_POINT
GAUSS_POINT
LEAST_SQUARES
WEIGHTED_LEAST_SQUARES
LEAST_SQUARES_QUADRATIC
~~~

The physical consequence is that downstream operators can request a gradient strategy without duplicating its implementation.

This is particularly important for the Laplacian: the corrected, limited and over-relaxed forms consume a selected cell gradient rather than embedding one fixed reconstruction algorithm.

## 11. Gradient versus face reconstruction

A CFDX gradient is not itself a face-value scheme.

These are separate mathematical operations:

~~~text
cell values
   |
   +----> gradient reconstruction
   |          |
   |          v
   |       grad(phi)
   |
   +----> face-value reconstruction
              |
              v
            phi_f
              |
              v
        convection/diffusion flux
~~~

This separation is required to prevent the verification of one reconstruction layer from accidentally qualifying another.

## 12. Interaction with the Laplacian

For an internal face, CFDX decomposes the diffusive contribution into an orthogonal two-point term and, for corrected schemes, a non-orthogonal correction.

Let

\\[
\\mathbf d=\\mathbf C_N-\\mathbf C_P.
\\]

The orthogonal coefficient is

\\[
\\alpha=
\\frac{\\mathbf S_f\\cdot\\mathbf d}{|\\mathbf d|^2}.
\\]

The basic contribution is

\\[
F_f^{orth}
=
\\alpha(\\phi_N-\\phi_P).
\\]

The residual face vector is

\\[
\\mathbf S_f^{corr}
=
\\mathbf S_f-\\alpha\\mathbf d.
\\]

A reconstructed face gradient then gives

\\[
F_f^{corr}
=
\\mathbf S_f^{corr}
\\cdot
\\mathbf g_f.
\\]

Thus the accuracy of a corrected Laplacian is coupled to the quality of the gradient supplied to it.

## 13. Verification hierarchy for N2

The N2 campaign must test progressively stronger properties.

### Level 1 — algebraic/unit checks

Verify:

- constant field gives zero gradient;
- affine field gives the exact analytical gradient where the stencil is linear-exact;
- translation of the mesh does not change the gradient;
- scaling coordinates scales the gradient with inverse length;
- non-finite and degenerate geometry is rejected;
- rank-deficient stencils are diagnosed.

### Level 2 — geometric families

Use at least:

1. Cartesian structured cells;
2. affine-distorted cells;
3. skewed polyhedra;
4. tetrahedral families;
5. boundary-truncated stencils;
6. mixed/polyhedral meshes.

The same analytic field must be evaluated across all families.

### Level 3 — convergence

For a smooth manufactured field, define

\\[
E_h=
\\left(
\\frac{\\sum_P V_P
|\\mathbf g_P-\\mathbf g_{exact,P}|^2}
{\\sum_P V_P}
\\right)^{1/2}.
\\]

For two successive meshes,

\\[
p_{obs}
=
\\frac{\\log(E_h/E_{h/2})}
{\\log(2)}.
\\]

A claimed order must be supported by a mesh family for which the refinement ratio and geometric similarity are controlled.

### Level 4 — coupling

The gradient must also be tested as an input to:

- corrected Laplacian;
- diffusion flux;
- higher-order face reconstruction;
- turbulence/wall-distance operators where applicable.

Passing a standalone gradient test does not qualify these coupled uses.

## 14. Why N2 is not automatically qualified

A gradient can pass constant/linear exactness while showing a lower observed order on a particular tetrahedral family.

Therefore the following are separate claims:

~~~text
implemented
    ↓
mathematically consistent
    ↓
verified on controlled manufactured fields
    ↓
verified on target mesh families
    ↓
qualified for CFDX production scope
~~~

The current code provides substantially more than a single gradient implementation: multiple schemes, explicit boundary policy, rank/conditioning diagnostics and higher-order fitting paths exist. The remaining qualification question is whether the V&V campaign demonstrates the requested order and robustness on the declared mesh families.

## 15. File-to-equation map

| File | Scientific role |
|---|---|
| `core/numerics/gradient.h` | Gradient algorithms and dispatch |
| `core/numerics/gradient_stencil.h` | Boundary policy and stencil-quality contract |
| `core/fvm/least_squares_gradient.h` | LS/WLS algebra and conditioning diagnostics |
| `core/numerics/laplacian.h` | Consumption of gradients by non-orthogonal diffusion |
| `core/geometry/geometry_cache.h` | Cell centres, face centres, areas and volumes consumed by reconstruction |
| `core/field/field.h` | Field storage and metadata |
| `core/mesh/ownership.h` | Owner/neighbour orientation required for face signs |

## 16. Improvement roadmap

The scientific roadmap is:

1. make the face-value reconstruction layer independently selectable and testable;
2. complete boundary-aware reconstruction using explicit Dirichlet/Neumann information;
3. add robust rank/conditioning thresholds expressed in dimensionless geometric terms;
4. compare normal-equation LS against more stable QR/SVD formulations where conditioning requires it;
5. qualify vertex Green–Gauss only for declared mesh classes;
6. establish controlled second-order evidence for the quadratic fit;
7. add gradient limiting as a separate, independently verified layer;
8. connect each qualified gradient scheme to the flux/reconstruction schemes that consume it.

The objective is not to make every scheme pass every mesh. The objective is to make the scope, mathematical assumptions, diagnostics and evidence explicit.

## References

- OpenFOAM, *Numerical schemes*, including separate gradient, interpolation, Laplacian and divergence scheme definitions. citeturn0search2
- OpenFOAM, *Gradient schemes*, including Gauss, least-squares and limiting variants. citeturn0search1turn0search4
- Versteeg, H. K. & Malalasekera, W., *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*, 2nd ed., Pearson, 2007.
- Mavriplis, D. J., unstructured finite-volume reconstruction literature; use the project bibliography for the exact cited edition/paper before assigning a qualification requirement.
