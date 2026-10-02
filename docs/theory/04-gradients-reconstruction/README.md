> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# Gradients and Reconstruction

This chapter is the first complete pilot for the CFDX scientific Theory format.

## Learning path

```{toctree}
:maxdepth: 1
:caption: Course

01-introduction
02-continuous-gradient
03-green-gauss
04-least-squares
05-weighted-least-squares
06-vertex-reconstruction
07-boundary-reconstruction
08-face-reconstruction
09-conditioning
10-gradient-limiting
11-accuracy-and-convergence
12-comparison-of-methods
13-cfdx-implementation
14-verification
15-worked-example
16-numerical-pitfalls
17-summary
```

## Executable science

The chapter includes an independent Python reference study for Green–Gauss, least-squares and weighted least-squares reconstruction on a deterministic skewed quadrilateral mesh.

```{toctree}
:maxdepth: 1
:caption: Executable experiments

experiments/README
```

## Scientific status

The reference study is a mathematical and teaching reference, not a CFDX qualification oracle. Formal N2 qualification remains governed by the V&V requirements and executable CFDX evidence.

## Traceability

Theory explains the numerical method. Developer documentation owns implementation contracts. V&V documentation owns evidence and qualification status.


## Repository-grounded scientific core

The implementation is in `src/cfdx/core/numerics/gradient.h`, `gradient_stencil.h` and `src/cfdx/core/fvm/least_squares_gradient.h`. The following equations are the contracts that the code implements.

### Cell-centred Green–Gauss

For a control volume (V_P),

$$
\int_{V_P}\nabla\phi\,dV
=
\oint_{\partial V_P}\phi\,\mathbf n\,dS,
$$

hence

$$
(\nabla\phi)_P
\approx
\frac{1}{V_P}\sum_f \phi_f\mathbf S_{f,P}.
$$

CFDX computes the face value for an internal face as

$$
\phi_f=(1-w)\phi_P+w\phi_N,
\qquad
w=\frac{|\mathbf C_f-\mathbf C_P|}
{|\mathbf C_N-\mathbf C_P|}.
$$

For the low-level Module-0 boundary operator, a boundary face uses (\phi_f=\phi_P), i.e. a zero-gradient treatment. This is not a replacement for the full physical boundary-condition machinery.

### Vertex Green–Gauss

The vertex path first reconstructs

$$
\phi_v=
\frac{\sum_c w_{cv}\phi_c}{\sum_c w_{cv}},
\qquad
w_{cv}=\frac1{|\mathbf P_v-\mathbf C_c|},
$$

then averages vertex values to faces and applies the same Green–Gauss integral. It is therefore a distinct reconstruction chain and requires independent verification.

### Least squares

CFDX minimises

$$
J(\mathbf g)=
\sum_{N\in\mathcal N(P)}
w_N
[\phi_N-\phi_P-\mathbf g\cdot(\mathbf x_N-\mathbf x_P)]^2.
$$

The normal equations are

$$
\mathbf A\mathbf g=\mathbf b,
$$

with

$$
\mathbf A=\sum_Nw_N\mathbf d_N\mathbf d_N^T,
\qquad
\mathbf b=\sum_Nw_N\mathbf d_N(\phi_N-\phi_P).
$$

The low-level implementation supports uniform, inverse-distance and inverse-distance-squared weighting and reports stencil rank and a pivot-based conditioning estimate.

### Rank and conditioning

A three-dimensional gradient requires three independent geometric directions. CFDX explicitly records sample count, rank, pivots, scale, condition estimate and full-rank status. A rank-deficient stencil is therefore diagnosed rather than silently promoted to a valid 3-D reconstruction.

### Boundary policies

The gradient stencil contract explicitly distinguishes:

```text
EXCLUDE_BOUNDARY
ZERO_GRADIENT_GHOST
REJECT_BOUNDARY_STENCIL
```

The policy is part of the numerical method and must be included in V&V.

### Higher-order reconstruction

The repository contains extended-stencil and quadratic-basis least-squares paths. The quadratic model contains the linear gradient terms plus the six independent quadratic terms:

$$
\phi_N-\phi_P\approx
g_xd_x+g_yd_y+g_zd_z+
m_{xx}d_x^2+m_{yy}d_y^2+m_{zz}d_z^2+
m_{xy}d_xd_y+m_{xz}d_xd_z+m_{yz}d_yd_z.
$$

The first three coefficients are the gradient. This is the relevant path for demonstrating genuine second-order behaviour on smooth unstructured families; linear exactness alone is insufficient evidence.

## Gradient versus face reconstruction

A gradient and a face-value reconstruction are separate operators:

```text
cell values
   ├──> gradient operator ──> grad(phi)
   │
   └──> face reconstruction ──> phi_f ──> numerical flux
```

The V&V of N2 must therefore avoid qualifying a face reconstruction merely because its gradient input passed a manufactured-solution test.

## Coupling to diffusion

CFDX's Laplacian implementation in `src/cfdx/core/numerics/laplacian.h` uses

$$
\alpha=
\frac{\mathbf S_f\cdot\mathbf d}{|\mathbf d|^2},
\qquad
\mathbf d=\mathbf C_N-\mathbf C_P,
$$

for the orthogonal contribution

$$
F_f^{orth}=\alpha(\phi_N-\phi_P).
$$

Corrected/limited/over-relaxed variants add a non-orthogonal contribution based on a reconstructed face gradient. Consequently, gradient accuracy and Laplacian accuracy cannot be treated as completely independent.

## N2 verification hierarchy

The minimum scientific ladder is:

1. **Algebraic properties:** constants, affine fields, translation invariance, coordinate scaling, finite/degenerate geometry checks.
2. **Stencil robustness:** rank-deficient, near-singular and boundary stencils.
3. **Mesh families:** Cartesian, affine-distorted, skewed, tetrahedral and general polyhedral families.
4. **Manufactured convergence:** weighted (L_2) gradient error and observed order
   $$
   p_{obs}=\frac{\log(E_h/E_{h/2})}{\log 2}.
   $$
5. **Coupled use:** corrected Laplacian, diffusion flux, face reconstruction and relevant turbulence/wall-distance consumers.

A passing unit test is evidence of implementation correctness for that property; it is not by itself qualification of N2.

## Code-to-equation map

| Source | Mathematical role |
|---|---|
| `core/numerics/gradient.h` | Gauss, vertex Gauss, LS/WLS and scheme dispatch |
| `core/numerics/gradient_stencil.h` | boundary policy and stencil-quality contract |
| `core/fvm/least_squares_gradient.h` | weighted normal equations, rank and conditioning |
| `core/numerics/laplacian.h` | gradient consumption in non-orthogonal diffusion |
| `core/geometry/geometry_cache.h` | cell/face geometry entering reconstruction |
| `core/mesh/ownership.h` | owner/neighbour orientation and face sign |
| `core/field/field.h` | cell-field storage and dimensions |

## Scientific references

OpenFOAM documents Gauss and least-squares gradients as distinct finite-volume schemes and separates gradient selection from interpolation and limiting. citeturn0search2turn0search1

The detailed numerical chapter should retain the project bibliography entry for Versteeg & Malalasekera (2007) and the exact benchmark papers used by the N2 V&V campaign.

## Status rule

N2 should only move from **PARTIAL** to **QUALIFIED** when the declared gradient schemes, boundary policies and mesh families have executable evidence covering accuracy, convergence, conditioning and coupled use. Implementation completeness and qualification remain separate statuses.
