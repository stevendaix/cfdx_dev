# 12 — Comparison of gradient methods

This page explains how each method actually obtains a cell gradient, why the formula follows from the divergence theorem or a local Taylor expansion, what geometric information it needs, where its error enters, and what CFDX must verify.

## 1. Problem

For a scalar field \(\phi(\mathbf x)\), CFDX needs

$$
\mathbf g_P=(\nabla\phi)_P.
$$

The numerical problem is to reconstruct this quantity from cell values and mesh geometry.

The methods differ mainly in the mathematical identity they discretise.

## 2. Green–Gauss

The divergence theorem gives

$$
\int_{V_P}\nabla\phi\,dV
=
\oint_{\partial V_P}\phi\mathbf n\,dS.
$$

With the oriented face-area vector

$$
\mathbf S_f=\mathbf n_fA_f,
$$

we obtain

$$
\int_{V_P}\nabla\phi\,dV
=
\sum_f\phi_f\mathbf S_f.
$$

Approximating the cell-average gradient gives

$$
\boxed{
(\nabla\phi)_P
\approx
\frac1{V_P}\sum_f\phi_f\mathbf S_f
}.
$$

This immediately shows that Green–Gauss requires both geometry and face values.

For linear face interpolation,

$$
\phi_f=(1-w)\phi_P+w\phi_N,
\qquad
w=
\frac{|\mathbf C_f-\mathbf C_P|}
{|\mathbf C_N-\mathbf C_P|}.
$$

Therefore the actual chain is

$$
\phi_P,\phi_N
\rightarrow\phi_f
\rightarrow\phi_f\mathbf S_f
\rightarrow\nabla\phi_P.
$$

For a constant field \(\phi=C\),

$$
\nabla C=
\frac{C}{V_P}\sum_f\mathbf S_f=0,
$$

because a closed cell satisfies

$$
\sum_f\mathbf S_f=0.
$$

This is a direct geometric verification test.

For an affine field

$$
\phi=a+bx+cy+dz,
$$

the exact gradient is

$$
\nabla\phi=[b,c,d]^T.
$$

Affine exactness is a consistency property; it is not by itself proof of second-order convergence on arbitrary polyhedra.

## 3. Least Squares

Start from the local Taylor expansion

$$
\phi_N
=
\phi_P+
\nabla\phi_P\cdot
(\mathbf x_N-\mathbf x_P)+O(h^2).
$$

Define

$$
\mathbf d_N=\mathbf x_N-\mathbf x_P,
\qquad
\Delta\phi_N=\phi_N-\phi_P.
$$

Then

$$
\Delta\phi_N\approx\mathbf g_P\cdot\mathbf d_N.
$$

With several neighbours, choose \(\mathbf g_P\) by minimising

$$
J(\mathbf g)=
\sum_Nw_N
(\Delta\phi_N-\mathbf g\cdot\mathbf d_N)^2.
$$

Setting the derivative to zero gives the normal equations

$$
\boxed{A\mathbf g=\mathbf b}
$$

where

$$
A=\sum_Nw_N\mathbf d_N\mathbf d_N^T,
$$

and

$$
\mathbf b=\sum_Nw_N\mathbf d_N\Delta\phi_N.
$$

So LS constructs a local small linear system and solves for the gradient.

### Rank

In 3-D, three independent geometric directions are required. If

$$
\operatorname{rank}(A)<3,
$$

the gradient is not uniquely determined from the stencil.

This is why neighbour count alone is not a sufficient quality criterion.

## 4. Weighted least squares

WLS uses the same Taylor model but assigns different importance to neighbours:

$$
J(\mathbf g)=
\sum_Nw_N
(\Delta\phi_N-\mathbf g\cdot\mathbf d_N)^2.
$$

Examples are

$$
w_N=1,
\qquad
w_N=|\mathbf d_N|^{-1},
\qquad
w_N=|\mathbf d_N|^{-2}.
$$

The matrix remains

$$
A=\sum_Nw_N\mathbf d_N\mathbf d_N^T.
$$

The weight law therefore changes the local metric of the fit. It is part of the numerical method and must be fixed during a controlled comparison.

The formal solution is

$$
\mathbf g=A^{-1}\mathbf b
$$

when A is nonsingular.

Perturbations are amplified approximately according to

$$
\frac{\|\delta\mathbf g\|}{\|\mathbf g\|}
\lesssim
\kappa(A)
\frac{\|\delta\mathbf b\|}{\|\mathbf b\|}
+\cdots.
$$

Thus conditioning is not merely a diagnostic number: it measures sensitivity of the reconstructed gradient to perturbations.

## 5. Vertex Green–Gauss

The vertex method introduces an additional reconstruction stage.

For vertex \(v\),

$$
\phi_v=
\frac{\sum_{c\in\mathcal C(v)}w_{cv}\phi_c}
{\sum_{c\in\mathcal C(v)}w_{cv}},
$$

for example with

$$
w_{cv}=\frac1{|\mathbf x_v-\mathbf C_c|}.
$$

Vertex values are then used to construct face values, followed by Green–Gauss:

$$
(\nabla\phi)_P
\approx
\frac1{V_P}\sum_f\phi_f\mathbf S_f.
$$

The chain is therefore

$$
\phi_c\rightarrow\phi_v\rightarrow\phi_f\rightarrow\nabla\phi_P.
$$

Its additional error sources are vertex interpolation, vertex connectivity and boundary treatment.

## 6. What the methods actually depend on

| Method | Mathematical construction | Required information | Main error mechanisms |
|---|---|---|---|
| Green–Gauss | divergence theorem | face areas, normals, volume, face values | face reconstruction, geometry, boundary closure |
| Least Squares | Taylor expansion + minimisation | cell centres, neighbour values | stencil geometry, truncation, conditioning |
| Weighted LS | weighted Taylor fit | cell centres, neighbour values, weights | conditioning, weight law, stencil geometry |
| Vertex Green–Gauss | vertex reconstruction + divergence theorem | vertices, connectivity, face geometry | vertex interpolation, topology, geometry, boundary treatment |

This table is only the summary. The equations above define what each entry means.

## 7. Why second order cannot be inferred from the method name

If

$$
E_h\sim Ch^p,
$$

then

$$
p_{obs}=\frac{\log(E_h/E_{h/2})}{\log2}.
$$

Observed order depends on mesh family, refinement ratio, field regularity, boundary treatment, stencil construction, geometry, norm and iterative error.

Therefore linear exactness does not imply a measured second-order slope on an arbitrary polyhedral family.

## 8. Controlled comparison

Keep fixed:

$$
\boxed{
\text{field}+
\text{mesh family}+
\text{refinement ratio}+
\text{error norm}+
\text{boundary policy}+
\text{convergence criterion}
}.
$$

For an exact manufactured gradient, use for example

$$
E_h=
\left[
\frac{
\sum_PV_P
\|\mathbf g_P-\nabla\phi_e(\mathbf C_P)\|^2
}{
\sum_PV_P
\|\nabla\phi_e(\mathbf C_P)\|^2
}
\right]^{1/2}.
$$

Report L2, L-infinity, observed order, rank failures, conditioning statistics, boundary/interior error, mesh quality, CPU cost and the exact numerical-method configuration.

## 9. Bibliography

The Theory chapter should cite the mathematical sources behind the methods:

1. Versteeg, H.K. & Malalasekera, W., *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*, 2nd ed., Pearson, 2007 — finite-volume formulation, divergence theorem and discretisation.
2. Jasak, H., *Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows*, PhD thesis, Imperial College London, 1996 — finite-volume error and gradient analysis.
3. Mavriplis, D.J., work on unstructured-grid gradient reconstruction and least-squares methods — reconstruction and accuracy on unstructured meshes.
4. Barth, T.J. & Jespersen, D.C., 1989, *The Design and Application of Upwind Schemes on Unstructured Meshes*, AIAA Paper 89-0366 — unstructured reconstruction and limiting.
5. OpenFOAM documentation on gradient schemes — practical reference for separating gradient, interpolation and limiting.

The project bibliography should contain stable identifiers for each reference rather than relying only on prose.

## 10. CFDX mapping

| Mathematical step | CFDX implementation |
|---|---|
| face vector and cell volume | core/geometry |
| face interpolation | core/numerics/interpolation.h |
| Green–Gauss | core/numerics/gradient.h |
| LS/WLS normal equations | core/fvm/least_squares_gradient.h |
| stencil policy | core/numerics/gradient_stencil.h |
| diffusion coupling | core/numerics/laplacian.h |
| field storage | core/field/field.h |

Implementation mapping is not qualification evidence; the V&V chapter must link each method to executable tests and benchmarks.
