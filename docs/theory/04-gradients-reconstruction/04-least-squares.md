# 4. Least-Squares Reconstruction

For cell `P` and neighbours `N_i`,

\[
\phi_{N_i}-\phi_P
\approx

\nabla\phi_P\cdot(\mathbf x_{N_i}-\mathbf x_P).
\]

Define

\[
\Delta\mathbf x_i=\mathbf x_{N_i}-\mathbf x_P,
\qquad
\Delta\phi_i=\phi_{N_i}-\phi_P.
\]

The gradient is obtained by minimising

\[
J(\mathbf g)=
\sum_i w_i
(\mathbf g\cdot\Delta\mathbf x_i-\Delta\phi_i)^2.
\]

With design matrix `A`, weights `W` and data `b`,

\[
J(\mathbf g)=|W^{1/2}(A\mathbf g-\mathbf b)|_2^2,
\]

leading formally to

\[
(A^TWA)\mathbf g=A^TW\mathbf b.
\]

The method is attractive for unstructured and polyhedral meshes because it uses local point geometry directly.

## Conditioning

The matrix `AᵀWA` must contain sufficient independent geometric information. In three dimensions, the displacement vectors must span the required space.

A production implementation therefore needs explicit diagnostics for insufficient neighbours, rank deficiency, near-singularity, excessive condition number and invalid weights. A successful linear solve is not evidence that the stencil is well conditioned.
