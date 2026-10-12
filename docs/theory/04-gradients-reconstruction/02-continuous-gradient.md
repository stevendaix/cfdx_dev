# 2. The Continuous Gradient

For a sufficiently smooth scalar field `φ(x)` in `Ω ⊂ ℝᵈ`,

\[

\nabla\phi =
\begin{bmatrix}
\partial\phi/\partial x_1\\
\ddots\\
\partial\phi/\partial x_d
\end{bmatrix}.
\]

For a direction `d`, the directional derivative is

\[
D_{\mathbf d}\phi=
\nabla\phi\cdot\mathbf d.
\]

Thus the gradient contains the complete first-order spatial variation of a scalar field.

## Local linearisation

Around `x_P`,

\[
\phi(\mathbf x_P+\Delta\mathbf x)
=
\phi_P+
\nabla\phi_P\cdot\Delta\mathbf x
+\mathcal{O}(|\Delta\mathbf x|^2).
\]

This relation is the basis of cell-centred reconstruction: neighbouring values provide information from which the unknown gradient can be inferred.

## Why CFD needs gradients

Gradients enter viscous stresses, heat conduction, diffusion operators, pressure/velocity reconstruction, higher-order convection schemes, turbulence models and several multiphysics closures.

A gradient implementation must therefore be judged by its mathematical and geometric assumptions, not merely by whether it returns a vector.
