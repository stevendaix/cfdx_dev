# 3. Green–Gauss Reconstruction

The divergence theorem gives

\[
\int_V
\nabla\phi\,dV
=
\oint_{\partial V}\phi\,\mathbf n\,dS.
\]

Approximating the gradient as sufficiently smooth over a control volume gives

\[

\nabla\phi_P
\approx
\frac{1}{V_P}
\oint_{\partial V_P}\phi\,\mathbf n\,dS.
\]

For a polyhedral cell,

\[

\nabla\phi_P
\approx
\frac{1}{V_P}
\sum_f \phi_f\mathbf S_f,
\qquad
\mathbf S_f=\mathbf n_f A_f.
\]

The practical method depends critically on how `φ_f` is reconstructed.

## Numerical interpretation

Green–Gauss behaviour depends on:

1. face-value reconstruction;
2. mesh orthogonality and skewness;
3. boundary-face treatment;
4. geometric consistency;
5. any non-orthogonal correction.

Linear-field exactness under controlled conditions does not by itself demonstrate second-order convergence on arbitrary polyhedral meshes. That requires a refinement study with an independently defined reference field.
