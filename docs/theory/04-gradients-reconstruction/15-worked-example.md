# 15. Worked Example: A Complete Gradient Reconstruction Study

This section illustrates the complete workflow on a smooth manufactured field.

We use

\[
\phi(x,y)=x^2+2xy+3y^2
\]

with exact gradient

\[
\nabla\phi(x,y)=
\begin{bmatrix}
2x+2y\\
2x+6y
\end{bmatrix}.
\]

The field is deliberately simple: the exact derivative is known analytically, while the field is not linear. A quadratic field therefore tests more than linear exactness.

## What is measured?

At each cell centre \(P\), the numerical method returns \(\nabla\phi_P^h\). We compare it with the exact gradient at the same location.

A relative discrete \(L_2\) error is

\[
E_2 =
\left(
\frac{
\sum_P V_P
\left\|
\nabla\phi_P^h-\nabla\phi(\mathbf x_P)
\right\|_2^2
}{
\sum_P V_P
\left\|
\nabla\phi(\mathbf x_P)
\right\|_2^2
}
\right)^{1/2}.
\]

The volume weighting makes the norm represent the physical domain rather than merely the number of cells.

## Refinement

For mesh sizes \(h_i\) and errors \(E_i\),

\[
p_i =
\frac{\log(E_i/E_{i+1})}
{\log(h_i/h_{i+1})}.
\]

The observed order is an experimental result. It must not be replaced by the order expected from a formal derivation.

## What this experiment does not prove

This experiment alone does not qualify the CFDX gradient implementation. It does not prove correct boundary conditions for every physical model, robustness on every polyhedral topology, behaviour for discontinuous fields, nonlinear solver convergence or parallel equivalence.

Those properties require separate V&V evidence.
