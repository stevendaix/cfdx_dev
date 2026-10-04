# 5. Weighted Least Squares

Weighted least squares uses

\[
J(\mathbf g)=
\sum_i w_i
(\mathbf g\cdot\Delta\mathbf x_i-\Delta\phi_i)^2.
\]

A common family uses distance-dependent weights such as

\[
w_i\propto|\Delta\mathbf x_i|^{-p}.
\]

The precise weighting law is part of the numerical-method definition and must not be treated as an undocumented implementation detail.

## Why weighting matters

Irregular meshes may contain large variations in neighbour distance and angular distribution. Weighting can reduce the influence of distant observations, but it cannot remove geometric degeneracy.

CFDX should therefore expose the weighting policy and report conditioning diagnostics independently from the final gradient.
