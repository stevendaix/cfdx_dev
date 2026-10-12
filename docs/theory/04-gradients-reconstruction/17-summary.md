# 17. Summary

Gradient reconstruction connects the continuous differential operator to the discrete geometry used by a finite-volume solver.

The central ideas are:

1. The gradient is the first-order spatial variation of a field.
2. Green–Gauss reconstructs the gradient from a surface integral.
3. Least-squares methods infer the gradient from local point geometry.
4. Weighting changes the local approximation and must be specified explicitly.
5. Vertex methods introduce an additional reconstruction layer.
6. Boundary treatment is part of the numerical method.
7. Face reconstruction should remain separate from gradient reconstruction.
8. Conditioning is a first-class diagnostic for geometric stencils.
9. Linear exactness is not equivalent to demonstrated second-order convergence.
10. Verification requires controlled fields, mesh families, norms and reproducible evidence.

For CFDX, the objective is not simply to implement several gradient algorithms. Their mathematical definitions, implementation contracts and verification evidence must remain traceable.
