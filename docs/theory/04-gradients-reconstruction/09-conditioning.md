# 9. Conditioning and Rank Deficiency

For least squares, the local geometry is represented by a matrix such as

[
M=A^TWA.
]

Its spectral properties determine sensitivity to perturbations in the input data. For a symmetric positive-definite matrix, a condition estimate can be expressed as

[
kappa(M)=
rac{lambda_{max}(M)}
{lambda_{min}(M)}.
]

Large values indicate sensitivity; rank deficiency means that the stencil cannot uniquely determine all gradient components.

A robust implementation should distinguish well-conditioned, poorly conditioned, rank-deficient, insufficient and invalid stencils.

The response policy must be deterministic and documented. Silent fallback to another reconstruction changes the numerical scheme without exposing that change to the user or V&V system.
