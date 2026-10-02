# 16. Numerical Pitfalls

Gradient verification is particularly vulnerable to apparently convincing but incomplete tests.

## Polynomial exactness

Exact reproduction of a linear field is useful, but it is not sufficient evidence of second-order asymptotic convergence.

## Overfitting a mesh

A method can look excellent on an orthogonal structured mesh and deteriorate on skewed or highly non-orthogonal cells. The mesh family must reflect the intended application.

## Hiding conditioning failures

Replacing a failed stencil with another method may make a test pass while changing the numerical scheme. Any fallback must be explicit and documented.

## Measuring only residuals

Gradient accuracy is a field reconstruction question. A solver residual cannot substitute for an independent gradient error measurement.

## Mixing errors

If the reference, geometry, iterative tolerance and reconstruction method all change simultaneously, the observed difference cannot be attributed to one source.

## Accidental limiter activation

A limiter changes the algorithm. Formal-order experiments must record whether limiting is enabled and should normally test unlimited reconstruction separately.
