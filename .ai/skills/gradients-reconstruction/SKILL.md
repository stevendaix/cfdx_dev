# Gradients and reconstruction

## Purpose

Guide implementation and audit of cell/vertex gradients and face reconstruction in CFDX.

## Required reasoning

1. Identify the reconstruction stencil and its support.
2. State whether the method is Green–Gauss, least-squares, weighted least-squares, or another formulation.
3. State geometric weights and normalization explicitly.
4. Detect rank deficiency and ill-conditioning in least-squares systems.
5. Define the boundary-neighbour policy separately from the interior stencil.
6. Separate gradient computation from face-value reconstruction and from limiting.
7. Verify formal order on appropriate mesh families rather than relying on a single geometry.

## Qualification evidence

Use constant and linear exactness tests first. Then use quadratic/polynomial or manufactured fields and controlled refinement on multiple mesh classes. Report the measured order, norm, stencil, boundary treatment, and geometry.

For polyhedral/irregular meshes, do not assume second-order accuracy merely because a least-squares or Green–Gauss formula is theoretically capable of it.

## Anti-patterns

- Do not conflate gradient order with face interpolation order.
- Do not conceal singular/near-singular systems with arbitrary regularisation.
- Do not accept a convergence slope without checking the asymptotic regime and refinement definition.
