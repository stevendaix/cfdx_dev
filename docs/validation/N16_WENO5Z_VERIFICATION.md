# N16.1 — WENO5-Z finite-volume reconstruction

## Implemented

This PR implements the first concrete N16 method: an isolated WENO5-Z
finite-volume reconstruction kernel operating on five consecutive **cell
averages on a uniform 1-D stencil**.

It includes:

- standard FV WENO5 candidate polynomials;
- WENO-Z smoothness indicators and nonlinear weights;
- finite-value and parameter validation;
- an explicit bounded reconstruction policy;
- constant and linear exactness tests;
- a smooth sine-wave refinement campaign;
- a discontinuous-stencil boundedness test.

## Verification

For \(\phi(x)=\sin(2\pi x)\), exact cell averages are generated and the
reconstruction is compared with the exact value at the cell face.

The executable campaign uses 20, 40, 80 and 160 cells and requires every
observed refinement rate to exceed 4.0. This verifies the high-order
reconstruction independently from mesh geometry, convection flux assembly and
nonlinear solver convergence.

## Deliberate scope boundary

This PR does **not** claim:

- arbitrary 2-D/3-D WENO;
- tetrahedral/polyhedral WENO stencils;
- production convection registration;
- solver-level WENO convergence;
- general multidimensional monotonicity;
- dual-time or multirate integration;
- advanced AMG/Krylov/preconditioner work;
- matrix-free high-order operators.

The current kernel is therefore **experimental/verified as a 1-D uniform-grid
FV primitive**, not a production CFDX convection scheme.

## Next steps

1. Define a geometry-aware multidimensional stencil contract.
2. Build WENO face stencils on Cartesian meshes.
3. Add conservative flux assembly and flux-antisymmetry tests.
4. Extend to skew/non-orthogonal and polyhedral populations.
5. Compare accuracy, boundedness, robustness and cost with the N4 schemes before
   any production registration.
