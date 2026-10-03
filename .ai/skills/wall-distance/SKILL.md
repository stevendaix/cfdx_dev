# Skill: Wall Distance

## Purpose
Guide implementation and verification of wall-distance fields used by turbulence and wall models.

## Scope
- Exact/geometric distance where available.
- PDE-based Poisson, Hamilton–Jacobi and hybrid approaches.
- Boundary conditions, corner behavior and reconstruction.
- Monotonicity, positivity and gradient quality.

## Method
1. Define the target distance field and boundary conditions mathematically.
2. Establish a trusted reference on simple geometries.
3. Verify the scalar distance value independently from its gradient/reconstruction.
4. Test positivity, monotonicity and wall-boundary behavior.
5. Use controlled geometry/refinement studies before complex meshes.
6. Report value error and gradient/reconstruction error separately.

## Required evidence
- Plane/channel analytical distance.
- Curved or spherical geometry reference.
- Refinement/convergence measurements.
- Monotonicity and positivity diagnostics.
- Corner/localisation diagnostics where relevant.

## Review rules
A small PDE residual is not sufficient evidence of accurate wall distance. Do not hide reconstruction-order or corner errors by changing the residual tolerance.
