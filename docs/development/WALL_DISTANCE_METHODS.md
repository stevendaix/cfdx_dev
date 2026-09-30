# Wall-distance methods and benchmark

This document records the wall-distance methods implemented for the turbulence
architecture and the first distance-only benchmark.

## Methods

CFDX exposes one method enumeration so that turbulence models do not own
independent wall-distance implementations:

- EXACT_GEOMETRIC: point-to-triangle minimum; reference implementation.
- SEARCH_BASED: two-stage vertex estimate followed by exact face distance in
  a configurable near-wall threshold. This follows the design direction of
  Zhang et al. (NASA/FUN3D, ICCFD12 2024).
- MESH_WAVE: graph/advancing-front propagation from wall-near seeds.
- DIRECTIONAL_MESH_WAVE: directional variant of graph propagation.
- POISSON: solve -laplacian(phi)=1 with phi=0 on cut wall faces and outer
  homogeneous Neumann treatment, then use the Tucker reconstruction
  d = sqrt(|grad(phi)|^2 + 2 phi) - |grad(phi)|.
- EIKONAL: monotone upwind/fast-sweeping solution of |grad(d)| = 1 with a
  numerical wall seed band.
- HAMILTON_JACOBI: pseudo-time relaxation of the implemented modified
  equation |grad(d)| = 1 + epsilon*d*laplacian(d), with Godunov gradient.
- ADVECTION_DIFFUSION: transport form U.grad(d) = 1 + Gamma*laplacian(d),
  using upwind advection and central diffusion.
- HYBRID_POISSON_EIKONAL: Poisson initialization followed by H-J refinement;
  no arbitrary weighted blend of two distance fields.

The PDE methods in this first implementation deliberately live in the
distance benchmark layer. They are not yet wired into the production CFD
linear-system backend. This prevents the benchmark from silently inheriting
solver-specific preconditioning or convergence behaviour.

## First complex-geometry benchmark

`apps/verification/wall_distance_complex_benchmark.cpp` builds a triangulated
wing-body-tail-like multi-component structure and samples a 3-D Cartesian
volume around it. The reference is the exact point-to-triangle distance on
the same triangulated surface.

Reported quantities:

- global relative L2 error;
- relative L-infinity error;
- near-wall relative L2 error;
- monotonicity violations against the reference ordering;
- wall-distance setup/evaluation time.

The exact geometric method is a mandatory self-consistency gate. No method is
declared physically validated by this benchmark; it only measures the distance
field against a geometric reference. For Poisson, the separate manufactured
benchmark and error-budget decomposition are mandatory because a small linear
residual does not establish correctness of the reconstructed distance.

## Bibliography / external references

1. X. Zhang, B. Diskin, A. Walden, E. J. Nielsen, "A Wall-Distance Method
   for Turbulence Modeling", ICCFD12, NASA NTRS document 20240006379, 2024.
   https://ntrs.nasa.gov/citations/20240006379
2. OpenFOAM documentation, "Wall distance calculation methods",
   including mesh-wave, directional mesh-wave and Poisson.
   https://doc.openfoam.com/2306/tools/processing/numerics/schemes/wall-distance/
3. SU2 geometry implementation uses an accelerated wall-surface search
   (`CADTElemClass`) and stores the nearest wall element/distance.
   https://github.com/su2code/SU2
4. The production turbulence verification/validation campaign remains
   separate and is to be connected later through #118. This benchmark does
   not qualify SST/SA/DES/DDES/IDDES.