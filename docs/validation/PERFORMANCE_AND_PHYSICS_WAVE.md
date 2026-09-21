# CFDX performance and efficient-physics wave

This branch is based on PR #9, feature/verification-analytical-benchmarks, and preserves its M1-M4 verification campaign.

Implemented foundations:
- first-class matrix-free linear-operator abstraction;
- matrix-free orthogonal FVM diffusion operator;
- reusable Krylov and contiguous GMRES workspaces;
- mixed-precision policy and FP32/FP64 reduction helpers;
- lightweight matrix-free V-cycle preconditioner;
- low-storage SSPRK3/RK2 integration;
- adaptive CFL and pseudo-transient CFL;
- Boussinesq buoyancy;
- Spalart-Allmaras local closure;
- shared ideal-gas thermodynamic state;
- streaming DOM accumulation with one active intensity field.

The fine-grid matrix is no longer a mandatory representation. CSR remains supported as a compatibility backend, while matrix-free operators are intended for large FVM systems. The memory planner can now be used to reason about solver-vector lifetime, and the new workspaces remove repeated heap allocations from iterative methods.

The next production integration steps are direct use of matrix-free operators in scalar equation drivers, pressure AMG, residual replacement for mixed precision, MPI halo overlap, pipelined Krylov reductions, and full low-Mach/compressible fluxes.
