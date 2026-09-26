# CFDX performance + efficient-physics program

This PR starts from the validated PR #9 M1-M4 baseline and keeps the CSR backend available for regression comparison.

## Implemented in this wave

### M5 — solver/runtime
- Matrix-free linear-operator abstraction with optional diagonal access.
- Matrix-free FVM diffusion and fused convection-diffusion operators.
- Reusable BiCGStab and contiguous GMRES workspaces.
- Adaptive GMRES restart and periodic true-residual replacement.
- Fused local Krylov reductions for future MPI all-reduce batching.
- Chebyshev polynomial smoothing suitable for matrix-free/GPU execution.
- Matrix-free agglomerated two-level V-cycle preconditioning.
- Solver-mesh representation separating construction topology from solver data.
- Explicit local-time-stepping and pseudo-transient controls.
- Inexact-Newton/Eisenstat-Walker forcing and adaptive pressure-correction controls.

### M7 — precision
The existing precision policy remains the control point for FP32 operators and FP64 reductions. Residual replacement is now part of GMRES control so reduced-precision paths can periodically restore an independently evaluated residual.

The implementation deliberately does not force FP16/BF16 globally: current CFD literature shows that selective precision is more robust than an unconditional precision downgrade.

### M8 — efficient convergence
- adaptive CFL;
- local time-step field;
- pseudo-transient CFL;
- physics-based linear tolerance;
- Eisenstat-Walker forcing;
- adaptive pressure correctors;
- reduction batching API.

The steady incompressible solver now consumes the inexact linear forcing and
adaptive pressure-corrector policies. They remain disabled by default so
existing cases are bit-for-bit policy compatible. Library callers enable them
through `IncompressibleSolverControls::acceleration`; the production executable
provides `--adaptive-convergence`. The first PISO/PIMPLE/fractional-step outer
iteration uses the configured maximum number of pressure corrections, then the
normalized continuity residual selects subsequent counts. Each iteration
records the effective linear tolerance and pressure-corrector count.

### M9 — efficient physics
- Boussinesq buoyancy;
- low-Mach thermodynamic density/preconditioning utilities;
- Spalart-Allmaras closure foundation;
- shared thermodynamic cache;
- existing ideal-gas state extended by cached reuse.

### M10 — radiation
- streaming DOM API with one active directional intensity;
- radiation workspace re-use in the participating-radiation driver;
- optical-thickness selector for Rosseland/P1/DOM;
- Rosseland diffusion coefficient and P1 coefficient helpers.

### Memory/runtime architecture
SolverMesh provides the compact representation intended to remain resident during iterations, while construction-only mesh data can be released after setup. The matrix-free path removes the fine-grid CSR matrix from the required solver representation. Krylov storage is contiguous and all large work vectors are allocated outside iteration loops.

## What remains intentionally backend-specific

The public interfaces for MPI halo exchange, GPU placement and block Navier-Stokes/Schur solvers are kept dependency-free. They can therefore be connected to the existing MPI/CUDA backends without imposing those dependencies on the core library.

The same applies to chemistry: reduced chemistry/tabulation is a later physics package because it changes the field working-set architecture substantially.

## Verification and benchmarking

test_performance_runtime exercises the new numerical/runtime controls. Existing PR #9 analytical, reference and coupled validation remains the correctness baseline.

A benchmark should compare at least:
- peak resident memory / cell;
- matrix storage bytes / cell;
- Krylov iterations;
- wall time / iteration;
- wall time / converged solve;
- residual history;
- MPI communication volume when a distributed backend is enabled.

No performance claim is considered validated until these quantities are measured on the same case, mesh, tolerance and hardware.

The convergence policy follows the same outer/inner accuracy principle described
in [SU2's linear solver guidance](https://su2code.github.io/docs_v7/Linear-Solvers-and-Preconditioners/)
and the residual-driven outer-corrector control exposed by
[OpenFOAM](https://doc.openfoam.com/2306/tools/processing/numerics/solvers/case-termination/).
