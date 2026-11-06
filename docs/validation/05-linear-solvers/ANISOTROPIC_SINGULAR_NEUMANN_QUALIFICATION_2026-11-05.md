# CFDX: qualification of anisotropic singular-Neumann pressure solvers

Date: 2026-11-05. Scope: native C++17 CFDX linear algebra and a small
steady-incompressible production-path smoke test.

## Mathematical contract

For a finite-volume pressure operator assembled using symmetric face
transmissibilities, the pressure gauge is arbitrary:

    A 1 = 0, 1^T A = 0, 1^T b = 0.

The new projected-GMRES path uses P = I - (1 1^T)/n. The Krylov operator is
P A P, and every preconditioner input and correction is projected. Both
right and left null-space conditions are verified before setup. An
incompatible right-hand side or an incorrect null-space declaration returns
NOT_APPLICABLE without modifying the initial pressure vector.

A requested GMRES or FGMRES with constant-null-space projection uses the
flexible, right-preconditioned Arnoldi recurrence; this supports variable
AMG corrections. The final convergence certificate is always computed on
the original, unscaled physical matrix. Explicit matrix scaling together
with constant-null-space projection remains deliberately unsupported.

For homogeneous b = 0 and nonzero initial pressure, the relative residual
uses ||A x0|| as its reference. This prevents skipping certification or
reporting division by zero.

## Qualification problems

1. Two-point finite-volume diffusion on graded 2-D Cartesian meshes with
   Neumann walls, anisotropy kx/ky or ky/kx of 1e3 and 1e6, grids 64x64
   and 128x128, with an exact discrete pressure field. The test checks
   the true residual, gauge-invariant pressure error, conservative face
   flux balance, and initial arbitrary gauge.
2. A rotated (NE-direction) algebraically dominant coupling with an exact
   constant null mode, to avoid a purely axis-aligned qualification.
3. A continuous manufactured Neumann diffusion source with exact
   cos(pi*x)*cos(2*pi*y) solution on graded 12x12, 24x24, 48x48 meshes.
   The source's discrete compatibility correction is reported explicitly,
   not silently performed inside the solver.
4. An actual 2-cell steady-incompressible equilibrium case with pure
   Neumann pressure and both anisotropic AMG preconditioners via FGMRES.
5. Fail-closed scenarios: incompatible right-hand side, an anchored
   (nonsingular) pressure matrix mislabeled as singular, wrong left
   null-space, disallowed CG/preconditioner combination, and zero RHS.

## Results

- Release: 14/14 targeted CTest tests pass.
- ASan/UBSan: 8/8 groups in the new dedicated test pass; this is not a
  sanitizer run of the entire CFDX project.
- All 24 repeated runs (2 grid sizes x 2 directions x 2 anisotropy
  ratios x 3 repeats) for EACH multigrid combination pass the 1e-9
  original-matrix relative residual gate.
- The line preconditioner ALONE fails 12/12 runs at anisotropy 1e3
  (6 per grid size), reaching 350 iterations, but passes 12/12 runs at
  anisotropy 1e6. Failures are retained as benchmark observations.
- Continuous manufactured-source solution errors:

  | Grid | Relative pressure error |
  |---|---:|
  | 12x12 | 0.0258576335519 |
  | 24x24 | 0.00634203109973 |
  | 48x48 | 0.00157770238702 |

  Each factor-of-two grid refinement decreases the pressure error by
  approximately four, consistent with second-order spatial convergence
  on this *specific* graded-grid test.

- Rotated-coupling test: FGMRES+anisotropic Boomer 3 iterations and true
  relative residual 4.7835e-11, FGMRES+anisotropic smoothed aggregation
  3 iterations and 7.5421e-11.

128x128 wall-time medians (seconds) from 3 repeated runs, including PC
setup and Krylov solve, on the current machine:

| Direction | Ratio | CG + native AMG | FGMRES + anisotropic Boomer | FGMRES + anisotropic SA | FGMRES + lines |
|---|---:|---:|---:|---:|---:|
| x | 1e3 | 0.082733 | 0.131878 | 0.116837 | FAIL (350 iters) |
| y | 1e3 | 0.084921 | 0.095383 | 0.102827 | FAIL (350 iters) |
| x | 1e6 | 0.076398 | 0.071820 | 0.052721 | 0.022009 |
| y | 1e6 | 0.068016 | 0.068160 | 0.057608 | 0.019208 |

A small true residual is not the same as an accurate pressure field in a
very ill-conditioned matrix. The benchmark records both. These timings
are not reference performance measurements of official PETSc or HYPRE,
which are not installed, and are not end-to-end CFD wall times.

## Reproduction

Build the repository with CMake in Release mode. Then:

```bash
ctest --test-dir build -R 'test_(null_space|gmres_solver|fgmres_solver|amg_dispatch_certificate|krylov_preconditioning|krylov_diagnostics|amg_preconditioner_qualification|n8_amg_anisotropy_scaling|amg_chebyshev|anisotropic_krylov_matrix|petsc_hypre_native_port_anisotropy|neumann_anisotropic_pressure_qualification|amg_sparse_coarse_regression|n9_s10_pressure_reference_policies)' --output-on-failure

python3 tools/benchmarks/qualify_singular_anisotropic_pressure.py \
  --binary build/test_neumann_anisotropic_pressure_qualification \
  --sizes 64 128 --repetitions 3 \
  --raw raw_pressure_benchmarks.csv \
  --medians median_pressure_benchmarks.csv
```

The benchmark runner fails closed if any multigrid method fails the
requested residual gate. It marks the line-only failure as a measured
counterexample rather than a passing qualification.

## Important limits and next qualification gates

The 2-cell steady-incompressible test proves dispatch compatibility, not
physical qualification of NACA0012 at N128/N256. Large 3-D stretched,
nonorthogonal, skewed and unstructured meshes, variable coefficients,
coupled Navier-Stokes/SST, thermal/radiation, MPI and GPU are not yet
qualified by this change. Running the *complete* CFD solution with mass
conservation, traction/force accuracy, grid convergence, and total wall
time remains required.

This code is an independent native CFDX implementation of a documented
null-space-projection technique; no PETSc or HYPRE source files are copied.
PETSc/HYPRE remain algorithmic references, not installed dependencies.