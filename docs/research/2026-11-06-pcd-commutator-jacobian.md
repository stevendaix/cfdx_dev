# CFDX — PCD du commutateur et certificat numerique du Jacobien

**Base** : sources du ZIP experimental PISO du 2026-11-06 (sans metadonnees Git).
**Statut** : implementation de recherche, opt-in. Ni qualification NACA, ni
validation des equations SST, ni equivalence SIMPLE/COUPLED n'est declaree.

## 1. Bibliographie verifiee

- D. Kay, D. Loghin, A. J. Wathen, *A Preconditioner for the Steady-State
  Navier-Stokes Equations*, SIAM J. Sci. Comput. 24(1):237-256 (2002),
  https://doi.org/10.1137/S106482759935808X.
- H. Elman et al., *Effects of Boundary Conditions on Preconditioning
  Strategies for the Navier-Stokes Equations*, transparency 3,
  https://www.maths.dur.ac.uk/lms/088/Talks/elman.pdf.
- H. Elman, V. Howle, J. Shadid, R. Shuttleworth, R. Tuminaro,
  *Block Preconditioners Based on Approximate Commutators*,
  SIAM J. Sci. Comput. 27(5):1651-1668 (2006),
  https://www.cs.umd.edu/users/elman/papers/spac.pdf.
- OpenFOAM 14, `incompressibleFluid/correctPressure.C`: consistent
  pressure-response and non-orthogonal correctors, available at
  https://github.com/OpenFOAM/OpenFOAM-14/blob/master/applications/modules/incompressibleFluid/correctPressure.C.

The FE derivation is

    (B F^-1 B^T)^-1 ~= Mp^-1 Fp Kp^-1.

For CFDX S = C - D Auu^-1 G and C = 0 with compatible block signs,
we use -Mp^-1 Fp Kp^-1. This sign must be verified against each real
assembled operator, not applied by analogy to a different convention.

**WARNING**: PCD needs pressure-space boundary operators consistent with
velocity BCs and with the discrete commutator; a valid algebraic composition
is not automatically a robust collocated finite-volume preconditioner.

## 2. Implemented opt-in

`src/cfdx/core/linalg/pcd_schur.h` now supports:

- `LegacySingleConvectionDiffusionInverse`: historical `-Fp^-1 rhs`, default;
- `LumpedMassCommutator`: `z=Kp^-1 rhs`, `y=Fp*z`,
  `pressure=-diag(Mp)^-1*y`.

The second mode rejects non-diagonal mass, nonpositive mass entries, a stale
or inaccurate Kp callback, nonfinite arithmetic and graph changes. `apply`
calculates into a temporary so failure cannot damage the caller's output.
`update_pressure_values_and_solves(...)` can atomically replace the pressure
matrices and their solve callbacks. An update to the matrix used by a nested
inverse without a new matching callback is rejected in BOTH modes.

Enable from a C++ caller explicitly:

```cpp
controls.algorithm = cfdx::physics::PressureVelocityAlgorithm::COUPLED;
controls.coupling.schur_model = cfdx::physics::CoupledSchurModel::PCD;
controls.coupling.pcd_lumped_mass_commutator_experimental = true;
```

The flag is rejected with non-COUPLED or non-PCD algorithms. There is no
silent fallback, no automatic mode selection, and no .cfdx.h5/CLI exposure
in this increment. Default controls and existing physical gates stay
unchanged. The experimental Kp true residual uses a 1e-9 relative check,
stronger than the existing 1e-8 default nested Kp solver request; physical
cases could reject until nested accuracy and gauge behaviour are qualified.

## 3. Added operator Jacobian certificate

`src/cfdx/physics/pressure_velocity_jacobian_certificate.h` provides
`certify_pressure_velocity_jacobian(assembled, state, residual, settings)`.
It compares columns of the assembled matrix [M G; D C] to **CENTRAL finite
differences of an independently supplied residual callback** at the same
state, including boundary and gauge equations. It returns the number of
columns tested, the maximum absolute discrepancy, maximum tolerance ratio
and worst (row, column). Input/nonfinite residual failures are fail-closed.
All columns are tested by default; deterministic sampling is an explicit
opt-in and never represented as exhaustive.

Tests create an independent affine residual, then inject a missing pressure
coupling C entry and a wrong-sign G entry. These are mathematical regression
oracles, **not yet finite-difference probes of production Rhie-Chow code**.
The next necessary step is to wire the true production residual callback
at the exact assembly state, with the same boundary conditions and gauge.

## 4. Verification

Execute from the project root:

```sh
cmake -S . -B build-pv -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON
cmake --build build-pv --target test_pcd_schur \
  test_pressure_velocity_jacobian_certificate \
  test_pressure_velocity_system test_steady_incompressible_solver \
  test_ghia_cavity
ctest --test-dir build-pv --output-on-failure \
  -R '^(test_pcd_schur|test_pressure_velocity_jacobian_certificate|test_pressure_velocity_system|test_steady_incompressible_solver|test_ghia_cavity_quick)$'
```

The unit suite covers 17 PCD cases, 5 Jacobian checks and 32 steady solver
cases. CTest targeted suites additionally exercise the quick Ghia Re=100
cavity case. Refer to the accompanying exact local test log for status.

## 5. Open qualification gates

1. Real FV Jacobian finite-difference comparison, including Rhie-Chow,
   skew faces, gauge and all pressure boundary types.
2. Nonzero-convection exact-Schur action comparison on matrices CAPTURED
   from the real COUPLED solver, old and new PCD with same solver budget.
3. Ghia Re=100/1000, Poiseuille, 3-D channel and NACA N64/N128/N256,
   with exact same mesh, gradients, BC, true nonlinear/linear residuals,
   mass, pressure and viscous-force gates; no tolerance relaxation.
4. Separate SST/energy tests to verify actual returned nonlinear state.
5. Expose the option in the user-facing configuration only after physical
   qualification. MPI/AMG parallel tests are not covered here.

## 6. Why no claim of speedup

No representative converged physical case has yet established that
`LumpedMassCommutator` outperforms the legacy preconditioner. The historical
512-unknown production matrix actually **rejected** the composite PCD mode,
so this should be investigated as a pressure-operator scaling/BC problem,
not hidden by loosening a threshold. The oracle and tests added here form
an auditable starting point, not evidence of physical parity or performance.
