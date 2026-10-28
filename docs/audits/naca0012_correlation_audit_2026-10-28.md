# NACA0012 Re=1000, alpha=0: correlation audit and implementation gates

Date: 2026-10-28

## Scope and provenance

This audit summarizes the experiment report generated from the user-provided
`cfdx_dev-fix-f03-rhie-chow-grid-qualification.zip`, identified as a snapshot
of the code around PR #1133. It is **not evidence for the current HEAD** and
its experimental patch must not be applied blindly to the current branch.
The values below are diagnostic observations, not qualification results.

## Evidence from the snapshot

- A symmetrized O-grid with 24 radial layers converged to Cd = 0.0901496.
- A 36-layer variant converged to Cd = 0.120912, but the 48- and 72-layer
  variants diverged; reduced relaxation only delayed stagnation.
- Pressure drag evaluated using owner-cell pressure has a demonstrated
  affine-field face-integration error (6.45%, 3.10%, 1.84%, 0.89% for
  24/36/48/72 layers in the reported manufactured test).
- Wall-viscous drag reconstructed from owner-cell Green-Gauss gradients differs
  materially from the diagnostic 1-D normal estimate. Neither estimate is a
  certified oracle on the reported highly oblique faces.
- The report describes a finite trailing-edge thickness of 0.00252c, while
  some cited reference configurations use a sharp trailing edge.
- The report did not isolate the farfield boundary-condition effect or certify
  second-order convection.

The numerical proximity of one grid to Cd ~= 0.12 is **not** validation:
there is no demonstrated grid convergence and the force post-processing is
not yet independently qualified.

## Work items and acceptance evidence

### M1 — Geometry and mesh quality

- Keep the outer-boundary reflection-symmetry regression in PR #1141.
- Add reproducible mesh-quality reporting: radial first-layer spacing, growth
  ratio, cell area positivity, edge-length aspect ratio and internal-face
  non-orthogonality. The reported spacing is currently Euclidean displacement,
  not a measured wall-normal distance; skewness and a true wall-normal
  distance remain to be implemented.
- Replace the current linear blend with an appropriately designed C/O-grid
  normal-layer construction and smoothing only after a positive-Jacobian
  regression and mesh-quality baseline exist.
- Document the exact trailing-edge geometry. A sharp-TE option requires a
  topology that does not create a degenerate O-grid cell.

**Gate:** all cells valid, reflection symmetry retained, and reported near-wall
resolution/quality for each refinement level. Do not call the mesh qualified
from topology validity alone.

### M2 — Nonlinear solver robustness

- Reproduce the 48/72-layer divergence on the current HEAD before attributing
  it to a solver defect.
- Test line-search/backtracking with complete rollback of fields and face
  fluxes; log physical momentum and continuity residuals before divergence.
- Check consistency of the COUPLED velocity update, H/A predictor and
  Rhie-Chow face-flux reconstruction.

**Gate:** repeatable convergence on all selected grids without relaxing
physical residual criteria; retries must restore the full state.

### M3 — Face-based force post-processing

- Integrate pressure using a face pressure reconstruction with documented
  boundary treatment, not the owner-cell value by default.
- Add manufactured affine-pressure tests, including a curved-wall case, and
  preserve the constant-pressure closed-surface zero-force test.
- Reconstruct velocity gradients/viscous traction at faces with a documented
  order and compatible no-slip treatment.
- Report Cp, Cf, Cd_pressure, Cd_viscous, Cl, Cm, and a control-volume momentum
  balance independently.

**Gate:** manufactured tests pass and force estimates demonstrate spatial
convergence; do not select a reconstruction because it produces Cd ~= 0.12.

### M4 — Physical correlation and qualification

- Document reference geometry, Reynolds number, angle, farfield size,
  boundary conditions, steady/unsteady interpretation, and coefficient
  convention before selecting an external reference.
- Compare boundary-condition and convection-scheme variants only on
  converged meshes.
- Require at least three converged grids and estimate continuum error before
  comparing Cd with a literature value.

**Gate:** all grids converge; physical residuals and conservation gates pass;
Cd_pressure and Cd_viscous are independently stable; the comparison accounts
for geometry/reference differences and uncertainty.


## Local verification evidence — 2026-10-28

The repository snapshot was configured in Release mode with CMake 3.31.6,
GCC 14.2 and HDF5 1.14.5. MPI was not available in this environment, so these
results are serial.

- The mesh symmetry/quality Python tests pass: **13/13**.
- The mesh generator emits parseable `naca0012_mesh_quality.json`; the report
  serialization newline defect is fixed.
- The C++ NACA0012 quick topology/BC contract passes.
- The selected thermal and turbulence C++ tests pass: **7/7**, covering
  `test_thermal_vv`, `test_thermal_zones_vv`,
  `test_phase10_turbulence_hardening`,
  `test_turbulence_qualification_equations`,
  `test_production_turbulence_contract`, the mesh-generation fixture and the
  NACA0012 quick contract.
- The turbulence qualification-matrix and qualification-registry scripts pass.

These are software verification results, not a completed physical campaign.
The current mesh report exposes substantial quality limitations:

| Surface/radial levels | Mean first-layer spacing (chord units) | Maximum edge-length ratio | 95th percentile non-orthogonality | Maximum non-orthogonality |
|---|---:|---:|---:|---:|
| 64 / 24 | 6.47e-2 | 102.43 | 20.09 deg | 53.57 deg |
| 128 / 48 | 1.86e-2 | 117.63 | 21.43 deg | 65.96 deg |
| 256 / 96 | 5.33e-3 | 157.06 | 21.93 deg | 73.24 deg |

All reported 2-D polygon areas are positive, but that alone does not qualify
the mesh. The first-layer spacing is not yet wall-normal, the maximum
non-orthogonality increases with refinement, and the maximum edge-length ratio
is large. The generator still needs a proper boundary-layer/O-grid construction
with orthogonality and skewness control; merely changing the radial stretching
law is not an acceptable fix.

A local single-grid solve on the 64/24 mesh converged in 456 iterations with
`Cd = 0.0901496`, `Cl = -9.15e-11` and normalized continuity
`8.84e-11`. This is reproducible diagnostic evidence, not correlation: the
drag remains below the cited approximate `Cd ~= 0.12`, and one grid does not
establish convergence. The three-grid `test_naca0012_qualification` campaign
did not complete within a 300-second local limit; therefore no three-grid
qualification result is claimed.

An exploratory exponential radial-stretching change was rejected and is not
in the implementation: on the 64/24 grid it caused pressure-correction CG to
return `NOT_APPLICABLE` at outer iteration 13 (relative residual about 0.673).
That negative result is why no radial-stretching default was changed to make
the first-layer spacing look better while destabilizing the actual solve.

The existing thermal/turbulence tests establish component-level contracts only.
The registry still marks `TUR-NACA0012`, `TUR-FLATPLATE` and
`TUR-CHANNEL` as planned. A production RANS NACA0012 run, wall-`y+`
evidence, independent force reconstruction, and thermal/CHT/radiation coupling
in a real external-flow qualification case remain open.

## PR #1141 boundary

PR #1141 corrects outer-boundary reflection symmetry and tests that geometry
invariant. It does not, by itself, resolve M1–M4. Keep the symmetry fix
independently reviewable; implement the mesh generator, nonlinear solver and
face-traction changes as separately testable commits/PRs rather than mixing
unverified physical changes into the geometry regression.
