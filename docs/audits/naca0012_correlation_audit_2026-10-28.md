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
- Add reproducible mesh-quality reporting: first-cell wall-normal distance,
  growth ratio, cell Jacobian/area positivity, aspect ratio, skewness and
  non-orthogonality.
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

## PR #1141 boundary

PR #1141 corrects outer-boundary reflection symmetry and tests that geometry
invariant. It does not, by itself, resolve M1–M4. Keep the symmetry fix
independently reviewable; implement the mesh generator, nonlinear solver and
face-traction changes as separately testable commits/PRs rather than mixing
unverified physical changes into the geometry regression.
