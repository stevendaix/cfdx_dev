# 3D incompressible solver numerical audit — 2026-10-16

## Scope and evidence boundary

Audit target: the general 3D finite-volume numerical path in `master` at base commit
`bfaf71071c72f6a9f37391039826a5964787a2a5`. The axisymmetric VMFL036 solver is
out of scope here except as motivation.

This is a source-level audit and evidence-gap report. No production numerical
operator, tolerance, convergence threshold, mesh, or scheme default is changed
by this document. The findings below separate source facts and already
documented measurements from hypotheses that still require a solver-level run.

## Findings

### F1 — High risk: default cell gradient is not linearly exact on general tetrahedra

**Source facts**

- `src/cfdx/core/numerics/gradient.h`: `compute_gradient_gauss` interpolates
  an internal-face value using a scalar weight based on distances from the
  face centroid and cell centres, then applies Green–Gauss.
- This interpolation reconstructs a linear field at the face centroid only
  when the face centroid lies on the owner–neighbour centre line (or for
  special affine geometries). A general tetrahedral face does not satisfy that
  condition.
- `compute_gradient_gauss` is used directly by several paths, including
  `src/cfdx/physics/incompressible.h` for pressure gradients and
  `src/cfdx/physics/finite_volume_transport.h` for TVD reconstruction.
- `compute_laplacian` and `compute_convection` default their gradient scheme
  to `GradientScheme::GAUSS_TWO_POINT`; an explicit configuration may select
  another registered scheme, so the actual case selection must be traced before
  attributing this default to a particular run.

**Existing measured evidence**

The numerical-method contract document already records a tetrahedral gradient
campaign in which the two-point Gauss gradient has observed order near 0.1 on
a smooth field, while the quadratic least-squares option reaches approximately
second order on the qualified tetrahedral family. It also explicitly states
that the two-point face interpolation is not linearly consistent on general
tetrahedra.

**Impact**

This is a mesh-dependent accuracy risk for gradients entering pressure forces,
convection reconstruction, and explicit non-orthogonal corrections. It is not
evidence that all 3D cases are wrong: affine/hexahedral meshes can behave
substantially better, and explicit scheme selection changes the path.

**Required closure**

1. Make the production selection report available in a solver-level test for
   each audited case; do not infer it from API defaults.
2. Apply linear and smooth manufactured fields on the actual mesh family and
   report L1/L2/Linf errors and observed order.
3. Gate the selected production path, not just the standalone gradient kernel.
4. If the default is retained for compatibility, reject or prominently qualify
   its use on mesh families where the required consistency is not met; do not
   silently switch schemes based on mesh shape.

### F2 — High risk: corrected polyhedral Laplacian has a documented accuracy plateau

**Source facts**

- `src/cfdx/core/numerics/laplacian.h`: the orthogonal contribution uses
  `(Sf·d)/|d|² * (phi_N - phi_P)`. Corrected/limited/over-relaxed variants add
  `Sf_corr · grad_f`, where `grad_f` is currently the arithmetic average of
  cell-centred gradients, i.e. evaluated at the cell-centre midpoint rather
  than reconstructed at the geometric face centroid.
- The default gradient argument is `GAUSS_TWO_POINT`. The corrected path
  therefore inherits the selected gradient's consistency limits and has a
  separate face-location error on non-affine meshes.

**Existing measured evidence**

`docs/development/NUMERICAL_METHOD_CONTRACTS.md` records that on tetrahedral
meshes the smooth-field Laplacian remains O(1) under refinement even when a
quadratic least-squares gradient is used. It also records a later test of
quadratic-Hessian extrapolation to the face centroid that still plateaued:
the two-point implicit term itself carries a curvature mismatch. The document
therefore correctly classifies this as a deeper discretisation issue, not a
gradient-scheme-only fix.

**Impact**

Conservation across a shared internal face is necessary but does not establish
consistency or convergence. A conservative corrected Laplacian can still fail
the expected smooth-field order on non-affine polyhedra.

**Required closure**

1. Preserve the existing skew-ladder conservation tests.
2. Add/retain smooth MMS refinement tests on tetrahedra and non-affine
   polyhedra, with observed order and a strict error-decrease gate.
3. Do not claim the polyhedral diffusion family qualified while the O(1)
   plateau remains.
4. Treat a new face-normal reconstruction including curvature terms as a
   separate discretisation change with an independent derivation and oracle.

### F3 — Medium/high risk: Rhie–Chow unit tests do not by themselves qualify the 3D solver loop

The production helper `make_rhie_chow_mass_flux()` has analytic operator tests
in `tests/validation/test_n9_s9_rhie_chow_consistency.cpp`. These are valuable
kernel checks, but they do not alone establish that every pressure-velocity
algorithm assembles its pressure-correction operator with coefficients and
boundary terms consistent with the authoritative face flux.

A concrete source-level mismatch was confirmed in the current pressure
assembly. The Rhie–Chow implicit flux uses the centre-line unit vector
`e = d/|d|`, the projected area `Sf·e`, and the momentum response projected
onto `e`. Before PR #923, the pressure matrices instead used `|Sf|/|d|` and
projected the directional response onto the face normal. Those expressions
coincide on orthogonal faces but differ on skew faces. The mismatch was present
in segregated SIMPLE/SIMPLEC, the COUPLED block, fixed-pressure boundaries, and
the extra SIMPLEC flux term. PR #923 now routes these coefficients through a
shared operator and adds an analytic skew-face regression; CI and physical
validation are still pending.

The separate method sheet `docs/audits/03-interpolation-rhie-chow.md` now records
this contract. Previous CFDX diagnostics also found flux/velocity mismatch and
non-zero continuity on stationary states in particular paths; those historical
measurements still need reproduction on the current HEAD before being treated
as current regressions.

**Required closure**

For each supported 3D pressure-velocity algorithm (start with SIMPLE/SIMPLEC;
then audit COUPLED and PISO/PIMPLE separately):

- capture the selected scheme and algorithm in the run artifact;
- compute the divergence of the exact face flux used by the pressure equation;
- compare that flux with an independently reconstructed velocity flux, by
  internal/boundary face and by mesh family;
- recompute momentum and continuity residuals on the final corrected state;
- test pressure-offset invariance, a linear pressure field, skewed faces,
  fixed-pressure boundaries, and the relaxation sweep;
- retain the existing physical and `validation-total` gates unchanged.

No conclusion is made here that every algorithm currently fails these tests.

### F4 — Medium risk: convergence metrics need a final-state contract

The source already reassembles the three momentum equations after the final
pressure correction and computes the final momentum residual; it also computes
continuity from the authoritative `mass_flux` and separately compares that
flux with a reconstructed cell-velocity flux. Therefore, the stale-momentum-
residual concern identified in the axisymmetric solver should **not** be
assumed to apply to the general 3D solver.

The remaining audit is to verify the meaning and scaling of these existing
metrics, ensure all algorithms (especially COUPLED) report the same final-state
contract, and assert the flux mismatch on skew and polyhedral meshes. A small
linear residual or pressure-correction increment remains distinct from a small
nonlinear residual.

## Priority and proposed implementation sequence

| Priority | Work | Acceptance evidence |
|---|---|---|
| P0 | Re-run existing tetrahedral gradient and polyhedral Laplacian campaigns on the exact audited HEAD | Full logs, mesh family, scheme IDs, errors and observed orders |
| P0 | Trace the actual production numerical selections for one hexahedral, one tetrahedral and one non-affine polyhedral 3D case | Case selection report plus source call path |
| P1 | Add a solver-level 3D manufactured test for pressure/flux/velocity consistency | Fail-closed assertions on the authoritative flux, boundary flux and final-state residuals |
| P1 | Audit SIMPLE/SIMPLEC pressure-correction coefficients against the flux operator | Face-by-face coefficient/flux identity, including fixed-pressure boundaries |
| P1 | Verify the existing final-state residual and flux-mismatch diagnostics on skew/polyhedral meshes and each algorithm | Residuals refer to the final corrected state; authoritative and reconstructed fluxes are reported separately |
| P2 | Design a consistent polyhedral diffusion discretisation | Independent derivation, MMS convergence, conservation, and skew/refinement matrix |

## Change policy

No global default is changed in this audit commit. Switching all 3D meshes to
least-squares or quadratic least-squares without tracing the production path
could change validated hexahedral results, boundary behavior, cost, and
convergence. The polyhedral Laplacian plateau also cannot be repaired merely by
changing the gradient argument. Each correction should be a small follow-up
change with a failing regression reproduced first, no tolerance inflation, and
the existing validation gates preserved.

## References within the repository

- `src/cfdx/core/numerics/gradient.h`
- `src/cfdx/core/numerics/laplacian.h`
- `src/cfdx/physics/incompressible.h`
- `src/cfdx/physics/finite_volume_transport.h`
- `src/cfdx/physics/steady_incompressible_solver.h`
- `tests/validation/test_polyhedral_gradient_campaign.cpp`
- `tests/validation/test_polyhedral_laplacian_campaign.cpp`
- `tests/validation/test_n9_s9_rhie_chow_consistency.cpp`
- `docs/development/NUMERICAL_METHOD_CONTRACTS.md`
- `docs/audits/03-interpolation-rhie-chow.md`
