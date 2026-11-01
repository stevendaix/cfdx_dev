# CFDX lot 5 — scale-invariant S2S, finite forces and wall-turbulence guards

Date: 2026-11-01. **Implementation tested locally; no CFD model promoted to qualified.**

## Provenance and integration boundary

The changes are based on the **local synthetic lot-4 snapshot** (`f5c2617`),
constructed from the 2026-10-31 archive and the prior local patch sequence.
It is **not** a verified checkout of the GitHub PR #1155 HEAD. The owner
reported PR #1155 at `354a6b95` after lot 3. Apply lot 4 first, then this lot,
only after comparing preimages with the actual remote HEAD. Neither a remote
push nor a remote CI result is claimed here.

## S2S-02 — reciprocity independent of surface-area units/scale

The existing view-factor guard compared
`abs(A_i*F_ij - A_j*F_ji)` to `tolerance * max(1, abs(A_i*F_ij), abs(A_j*F_ji))`.
The implicit **1 m²** scale made an arbitrarily large relative reciprocity
error acceptable for micro-scale patches. Example:
`A = [2e-10, 4e-10] m²`, `F12 = F21 = 0.5` gives unequal
products `1e-10` vs `2e-10` but passed with tolerance `1e-10`.

The new check evaluates both products after division by
`max(A_i,A_j)`, then compares with relative tolerance against the larger
scaled exchange. This is invariant under scaling all areas and does not form
possibly overflowing `A*F` products. Zero exchange in both directions passes;
one-sided nonzero exchange fails. **No numerical tolerance was increased.**

An additional physical oracle solves the same open two-surface enclosure after
uniformly scaling both areas by `1e-10`. It requires unchanged radiosities and
heat flux densities and power proportional to the area scaling.

## FORCE-01 — finite aerodynamic forces and moments

`forces::integrate` accepted a NaN moment center, area/pressure products that
became infinite, and invalid `q*A` / `q*A*L` coefficient normalizations.
`integrate_axisymmetric` similarly accepted force overflow and nonrepresentable
`q*A`. The contract now validates the moment reference, both coefficient
normalizations, per-face force/moment terms, running sums, and all returned
projections and coefficients. No clipping, silent substitution or change of
physical sign convention was introduced. Checks are O(n_faces), no additional
field-sized allocations.

The general surface-force integration helper is not proof of a corrected
production NACA pressure/viscous-force reconstruction: that validation still
requires the production face state and converged three-mesh campaign. The
axisymmetric force routine is used by `axisymmetric_vmfl036_solver.h`.

## WALL-02 — representative y+ and DES diagnostics

`wall_distance_y_plus` previously calculated `d*u_tau/nu` in that order.
It could return `inf` for a mathematically representable y+ due to an
intermediate overflow, or zero due to an intermediate underflow. The common
fast double path remains; only problematic cases try `(d/nu)*u_tau` and an
extended-precision fallback when available. Truly non-representable y+ throws.
`des_wall_transition_ratio` now rejects NaN/Inf `Cdes` and non-finite results.
DDES rejects a non-finite `kappa` or viscosity ratio. The wall-distance grid
contract also checks the `solid` mask size before indexing.

**This hardens diagnostic inputs, not the actual SST transport solution**. A
valid physical SST solution still needs canal/flat-plate/BFS benchmark evidence,
independent wall-distance accuracy, and local y+ resolution criteria.

## Executable negative and positive evidence

- `test_s2s_area_scale_reciprocity`: invalid micro/macro reciprocity, scaled
  valid exchange, physical power/flux scaling.
- `test_force_finite_contract`: invalid origin, nonrepresentable force and
  moment normalizers, planar/axisymmetric face force overflow, nominal drag.
- `test_wall_turbulence_contract_guard`: y+ intermediate overflow/underflow,
  nonfinite DES constant, invalid wall-grid mask, valid diagnostic samples.

All three suites failed on the preceding code before their targeted fixes.
Selected downstream Release, sanitizer and production regression evidence is
recorded with the integration artifact; do not infer other suites passed.

## Remaining blockers (unchanged)

1. COUPLED/Rhie–Chow and dual-time first-step contract: inspect real residuals,
   mass conservation and NACA nonlinear failure, after reviewing relevant PRs.
2. NACA0012: converge **three** actual meshes, qualify Cp, Cf, Cd_p, Cd_v,
   Cd_total and Cl independently, apply GCI/observed-order tests and references.
3. Turbulent RANS SST: end-to-end channel, flat plate and backward-facing
   step; y+, wall stress, turbulence quantities and separation metrics.
4. Thermal/CHT: conjugate interface flux identity, energy balance by component,
   transient accuracy and coupled radiation/energy convergence.
5. Reliability/performance: actual MPI/GPU/CPU equivalence, peak memory and
   wall-clock time to the **same physical convergence**, not a per-iteration
   speed claim.

This is an implementation/verification lot, not a release qualification claim.
