# Lot 3B - anchored connected components in steady conduction (2026-11-01)

Status: reproduced on a topologically valid, disconnected Cartesian mesh,
fixed in the local snapshot, and covered by an independent CTest target.
No PR was modified and no physical model was promoted.

## Counterexample

Three cells in two closed Cartesian regions: one isolated cube with a
fixed-temperature boundary and a second region containing two connected
cubes with adiabatic boundaries only. The previous `has_fixed_temperature`
boolean was global. It let the second region through even though its
independent constant-temperature gauge remained free.

The assembled 3x3 matrix admitted the null vector `[0,1,1]`:
`A * [0,1,1] = [0,0,0]` (actual baseline result, not theoretical only).

## Correction

After validating face conductances and checking at least one global Dirichlet
condition (existing contract), perform an O(n_cells + n_faces * alpha(n_cells))
union-find traversal of internal faces. Mark each connected component with at
least one non-empty fixed-temperature face. Reject any component that is
unanchored with a diagnostic identifying a representative cell.

This does **not** introduce a pressure-like arbitrary gauge; it preserves the
repository's requirement of an explicitly prescribed steady temperature.
No conductivity, source, interface-resistance or boundary-flux equation was
changed. The thermal interfaces PR (#1154) may require a rebase against this
file; its remote file list could not be inspected from this environment.

## Verification

- `test_thermal_component_anchoring`: disconnected unanchored component
  rejected; separate components with Dirichlet assembled; connected two-cell
  domain with one Dirichlet assembled; pure Neumann still rejected.
- Existing `test_thermal_zones_vv` to be rerun after compilation.
- The pre-fix reproducer `CFDX_REPRO_THERMAL_DISCONNECTED_2026-11-01.cpp`
  exited 1 with a zero matrix-vector product for a nonzero vector.

This is a topological/discrete well-posedness test, not a benchmark asserting
full thermophysical or 3D CHT qualification.
