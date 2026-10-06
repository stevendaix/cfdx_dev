# N9-S9 — Rhie–Chow pressure–velocity consistency

N9-S9 closes the pressure–velocity consistency boundary with executable operator-level evidence. The tests are deliberately independent of a particular converged CFD solution.

## Qualification gates

1. **Linear-pressure consistency** — when the cell pressure gradient exactly matches the owner-to-neighbour pressure difference, the Rhie–Chow pressure correction vanishes and the authoritative face flux equals the velocity-interpolated flux.
2. **Pressure-offset invariance** — adding a constant to every pressure value does not change any conservative face mass flux.
3. **Alternating-pressure coupling** — an alternating collocated pressure mode is visible to continuity through the pressure-corrected face flux; it cannot become invisible merely because cell-centred velocities interpolate to zero.

The checks use the production `make_rhie_chow_mass_flux()` path rather than a duplicated test implementation.

## Total Validation

`test_n9_s9_rhie_chow_consistency` is labelled `n9;qualification;validation;long;validation-total`.

The test therefore belongs to the expensive `validation-total` campaign and must pass on the exact PR HEAD before N9-S9 is considered qualified.

## Non-goals

This step does not claim N9 closure. It does not replace the physical qualification matrix (N9-S7), the genuine external-flow campaign (N9-S8), or the pressure-reference policy campaign (N9-S10).

No tolerance relaxation, disabled test, silent fallback, or alternate implementation is permitted.
