# N9-S9 — Rhie–Chow pressure–velocity consistency

N9-S9 closes the pressure–velocity consistency boundary with executable operator-level evidence. The tests are deliberately independent of a particular converged CFD solution.

## Qualification gates

1. **Linear-pressure consistency** — when the cell pressure gradient exactly matches the owner-to-neighbour pressure difference, the Rhie–Chow pressure correction vanishes and the authoritative face flux equals the velocity-interpolated flux.
2. **Pressure-offset invariance** — adding a constant to every pressure value does not change any conservative face mass flux.
3. **Alternating-pressure coupling** — an alternating collocated pressure mode is visible to continuity through the pressure-corrected face flux; it cannot become invisible merely because cell-centred velocities interpolate to zero.

The checks use the production `make_rhie_chow_mass_flux()` path rather than a duplicated test implementation.

## Majumdar relaxation treatment (derivation, 2026-10-16)

The production RC term of PR #895 weights the face pressure flux with the
RELAXED momentum diagonals: D^RC = V/a_P^r where a_P^r = a_P/alpha_u, while
the pressure-correction velocity update uses the unrelaxed operator
D^pc = V/a_P. Hence **D^RC = alpha_u * D^pc** (alpha_u = 0.65 in the
campaigns: the RC coefficients are ~35 % smaller than the p' divergence
coefficients).

Derivation: with implicit under-relaxation the momentum equation solved is
u_P = alpha_u * (H_P/a_P - V/a_P * grad p) + (1 - alpha_u) * u_P^old, so the
velocity correction inherits the relaxed matrix,
u' = -alpha_u * V/a_P * grad p'. The Majumdar-consistent collocated flux must
therefore interpolate the SAME relaxed inverse diagonal that produced the
cell velocities:
phi_f = rho [ (H/a)_f . S_f - D^RC_f * (dp/dn)_f ], with D^RC_f = V/a_P^r.
The pressure-correction divergence D^pc only weights the p'-gradient on the
ITERATION path: p' -> 0 at convergence, so it cannot move the fixed point -
its residual effect is O(h^2) on smooth modes. Consequences, both testable:

- operator level: scaling rAU (equivalently alpha_u) must scale the RC
  correction linearly (tests below);
- fixed point: the CONVERGED solution must be alpha_u-independent
  (empirical campaign v4 on the diag branch, alpha_u in {0.5, 0.65, 0.8,
  0.95} on the VMFL036 coarse and fine cases; diagnostic, not a gate).

The unrelaxed variant (D^RC = D^pc) was tested empirically on 2026-10-15 and
rejected: fine_err 10.06 % FAIL vs 9.93 % PASS for the relaxed form
(cf. VMFL036_FINE_GAP_STUDY.md).

Empirical campaign v4 (2026-10-16, twin `port4_rhie_chow.ts`, diag branch,
evidence `diag-evidence/logs/alpha_sensitivity_2026-10-16.json`): the
alpha_u-spread of the converged Cd is 1.44 % on the coarse 24x36 mesh and
0.19 % on the fine 48x72 mesh, with the alpha_u = 0.65 coarse run
reproducing the reference Cd bit-exact (selfCheckDelta = 0) and distinct
iteration counts per alpha_u (plumbing verified: 4/4). The residual
alpha_u-sensitivity decays with refinement, consistent with the O(h^2)
iteration-path effect predicted by the derivation: the fixed point is
alpha_u-independent up to discretisation, and the fine-level spread
(0.19 %) is two orders of magnitude below the drag band (10 %).

## Analytical tests (operator level, production path)

Five additional operator-level tests, all through the production
make_rhie_chow_mass_flux() (no duplicated implementation):

1. the RC correction scales linearly with rAU - the operator-level statement
   of D^RC = alpha_u * D^pc (Majumdar);
2. the correction is linear in p (superposition of an alternating field and a
   linear field);
3. directional rAU anisotropy: only the face-normal projection of the
   directional inverse diagonals enters rfn;
4. stretched mesh with unequal cell volumes: the face coefficient is the
   volume-weighted dAU = mean(V*rAU) interpolation;
5. uniform-mesh limit: the production correction collapses exactly to
   rho * dAU_f * (p_N - p_P)/d * (S_f . e).

## Total Validation

`test_n9_s9_rhie_chow_consistency` is labelled `n9;qualification;validation;long;validation-total`.

The test therefore belongs to the expensive `validation-total` campaign and must pass on the exact PR HEAD before N9-S9 is considered qualified.

## Non-goals

This step does not claim N9 closure. It does not replace the physical qualification matrix (N9-S7), the genuine external-flow campaign (N9-S8), or the pressure-reference policy campaign (N9-S10).

No tolerance relaxation, disabled test, silent fallback, or alternate implementation is permitted.
