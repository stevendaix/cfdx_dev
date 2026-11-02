# CFDX lot 4 — S2S radiosity and Boussinesq body-force contracts

Date: 2026-11-01. Status: **locally implemented and unit/integration-tested**.
Base: local synthetic reconstruction of PR #1155 lots 1+2, plus the unpushed
lot 3 (CFL/transport + disconnected thermal components). This source tree is
**not** a verified checkout of the current GitHub PR head. This patch must be
checked against the real branch preimages before remote application.

## R01: true convergence and finite radiative power in S2S

The original `solve_s2s_radiosity` always set `S2SResult::converged=true` after
the direct Gaussian elimination. For finite positive input temperatures whose
fourth power cannot be represented, the solver could return `inf`/`nan`
radiosity and radiation power while claiming convergence. Finite values such
as `area=1e308 m2` and `T=1000 K` also overflow the integrated surface power.

The direct solver solves `M J=b` where

- `M_ij = delta_ij - (1-epsilon_i) F_ij`,
- `b_i = epsilon_i sigma T_i^4 + (1-epsilon_i) G_external,i`,
- the net wall flux is `q_i=J_i-G_i` [W/m2].

After solving, evaluate **the original matrix** (not the triangular factors):

`eta = max_i |(M J-b)_i| / max(sum_j |M_ij J_j| + |b_i|, DBL_MIN)`.

Expose `radiosity_backward_error` and `worst_radiosity_row` on `S2SResult`.
`converged` now requires `eta <= controls.tolerance`. Reject non-representable
blackbody emission, RHS, radiosity, irradiation, surface power, accumulated
power and energy balance rather than returning a purportedly converged result.
Do not use the net power of an **open** enclosure as a zero-sum gate: an open
system exchanging with the ambient can and should have a nonzero net power.

Before correction (four-case test): 2/4 passed; two overflow scenarios failed.
After correction (five-case test): 5/5 passed, including an independent nominal
identity check and a strict `1e-20` threshold that rejects the direct solve at
observed true backward error ~`7.24e-17`.

These are **numerical safety/contract tests**, not a Monte-Carlo view-factor
qualification, a physical 3D enclosure verification or a timing benchmark.

## B01: transactional and finite Boussinesq acceleration

Production uses `a_c = -beta (T_c - T_ref) g_vector` (m/s2). Previously,
`build_cell_buoyancy_acceleration` resized and began overwriting `acceleration`
before checking every cell. A later invalid temperature left a partly updated
force field. The function also accepted vector-valued temperatures and could
silently publish `inf`/`nan` from finite but overflowing products.

Validate the scalar field shape, coupling coefficients and entire temperature
field including each derived 3D acceleration **before** modifying the output.
Only after full validation resize and fill the caller's vector. Use two passes
to reuse its storage instead of allocating a second O(n_cells) vector for each
external buoyancy iteration.

Before correction (five-case test): 1/5 passed. After correction: 5/5 passed.
The original positive hot/cold/reference-temperature sign test is preserved.
No arbitrary Boussinesq validity limit was introduced; separate verification
of the small density-variation approximation remains necessary.

## Integration and test instructions

```bash
cmake -S . -B build-lot4 -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON
cmake --build build-lot4 --target \
    test_s2s_radiosity_verified test_buoyancy_acceleration_contract \
    test_s2s_radiation_vv test_thermal_radiation_regression \
    test_production_s2s_energy_stage_vv \
    test_production_s2s_energy_outer_loop_vv \
    test_production_buoyancy_outer_loop_vv -j 3
ctest --test-dir build-lot4 --output-on-failure -R \
  '^test_(s2s_radiosity_verified|buoyancy_acceleration_contract|s2s_radiation_vv|thermal_radiation_regression|production_s2s_energy_stage_vv|production_s2s_energy_outer_loop_vv|production_buoyancy_outer_loop_vv)$'
```

Run the new test targets with UBSan as a separate sanity check. Do not assert
MPI validation without MPI-enabled builds. Neither this patch nor its tests
promote NACA, COUPLED, SST or the CHT physical qualification to PASS. Existing
PR #1154 (thermal interfaces) and #1155 should be reviewed for changed file
preimages, CI and compatibility before this patch is applied remotely.
