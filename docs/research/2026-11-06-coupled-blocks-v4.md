# CFDX V4 — Continued pressure–velocity implementation and qualification

## Baseline and scope

Base: verified V3 full-source archive `CFDX_couplage_PV_briques_v3_sources_completes_2026-11-06.zip`. This is the V4 continuation delivered on PR #1183, based on the verified V3 source archive; it does not modify PR #1182. Existing algorithms and physical tolerances remain unchanged by default.

## Changes to executable numerical code

1. **Consistent pressure-mass PCD requires an SPD certificate** (`src/cfdx/core/linalg/pcd_schur.h`). Symmetry and positive diagonal alone do *not* imply positive definiteness: the symmetric matrix `[2,3;3,2]` has eigenvalues `5` and `-1`. For up to 96 pressure unknowns, validate via Cholesky (with a relative pivot safety floor). For larger masses, only accept the sufficient condition of symmetric *strict* diagonal dominance with positive diagonal; valid but non-dominant SPD masses can be conservatively rejected. The legacy inverse action remains unchanged.

2. **Schwarz is now selectable in production COUPLED via an opt-in** (`pressure_velocity_algorithms.h`, `steady_incompressible_solver.h`). Use `CouplingControls::coupled_schwarz_experimental` with a BlockLocal Schur model, the resolved CoupledBlockSchur linear plan and a serial, bounded local LU. Controls specify the core cells per subdomain, graph-overlap hops, max local unknowns (4–256, multiple of 4), and restricted or multiplicity-weighted assembly. Every unsupported configuration or failed local factorization fails closed with no silent AMG fallback. The exact existing default is retained; no MPI or mesh-independent performance claim.

3. **Small non-SPD Schur reference oracle** (`src/cfdx/core/linalg/small_dense_schur_oracle.h`). A row-scaled, pivoted dense LU with recomputed true residual builds and inverts `S=C-D Auu^-1 G` for small systems (maximum 64 by default). Both nonsymmetric and negative-definite Schur systems are supported, avoiding an unconditional CG/SPD assumption. Setup/update are transactional for normal invalid inputs, with graph-preserving numeric refresh. This is a diagnostics-only reference and NOT a scalable solver.

## Regressions and measured evidence

- The new consistent-mass PCD indefinite matrix counterexample is rejected; a previous valid SPD mass still operates after a rejected numeric refresh.
- `test_pcd_schur`: **22/22 cases**.
- `test_small_dense_schur_oracle`: **5/5 cases**, including negative-definite pressure Schur and an independently rechecked true residual.
- `test_steady_incompressible_solver`: **33/33 cases**, including invalid COUPLED Schwarz configuration rejection.
- `test_n9_s9_rhie_chow_consistency`: **19/19 cases**, including a four-cell manufactured nonzero-flow example solved through the opt-in **two overlapping Schwarz subdomains** and the real post-gauge Schur reference. Measured conservation `Linf=1.11022e-15`, Schur defect `3.68219e-16`, independent velocity solve defect `4.30211e-16`.
- The seven CMake-built unit executables and the two separately compiled `-O0` validation executables pass under **CTest 9/9**. The full `-O2` build of the two large validation translation units hit local time limits and was **not** counted as a successful optimized CMake build. The `cfdx_core` V4 library and seven unit binaries were compiled successfully by CMake.
- The PCD (22), dense oracle (5), and Schwarz (9) suites also pass under ASan+UBSan, **36/36 cases**. Those builds sanitize the unit executables, not all libraries or MPI code.

## Release blockers, still open

- Physical convergence on NACA N64/N128/N256 and three-grid force/skin-friction/GCI qualification.
- N9 extended/Ghia Re=1000, full Release & DebugSanitizers, PyFluent and MPI validation.
- Production nonlinearity: a passing linear Schur inverse probe is neither proof of a correct nonlinear Jacobian nor proof of pressure–velocity equivalence across algorithms.
- Compare Schwarz, PCD, LSC, BFBt and existing AMG on *the same* stored nonzero-flux production matrices and with true residual and wall-clock acceptance gates. RAS without a coarse space is not expected to scale to large or highly anisotropic meshes.
- Future work: link the directional Jacobian oracle to an independently assembled production FV residual with correctly treated deferred fluxes and boundary conditions; extend small Schur comparison to actual PCD/LSC inverse actions; prohibit any automatic selection until validated.

## Usage (experimental; not a recommended production default)

```cpp
IncompressibleSolverControls controls;
controls.algorithm = PressureVelocityAlgorithm::COUPLED;
controls.coupling.schur_model = CoupledSchurModel::BlockLocal;
controls.coupling.coupled_schwarz_experimental = true;
controls.coupling.coupled_schwarz_cells_per_subdomain = 2;
controls.coupling.coupled_schwarz_overlap_hops = 1;
controls.coupling.coupled_schwarz_max_local_unknowns = 16;
controls.coupling.coupled_schwarz_weighted = false; // restricted RAS
```

For a small, captured gauge-fixed block matrix, `SmallDenseExactSchurInverse` provides a reference action to `certify_schur_inverse_action()`. Do not use it for NACA-scale solves.
