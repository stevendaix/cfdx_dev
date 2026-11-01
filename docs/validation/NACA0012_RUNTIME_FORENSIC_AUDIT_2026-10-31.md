# NACA0012 runtime forensic audit — 2026-10-31

## Scope and provenance

- Pull request: #1149
- Observed head: `a1b1d490cdf390c805c87bdd557cc24102d54044`
- Case: incompressible laminar NACA0012, Re=1000, zero incidence, SIMPLE,
  UPWIND convection, generated O-grids N64/N128/N256.
- This is a failure-attribution record, not a physical validation claim.
- No solver tolerance, iteration budget, relaxation, fallback, or physical gate
  is changed by this audit.

## Exact-head long-campaign evidence

The long validation job generated all three meshes successfully:

| Level | Cells | Maximum non-orthogonality | P95 non-orthogonality |
|---|---:|---:|---:|
| N64 | 3,096 | 53.57 deg | 20.09 deg |
| N128 | 12,336 | 65.96 deg | 21.43 deg |
| N256 | 49,248 | 73.24 deg | 21.93 deg |

All three generated meshes reported zero inverted or degenerate cells.
The increasing maximum non-orthogonality and aspect ratio remain relevant
geometry covariates; they are not by themselves proof of the failure cause.

N64 completed in 456 nonlinear iterations:

- Cd = 0.0901496
- Cl = -9.15e-11
- Cd_pressure = 0.0344606
- Cd_viscous = 0.055689
- continuity = 8.84e-11

N128 failed during the momentum predictor. The failure was not a mesh-file
generation error and did not reach N256.

## Linear-system failure snapshot

At the N128 Ux predictor failure, the emitted scalar linear diagnostics were:

| Quantity | Value |
|---|---:|
| Matrix rows | 12,336 |
| Matrix nonzeros | 61,166 |
| Matrix infinity norm | 2.52932e13 |
| Matrix symmetry-defect infinity norm | 2.26008e12 |
| RHS infinity norm | 4.15478e25 |
| Diagonal minimum | 4.78402e-2 |
| Diagonal maximum | 2.29938e13 |

The fallback attempts reported:

- GMRES: status 1, 4,000 iterations, residual 6.53586e20.
- Gauss-Seidel: status 1, 4,000 iterations, residual 2.63031e17.
- CG: status 3 on its first iteration; pAp = -9.40115e45 and alpha = NaN.
- Final Ux solve: status 3, one iteration, residual 1.329e25.

The symmetry defect divided by the matrix infinity norm is approximately 8.94%.
This is a convection-diffusion momentum matrix, so symmetry is not generally
expected and CG is not an appropriate method for this operator. The negative
pAp is therefore a useful incompatibility signal, not evidence that CG caused
the physical divergence. GMRES and Gauss-Seidel also fail after the RHS and
matrix scales have become extreme.

The large RHS (4.15e25) and coefficient range are consistent with a nonlinear
state that has already become numerically damaged before the final linear
breakdown. The current log does not identify the first iteration or the first
assembled face/cell contribution responsible. Do not infer a discretization
bug from the final matrix statistics alone.

## Wall-shear audit classification

N64 wall-shear values were:

| Reconstruction | Cd_visc |
|---|---:|
| Gauss-Green | 0.055689 |
| Weighted least squares | 0.0580625 |
| 1-D wall-normal estimate | 0.106308 |

The relative spread is 0.690069. The spread is an observation, not a gate.
The wall-shear campaign then failed while solving the next NACA level with the
same Ux momentum breakdown. Therefore this run does not establish an
independent wall-shear implementation failure: the full three-level
post-processing comparison is blocked by the flow-solver failure.

The Ghia Re=100/400/1000 cavity campaign passed in 2,742.19 seconds. This is
independent positive evidence for that campaign only; it does not validate
NACA or the momentum operator on the O-grid.

## Diagnostic instrumentation added for the next exact-head run

The opt-in `NacaMomentumTrace` now records, at every actual momentum predictor
attempt (including nonlinear retries):

- maxima and first non-finite indices for pressure, pressure gradient,
  velocity, face mass flux, cell flux divergence, and effective viscosity;
- the component and cell associated with maximum pressure-gradient and
  velocity magnitudes;
- the maximum-flux face together with owner and neighbour cell IDs;
- the cell with maximum integrated flux divergence;
- boundary-flux sum and absolute boundary-flux sum.

The long-validation workflow sets `CFDX_NACA_FORENSICS=1`, so the next
run retains this trace without changing the numerical path. The callback is
read-only and is null unless the forensic environment flag is enabled.

## Required next analysis

1. Compare the first N64 and N128 iterations where pressure, velocity,
   gradient, or flux magnitudes depart from their expected O(1) scales.
2. Correlate the first abnormal cell/face with O-grid geometry: face
   non-orthogonality, owner-neighbour centre vector, face area, and wall
   distance.
3. If the input fields remain bounded but the assembled matrix/RHS jumps,
   add an equation-level snapshot after `assemble_momentum_component` to
   separate convection, orthogonal diffusion, non-orthogonal correction,
   pressure-gradient source, and boundary contributions.
4. Only after that attribution, prepare a minimal A/B test of the implicated
   term. Keep the current discretization and all qualification gates unchanged
   until the evidence isolates a defect.

## Explicit non-claims

- N64 passing is not mesh convergence or physical reference validation.
- The wall-shear spread is not a reason to select the largest estimate.
- Ghia passing does not clear the NACA failure.
- The negative CG curvature is not the primary root cause of the NACA run:
  the nonsymmetric momentum operator and the failed nonsymmetric-capable
  fallbacks must be considered first.


## Follow-up pressure-equation forensics (2026-11-01)

The next diagnostic pass adds a per-cell decomposition of the pressure RHS.
For each pressure corrector, `NACA_PRESSURE_RHS` reports the maximum RHS
cell and the signed contributions at that cell from:

- predictor mass flux `phiHbyA`;
- implicit orthogonal pressure flux evaluated from the current pressure;
- deferred non-orthogonal pressure flux;
- fixed-pressure boundary source.

The reported `decomposition_error_linf` checks that these independently
reconstructed contributions sum to the assembled RHS. This is an internal
diagnostic consistency check, not a physical acceptance gate.
`NACA_PRESSURE_FACE` records the incident faces of the maximum-RHS cell
during the first system and every tenth nonlinear iteration, including owner,
neighbour, patch, predictor flux, implicit/deferred pressure flux and boundary
source. The summary is emitted on every pressure corrector.

### Exact pressure matrix snapshot

The pressure solver uses a gauge-reduced principal submatrix for the default
pure-Neumann/reference-cell path. Earlier diagnostics inspected the full-size
local matrix object, which is not the matrix passed to Krylov in that path.
The forensic hook now retains the actual solve matrix and RHS, including the
reduced objects, and all matrix summaries refer to those exact arguments.

To save the first pressure system for each mesh, enable the forensic hook with
`CFDX_NACA_FORENSICS=1` and set `CFDX_NACA_MATRIX_DUMP_PREFIX` to a writable
file prefix when running the NACA validation executable. The hook appends the
cell count to the prefix so sequential N64/N128/N256 runs cannot overwrite one
another. It writes:

- `<prefix>.cells<cell-count>.mtx`: the exact CSR pressure operator passed to
  Krylov, in Matrix-Market coordinate format;
- `<prefix>.cells<cell-count>.rhs.csv`: the matching RHS and the reduced-row
  to physical-cell mapping.

For example, use a distinct prefix for each mesh:
`CFDX_NACA_FORENSICS=1` and
`CFDX_NACA_MATRIX_DUMP_PREFIX=artifacts/naca-pressure-first`.
The output directory must already exist. No matrix file is written when either
the forensic hook is disabled or the prefix variable is unset.

Run the Jacobi-scaled spectrum diagnostic with:

```bash
python3 scripts/validation/analyze_naca_pressure_spectrum.py \
  artifacts/naca-N256-pressure-first.mtx \
  --rhs artifacts/naca-N256-pressure-first.rhs.csv
```

The analyzer reports matrix symmetry, graph connectivity, diagonal range,
and (only for a connected symmetric matrix with positive diagonal) the
extreme eigenvalues of `D^{-1/2} A D^{-1/2}`, their residuals, and the
resulting condition estimate. It does not modify the solver or promote any
qualification status. NumPy and SciPy are required for the spectrum step.

### Interpretation protocol

Compare N64, N128 and N256 snapshots from the same solver revision and
initialization. First inspect `decomposition_error_linf`; it should be near
round-off relative to the RHS scale. Then compare the four RHS contributions
at the maximum-RHS cell and its incident faces. A large predictor term points
upstream toward HbyA/face interpolation; a large deferred term points toward
the explicit non-orthogonal pressure response; a large boundary source points
toward fixed-pressure boundary reconstruction. Only if these terms are
well-scaled while the correction remains anomalous should the matrix spectrum
and local diagonal/RHS ratios become the primary focus.

These diagnostics remain evidence-only. No relaxation, discretization,
iteration limit, tolerance, or qualification gate is changed.


The two long NACA CTest entries use test-local environment overrides
`CFDX_NACA_AB_ALPHA_U=0.60` and `CFDX_NACA_AB_ALPHA_P=0.10` to reproduce
the controlled N256 experiment reported on 2026-11-01. These are diagnostic
run settings only: the solver defaults are unchanged, the physical gates are
unchanged, and this run cannot promote NACA qualification. The overrides are
kept local to the two NACA long tests; other CFD validation cases do not inherit
them.
