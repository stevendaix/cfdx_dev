# CFDX numerical-method contracts

This is the first architecture slice of Issue #461. The goal is to make a numerical method an explicit, traceable object rather than an implicit collection of defaults spread across solver code.

## Contract layers

A numerical method is described by five independent layers:

1. Mathematical formulation — what operator or algorithm is actually applied.
2. Configuration identity — the case/numerics key that selects it.
3. Numerical invariants — conservation, boundedness, monotonicity, consistency and null-space requirements where applicable.
4. Verification evidence — exact/operator tests, MMS, refinement and regression requirements.
5. Validation evidence — physical/reference cases. Validation never replaces verification.

The C++ contract header contains metadata and metadata validation only. It intentionally does not make an implementation PASS because it has been registered.

## Required fields

| Field | Meaning |
|---|---|
| id | Stable machine-readable method identifier |
| name | Human-readable name |
| family | Gradient, interpolation, convection, diffusion, temporal, solver, etc. |
| status | Planned / implemented / verified / validated |
| conservation | Declared conservation property |
| bounded | Whether a boundedness contract is claimed |
| monotone | Whether monotonicity is claimed |
| linear_exact | Whether linear-field exactness is claimed |
| formal order | Declared spatial/temporal order |
| formulation | Mathematical definition |
| configuration key | Explicit selection path |
| verification requirements | Minimum evidence required for promotion |

## Promotion rule

Implemented means executable code exists.

Verified requires executable quantitative evidence appropriate to the method, including an independent oracle or MMS/refinement evidence where applicable.

Validated requires an additional application or physical/reference comparison. A unit test of a limiter coefficient, for example, is not sufficient to promote a convection scheme to verified.

## Selection resolution

`core/numerics/numerical_method_selection.h` turns the registry `configuration key` from
metadata into an enforced selection path. A case resolves a scheme with
`select_numerical_method(family, configuration_key)`:

- the key must be non-empty — an empty key throws instead of applying a hidden default;
- the key must name exactly one registered method — unknown or ambiguous keys throw;
- the key must belong to the requested family — cross-family keys throw;
- the resolved `NumericalSchemeSelection` carries the method id and verification status, and
  `format_scheme_selection` renders a deterministic report line.

The resolver is total over the registry: every registered method round-trips from its own key.
This is the machine-checkable form of the rule "configuration must expose the selected method".

For a whole case configuration, `resolve_case_numerics(case_name, entries)` resolves every
`(family, key)` entry and produces a `CaseNumericsReport`: accepted entries give auditable
`NumericalSchemeSelection` records, rejected entries (unknown / empty / wrong-family key) are
recorded as errors — **no silent fallback, no throw**. `format_numerics_report` renders a
deterministic plain-text report (one `scheme[...]` line per resolved scheme plus `error:` lines),
which is what a case-loading or validation pipeline records. `tests/unit/test_case_numerics_report.cpp`
covers determinism, error trapping, cross-family rejection and a full-registry round trip at the
case level. The native CFDX HDF5 case path now supports a canonical `numerics.selection` block. It contains explicit registry configuration keys plus a case-local list of required families. `read_case_cfdx_h5` resolves that block before returning a loaded case and rejects empty, unknown, cross-family, missing-required or duplicate-family selections. The same contract is enforced by `write_case_cfdx_h5`, and the deterministic resolved report is persisted as `numerical_selection_report` in the case artifact. Legacy HDF5 files without the explicit selection block remain readable during migration; they are not considered N1-qualified until migrated to the explicit path.

The Python adapter normalization boundary now populates `numerics.selection` for the numerical families currently consumed by the production runtime. The production solver loads `.cfdx.h5` through the case loader and consumes the resolved pressure-velocity, convection and supported gradient selections; explicit linear-solver and preconditioner selections are consumed when present. Unsupported runtime selections are rejected rather than silently falling back. Remaining qualification work is end-to-end CI evidence and source-specific extraction for linear-solver choices.


## Common spatial-scheme contract

`core/numerics/spatial_scheme_contract.h` provides the common N1 vocabulary for gradient, reconstruction, interpolation, convection, diffusion and flux roles. It deliberately separates scheme identity (the registry configuration key) from the concrete operator implementation. A spatial contract declares formulation, formal order, boundedness/monotonicity and conservation claims; invalid bounded-convection claims are rejected.

The unified Python conversion pipeline is the single adapter boundary for Fluent, SU2, OpenFOAM, Star-CCM+ and Code_Saturne. It canonicalises the normalised `CaseSetup.numerics` into `numerics.selection` for every adapter. A source-specific numerical value that has no canonical mapping is now a blocking conversion gap rather than an implicit default. The HDF5 writer persists the deterministic selection report, and the native case loader resolves the same keys against the C++ registry before solver startup.

The spatial contract is architectural evidence, not a numerical qualification: quantitative verification remains owned by N2/N3/N4 and the corresponding V&V campaigns.

## Initial migration targets

The next implementation waves should register the existing numerical families without changing their algorithms first:

- Gauss gradient and future least-squares gradients;
- linear/upwind/limited interpolation;
- TVD limiter family;
- orthogonal/corrected/limited/uncorrected diffusion;
- Euler/CN/BDF2 temporal schemes;
- CFL and pseudo-time controllers;
- Krylov solvers and native AMG/FieldSplit/Schur preconditioners;
- SIMPLE/SIMPLEC/PISO/PIMPLE/fractional-step/coupled pressure-velocity algorithms.

This ordering keeps architecture work separate from numerical changes and makes subsequent V&V failures attributable to the method under test.

## Non-negotiable evidence rules

- No tolerance relaxation to promote a method.
- No disabled test to promote a method.
- No analytical value injected into the solver path.
- True residuals and conservation are recomputed independently where applicable.
- Configuration must expose the selected method; hidden defaults are not accepted as a final contract.
- Performance and scaling are separate claims from correctness.


## Repository audit baseline (Issue #461)

The registry is intentionally based on the production implementation rather than on the desired future method list.

- **Gradients:** Green-Gauss (cell- and vertex-based) and least-squares are implemented. Least-squares is currently a field-level wrapper around the existing rank-aware kernel in `core/fvm/least_squares_gradient.h`. The wrapper must resolve the adjacent cell from whichever side of a face the cell lies on; using only faces the cell owns reduces the stencil to a one-sided, first-order operator. The vertex-based Green-Gauss interpolates cell values to vertices with inverse-distance weights, then constructs face values and the usual Green-Gauss sum; it is linear-exact on affine meshes but not quadratic-exact. `tests/validation/test_gradient_verification.cpp` exercises all three schemes on orthogonal, affine-sheared and stretched mesh families (constant, linear, quadratic and smooth fields, measured observed order), and pins the zero-gradient boundary policy (boundary face interpolation weight 0) and the boundary-cell approximation. `tests/unit/test_gradient.cpp` pins the least-squares face-orientation behaviour. `tests/validation/test_polyhedral_gradient_campaign.cpp` runs the schemes on a tetrahedral (Kuhn) grid with a 2:1 refinement sweep (n = 4/8/16): constant exact for all; linear exact for LS, vertex GG and point-linear GG; smooth-field observed order — two-point GG ~0.25/0.08 (not linear-consistent), vertex GG ~1.06/0.98, point-linear GG ~1.45/1.08, least-squares ~0.77/0.96, 2-ring LS ~1.09/1.01, and **quadratic least-squares (~2.13/2.01) is the scheme that reaches second order on tetrahedra**. The quadratic-scheme order is asserted (floor 1.5) in the campaign.

> **Investigation (2026): why cell Green-Gauss is not even linear-consistent on tetrahedra.** The two-point Gauss face value `(1-w)φ_owner + w φ_nb` evaluates the field at a point on the owner→neighbour centre line. The `geometric_interpolation_weight` therefore only yields `φ(Cf)` when the face centroid `Cf` lies on that line (affine meshes). For a tetrahedron the face centroid is off the line, so a *linear* field is not reproduced (the interpolation error is O(1) with respect to `∇φ` and does not vanish under refinement — measured observed order ≈ 0.1). Neither the magnitude-ratio nor the projection weight fixes this; the divergence-theorem identity `Σ_f Cf⊗Sf = V·I` holds exactly, so the geometry is consistent and the residual is entirely the interpolation bias. This is the documented skewness issue of [Greenshields & Weller, Notes on CFD (CFD Direct, 2022), §3.15 "Gradient discretisation"], which recommends the **point-linear (skew-corrected)** interpolation (face-point value + inverse-distance vertex values, area-weighted over the face triangles) or least-squares for bad/elongated-tet meshes; the OpenFOAM v11 User Guide echoes this (`gradSchemes: cellLimited Gauss linear`, `leastSquares`). A nested fixed-subregion measurement confirms even least-squares/vertex-GG reach only ~1st order at these resolutions on tets. The **point-linear (skew-corrected)** face reconstruction is implemented as `compute_gradient_gauss_point` (`gradient.gauss_point`) and is linear-exact on tetrahedra (~1st-order-converging, ~1.1-1.45 observed order); `gradient.least_squares_quadratic` (2-ring linear+Hessian fit) is the scheme that reaches **~2nd order on tetrahedra** (observed 2.13/2.01 with a 2:1 sweep, asserted floor 1.5). On affine meshes the two-point GG remains quadratic-exact/second-order and is the default; the point-linear/quadratic options are the linear-consistent choices for polyhedral-heavy meshes (this is the N1 explicit-selection model).
- **Interpolation/limiters:** linear and upwind interpolation are production paths. MinMod, Van Leer, Superbee, Van Albada and MC exist as limiter coefficient kernels. **Blended central/upwind**: `InterpScheme::BLENDED` (`numerics.convection.blended`) computes `φ_f = β·0.5(φ_P+φ_N) + (1−β)·upwind` for β ∈ [0,1]. On the 3-D Cartesian campaign `tests/validation/test_convection_blending.cpp` verifies: face-level boundedness for every β (the convex blend of two cell values stays in the adjacent envelope), β=0 is exactly upwind, and the smooth-profile observed order goes ~1.0 (β=0) → ~0.95 (β∈[0.25,0.75]) → **exactly 2.0 (β=1, central)**, so β is an accuracy/robustness knob at the face layer. Production `compute_convection` wiring is pending the open #514 (convection.h).
- **Convection:** the scalar FVM transport path exposes upwind, bounded second-order upwind, limited (TVD/MUSCL) and central. `compute_convection` now accepts a `LimiterType` and computes the cell gradients internally for the limited path, so every limiter is selectable through `numerics.convection.tvd.<limiter>`. `tests/validation/test_convection_scheme_verification.cpp` verifies on a uniform 1-D slab grid: boundedness (no new extrema) for all five limiters, linear exactness on deep interior faces, observed order (upwind ≈ 1, TVD ≈ 2), and discrete conservation. QUICK-equivalent, central/blended hybrid and multidimensional limiter policies remain open.
- **Diffusion:** orthogonal, corrected, limited, uncorrected and over-relaxed variants exist in the finite-volume path. The over-relaxed formulation uses `alpha = |Sf|^2/(Sf.d)` so the residual correction `Sf - alpha d` is orthogonal to `Sf`; it is rejected on a face where `Sf.d = 0` because the coefficient is singular. The non-orthogonal correction reads a cell gradient; `compute_laplacian` takes a selectable `GradientScheme` (`gradient.h`) defaulting to the two-point Gauss gradient (affine-optimal), with `GAUSS_POINT`, `LEAST_SQUARES` and `LEAST_SQUARES_QUADRATIC` as linear-consistent options. `tests/validation/test_nonorthogonal_laplacian_campaign.cpp` quantifies the operators on an affine skew ladder (shear 0 → 2): corrected/limited/over-relaxed are linear-exact and converge at order ≈ 2 on a smooth manufactured field, while the uncorrected two-point operator does not converge under skew. `tests/validation/test_polyhedral_laplacian_campaign.cpp` documents the tetrahedral behaviour (constant exact; the correction strictly reduces the error vs uncorrected, but operators remain O(1)/non-convergent on a smooth field because the underlying cell gradients are not accurate on tets — see the gradient polyhedral accuracy gap). Measured operators: on interior tetrahedra the corrected Laplacian annihilates a *linear* field only when a least-squares-family gradient (linear-consistent at every cell) is used (`Linf ≈ 5e-13`); with the default two-point Gauss gradient the operator is not even linear-consistent on tets. A smooth-field Laplacian on tets stays O(1) even with linear-consistent gradients, so restoring polyhedral diffusion accuracy is a deeper discretisation problem, not just a gradient swap. **Precise characterisation (2026):** with a genuinely second-order quadratic-LS gradient in the correction, the corrected/over-relaxed smooth Laplacian still **plateaus at L2 ≈ 2.1 (observed order ≈ 0) even in the deep interior** `[0.3,0.7]³` of the Kuhn grid — the face gradient `grad_f = ½(g_P+g_N)` is evaluated at the cell-centre midpoint, which is off the face centroid on any non-affine mesh; the correction residual then does not vanish under refinement (exact only for linear fields). **Refined measurement (Oct-2026):** even evaluating the face gradient AT the face centroid by first-order Hessian extrapolation (from the quadratic fit, exposed as `compute_quadratic_fit_least_squares`) leaves the divergence **O(1)** (corrected L2 ≈ 1.88/1.80/1.94, over-relaxed ≈ 2.09/2.04/2.20 at n=8/16/24, order ≈ 0): the two-point term `alpha(φ_N−φ_P)` carries an O(h²·H) curvature mismatch per face which accumulates to O(H) in the divergence and no gradient evaluation choice removes it. Closing the polyhedral diffusion gap therefore needs a different discretisation (e.g. a fully reconstructed face-normal gradient with quadratic terms), tracked as a research-grade item.

### Diffusion skewness-correction contract

For an internal face with centre-distance vector `d = C_N - C_P` and area vector `Sf`, the four non-orthogonal policies are explicit:

| scheme | implicit coefficient | correction residual | property |
|---|---|---|---|
| UNCORRECTED | `(Sf.d)/|d|^2` | — | two-point only; exact only on orthogonal-ish meshes |
| CORRECTED | `(Sf.d)/|d|^2` | `Sf - alpha d` | residual orthogonal to `d` |
| OVER_RELAXED | `|Sf|^2/(Sf.d)` | `Sf - alpha d` | residual **orthogonal to `Sf`** (minimum magnitude); singular at `Sf.d = 0` |
| LIMITED | `(Sf.d)/|d|^2` | `Sf - alpha d`, limiter in [0,1] | correction bounded by the orthogonal part |

This decomposition table is the explicit skewness-correction contract; the campaigns above are its quantitative evidence.
- **Temporal/CFL:** Euler explicit/implicit, Crank-Nicolson, BDF2, local time stepping, adaptive CFL and pseudo-transient CFL are present. Explicit RK2 (midpoint) and RK3 (Williamson low-storage) are provided by `cfdx/physics/low_storage_time_integration.h` and registered (`numerics.temporal.rk2`/`.rk3`). `tests/validation/test_temporal_order_matrix.cpp` asserts the measured temporal order matrix on u'=-u with a 2:1 sweep (Euler~1, CN~2, BDF2~2, RK2~2, RK3~3). NOTE: the RK2 helper was previously implemented with `dt/2` weights on BOTH stages (violating `b·c = 1/2`), which made it only first order; it is fixed to the true midpoint rule (see #461 N5). **Variable-step BDF2**: `advance_time`/`TimeIntegrationContext` now carry `dt_prev`; the implicit solve uses the `bdf2_coefficients(dt_prev, dt)` ω-ratio coefficients (reducing to the constant-step formula at ω=1), and `tests/unit/test_temporal_history.cpp` verifies the coefficients, a single analytic variable-step step, the history lifecycle (`initialize`/`shift`/`has_prev`), ~2nd order on a scaled variable-step pattern, and the restart semantics (a fresh context re-seeded via `initialize(φ)` degrades only the first post-restart step to an implicit-Euler bootstrap, then BDF2 resumes; the restart/continuous gap is O(dt)).
- **Linear algebra:** CG, BiCGStab, GMRES and FGMRES are available; native AMG/Smoothed Aggregation AMG and block FieldSplit/Schur infrastructure are present. Planned catalog entries are not registered as implemented methods.
- **Pressure-velocity:** SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and coupled paths are present. Their registration remains `implemented` until method-specific quantitative evidence is attached.

This audit rule prevents capability inflation: a reusable numerical primitive is not promoted to a production scheme merely because its formula exists.
