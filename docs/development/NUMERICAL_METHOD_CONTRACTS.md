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

- **Gradients:** Green-Gauss and least-squares are implemented. Least-squares is currently a field-level wrapper around the existing rank-aware kernel in `core/fvm/least_squares_gradient.h`. The wrapper must resolve the adjacent cell from whichever side of a face the cell lies on; using only faces the cell owns reduces the stencil to a one-sided, first-order operator. `tests/validation/test_gradient_verification.cpp` now exercises both gradients on orthogonal, affine-sheared and stretched mesh families (constant, linear, quadratic and smooth fields, measured observed order), and `tests/unit/test_gradient.cpp` pins the face-orientation behaviour.
- **Interpolation/limiters:** linear and upwind interpolation are production paths. MinMod, Van Leer, Superbee, Van Albada and MC exist as limiter coefficient kernels; their presence does **not** imply that all five are selectable production convection schemes.
- **Convection:** the scalar FVM transport path currently exposes upwind, bounded second-order upwind and TVD; the production TVD branch is currently wired to **MinMod**. Other limiter kernels remain implementation primitives until a complete transport-path integration and V&V campaign exists.
- **Diffusion:** orthogonal, corrected, limited, uncorrected and over-relaxed variants exist in the finite-volume path. The over-relaxed formulation uses `alpha = |Sf|^2/(Sf.d)` so the residual correction `Sf - alpha d` is orthogonal to `Sf`; it is rejected on a face where `Sf.d = 0` because the coefficient is singular.
- **Temporal/CFL:** Euler explicit/implicit, Crank-Nicolson, BDF2, local time stepping, adaptive CFL and pseudo-transient CFL are present.
- **Linear algebra:** CG, BiCGStab, GMRES and FGMRES are available; native AMG/Smoothed Aggregation AMG and block FieldSplit/Schur infrastructure are present. Planned catalog entries are not registered as implemented methods.
- **Pressure-velocity:** SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and coupled paths are present. Their registration remains `implemented` until method-specific quantitative evidence is attached.

This audit rule prevents capability inflation: a reusable numerical primitive is not promoted to a production scheme merely because its formula exists.
