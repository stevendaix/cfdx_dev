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
