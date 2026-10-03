# Numerical Methods Maturity Audit — Issue #461

## Purpose

Repository-level audit of the numerical-method maturity work tracked by #461.

Status semantics:
- implemented = code path exists;
- verified = executable verification evidence exists;
- validated = quantitative physical/reference validation exists;
- qualified = all applicable gates are satisfied;
- partial = some required elements exist, but the package is incomplete;
- missing = no sufficient implementation/evidence found;
- planned = intentionally deferred.

An implementation registry entry never implies validation.

## Executive result

| Package | Status | Main evidence / gap |
|---|---|---|
| N1 contracts | PARTIAL | registry, explicit resolver + case-level resolution report (`resolve_case_numerics` / `format_numerics_report`, deterministic) exist; wiring every case loader and the common-scheme/contracts items remain |
| N2 gradients/reconstruction | QUALIFIED* | production-qualified point-linear Green–Gauss + LS/WLS/quadratic-LS families with stencil-quality/rank/conditioning diagnostics, explicit Dirichlet/Neumann/mixed boundary reconstruction, and orthogonal/skewed/stretched/distorted/non-affine tetrahedral V&V; quadratic LS has an executable ~2nd-order gate on the non-affine family. Two-point cell Green–Gauss is explicitly restricted to its affine/centre-line-consistent domain; vertex GG remains a secondary method with a narrower verified domain. Limiter qualification belongs to N4, not N2. *Promotion after exact CI run for the new non-affine campaign. |
| N3 diffusion | QUALIFIED | orthogonal/corrected/limited/over-relaxed diffusion has executable unit/V&V coverage, explicit skewness-correction contract, quantitative affine non-orthogonality ladders, singular-face rejection, and independent mesh-quality diagnostics (mean/max non-orthogonality, skewness and threshold counts). The qualification scope is the N3 acceptance contract in #461: controlled orthogonal/non-orthogonal diffusion, quantified degradation/robustness and no silent scheme change. General polyhedral higher-order diffusion remains a separately tracked research item under N10/N16 and is not claimed by N3. |
| N4 convection | PARTIAL | production wiring now covers upwind, SOU, MUSCL/TVD with all 5 limiters, central, explicit blending, QUICK and bounded QUICK; N2 gradient selection remains independent; 3-D Cartesian boundedness/order/conservation and production flux reconstruction are covered. Remaining qualification is general-mesh V&V (skew/non-orthogonal/tetra/polyhedral), steep-gradient/discontinuity matrix and final measured scheme-level gates |
| N5 temporal | PARTIAL | measured order matrix verified (Euler~1, CN~2, BDF2~2, RK2~2, RK3~3), variable-step BDF2 + history lifecycle + restart semantics (O(dt) gap) implemented, and a transient manufactured run on real operators delivering spatial order 2.0 (spectral heat) and mesh-locked temporal order 2.0 (CN/BDF2, advection); remaining: 3-D transient MMS and a stiff Newton-based temporal sweep |
| N6 CFL/time-step | PARTIAL | deterministic global/local convective CFL, geometry-aware h=2V/A characteristic length, bounded adaptive controller, auditable decision history, production nonlinear rollback/retry, bounded stagnant-cell pseudo-time policy, dual-time stepping + adaptive field integration, and a transactional physical-time retry/RESTART lifecycle (self-start on rejection, deterministic restart) are implemented and unit-verified; production PDE coupling of the lifecycle and a Newton-based stiff implicit solve remain |
| N7 nonlinear convergence | PARTIAL | convergence infrastructure exists; unified conservation/QoI/continuation qualification incomplete |
| N8 linear/preconditioners | PARTIAL | AMG/FieldSplit/Schur advanced; exact-Schur oracle + SIMPLE/SIMPLEC algebraic Schur (vs-oracle error qualified); LSC/BFBt/PCD, null-space and scaling remain |
| N9 pressure-velocity | PARTIAL | multiple algorithms are registered; common quantitative acceptance matrix is incomplete |
| N10 difficult meshes | PARTIAL | #423 supplies mesh robustness infrastructure; numerical degradation campaign is not fully integrated |
| N11 conservation/boundedness | PARTIAL | #479/#480 provide strong diagnostics; full benchmark campaign remains |
| N12 MMS | PARTIAL | scalar diffusion and analytical campaigns exist; reusable all-equation MMS laboratory is incomplete |
| N13 cross-backend | IMPLEMENTED — PENDING EXECUTION | Equivalence contract, serial/MPI campaign, assembled/matrix-free operator campaign, deterministic/normal reduction campaign, restart campaign, N-to-M evidence wiring and explicit CUDA BLOCKED semantics are implemented; final status remains execution-dependent and CUDA requires real hardware |
| N14 reference matrix | IMPLEMENTED | Machine-readable method matrix, five reference dossiers, explicit comparison categories and benchmark cross-reference are now present; quantitative cross-code campaigns remain external evidence work |
| N15 performance | PARTIAL | reproducible benchmark runner, result schema, provenance, strong/weak scaling orchestration and regression analysis added; real MPI/GPU/AMG/matrix-free/mixed-precision campaigns still require execution on appropriate hardware and solver instrumentation |
| N16 advanced methods | PARTIAL | N16.1 adds a verified WENO5-Z finite-volume reconstruction primitive (standard candidate polynomials, WENO-Z smoothness indicators and nonlinear weights, boundedness policy, constant/linear exactness and a smooth refinement campaign on a uniform 1-D cell-average stencil). It is an isolated kernel: it is not wired into any production scheme, and multidimensional/polyhedral higher-order reconstruction remains open |
| N17 turbulence | PARTIAL | 13-model registry/qualification work exists; physical qualification is incomplete |

## Critical findings

### Registry status is not qualification status

The existing numerical-method registry and capability matrix contain useful implementation metadata.
They must not be interpreted as proof of verification, validation or qualification.

The turbulence work in #473/#474 already uses the stronger separation:
equation verification, solver verification and physical validation.

The same discipline must be applied to the complete #461 matrix.

### N2 qualification boundary

N2 is considered closed only for gradient/reconstruction methods whose applicability domain is explicitly stated and covered by executable evidence. The two-point cell-centred Green–Gauss operator is **not** a general polyhedral gradient: on tetrahedra the owner/neighbour centre-line interpolation is not face-centroid consistent and is therefore retained as an implemented, domain-restricted method rather than being artificially promoted. The point-linear Green–Gauss operator is the polyhedral Green–Gauss path used for skew/non-affine meshes.

Gradient limiting is deliberately not a remaining N2 gate. Limiters are part of the N4 convection/reconstruction qualification and must be qualified at scheme level, not as a gradient-kernel property.

The new non-affine tetrahedral campaign perturbs the interior geometry with a smooth boundary-preserving deformation whose amplitude scales with h. It checks point-linear Green–Gauss, extended WLS and quadratic LS under non-affine geometry, with a strict second-order gate for quadratic LS.

### Limiter kernels are not complete convection schemes

MinMod, Van Leer, Van Albada, Superbee and MC coefficient implementations exist.
This is useful unit-level coverage, but it does not prove that every limiter is a complete
production convection scheme with face reconstruction, conservation, multidimensional
boundedness, smooth-field order, steep-gradient behaviour and skew/non-orthogonal robustness.

### Schur terminology must remain exact

The current coupled Schur implementation from #482 uses a block-local approximation:

S_tilde = C - D M_b^{-1} G

It must not be described as the exact Schur:

S = C - D A_uu^{-1} G

#481 already identifies this distinction. A small-system exact algebraic oracle remains useful.

### Pressure-velocity algorithm presence is not V&V

SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and coupled U-p paths are represented,
but the common Couette/Poiseuille/cavity/external/skew-mesh acceptance matrix is not complete
for every algorithm.

### Conservation instrumentation is ahead of qualification

#479/#480 now provide the correct direction: conservation and boundedness are independently
recomputed instead of inferred from solver residuals. The remaining work is the benchmark
campaign and quantitative acceptance evidence.

### Cross-backend support is not equivalence

MPI, matrix-free, CUDA, mixed precision and restart foundations exist. They need separate
functional, numerical-equivalence, performance and memory evidence.

## Remediation priorities

### P0 — eliminate false-positive maturity signals
- Keep implemented separate from verified and validated.
- Add evidence links to promoted registry entries.
- Keep exact/approximate Schur terminology explicit.
- Do not treat limiter coefficient tests as scheme qualification.

### P1 — complete numerical V&V
- gradient/reconstruction matrix;
- convection scheme matrix;
- non-orthogonal diffusion campaign;
- temporal refinement/order campaign;
- convergence/conservation/QoI gates;
- pressure-velocity benchmark matrix;
- skew/high-aspect-ratio/polyhedral campaign.

### P2 — scalable equivalence
- serial/MPI;
- assembled/matrix-free;
- CPU/CUDA;
- restart and N-to-M MPI;
- reproducible performance and memory measurements.

### P3 — turbulence physical qualification
Use #473/#474 as the authority for 13-model qualification semantics and keep equation,
solver and physical validation separate.

## Repository evidence inventory

Key sources audited:
- src/cfdx/core/numerics/numerical_method_contract.h
- src/cfdx/core/numerics/numerical_method_registry.h
- src/cfdx/core/numerics/gradient.h
- src/cfdx/core/fvm/least_squares_gradient.h
- src/cfdx/core/numerics/interpolation.h
- src/cfdx/core/numerics/convection.h
- src/cfdx/core/numerics/convection_assembly.h
- src/cfdx/core/numerics/laplacian.h
- src/cfdx/core/numerics/temporal.h
- src/cfdx/physics/temporal.h
- src/cfdx/physics/advanced_convergence.h
- src/cfdx/core/numerics/conservation.h
- docs/validation/CONSERVATION_BOUNDEDNESS.md
- docs/validation/NUMERICAL_METHOD_CAPABILITY_MATRIX.json
- docs/validation/NUMERICAL_MODEL_VERIFICATION_MATRIX.md
- docs/validation/CFDX_QUALIFICATION_CAMPAIGN.md
- docs/validation/CFDX_QUALIFICATION_REGISTRY.json
- docs/validation/TURBULENCE_QUALIFICATION_MATRIX.json
- docs/development/NUMERICAL_METHOD_CONTRACTS.md
- tests/validation/test_numerical_model_verification.cpp
- tests/validation/test_discretization_verification.cpp
- tests/validation/test_benchmark_matrix.cpp
- tests/unit/test_conservation_boundedness.cpp

Related work audited for overlap:
#60, #405, #423, #473, #479, #480, #481, #482, #484.

## Scope of this PR

This PR does not claim that #461 is complete.

It adds a durable, machine-checkable statement of what is currently implemented versus
verified/validated and makes the remaining gaps explicit. Future implementation PRs should
update this evidence rather than simply checking boxes in #461.
