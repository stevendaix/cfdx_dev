# Numerical Methods Maturity Audit — Issue #461

> V&V status semantics and promotion gates are governed by `docs/validation/CFDX_VV_GOVERNANCE.md`. This audit interprets repository evidence against that governance; it does not promote a package by implementation presence alone.

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

The audit uses the following evidence chain: requirement → implementation → executable test → code verification → solution verification → physical validation (where applicable) → qualification.

| Package | Status | Main evidence / gap |
|---|---|---|
| N1 contracts | PARTIAL | registry, explicit resolver + case-level resolution report (`resolve_case_numerics` / `format_numerics_report`, deterministic) exist; the remaining item is operational confirmation of the exact current HEAD in the dedicated N1 CI gate |
| N2 gradients/reconstruction | PARTIAL | point-linear Green–Gauss + LS/WLS/quadratic-LS families have substantial executable verification, including stencil-quality/rank/conditioning diagnostics, explicit boundary reconstruction and non-affine tetrahedral coverage; quadratic LS reaches the expected second order on the controlled non-affine family. However, the package is not yet qualified because independent face reconstruction and the complete polyhedral/skew/non-orthogonal acceptance matrix remain incomplete in the requirement audit. Two-point cell Green–Gauss is explicitly restricted to its affine/centre-line-consistent domain; vertex GG remains a secondary method with a narrower verified domain. Limiter qualification belongs to N4, not N2. |
| N3 diffusion | QUALIFIED | orthogonal/corrected/limited/over-relaxed diffusion has executable unit/V&V coverage, explicit skewness-correction contract, quantitative affine non-orthogonality ladders, singular-face rejection, and independent mesh-quality diagnostics (mean/max non-orthogonality, skewness and threshold counts). The qualification scope is the N3 acceptance contract in #461: controlled orthogonal/non-orthogonal diffusion, quantified degradation/robustness and no silent scheme change. General polyhedral higher-order diffusion remains a separately tracked research item under N10/N16 and is not claimed by N3. |
| N4 convection | PARTIAL | production wiring now covers upwind, SOU, MUSCL/TVD with all 5 limiters, central, explicit blending, QUICK and bounded QUICK; N2 gradient selection remains independent; 3-D Cartesian boundedness/order/conservation and production flux reconstruction are covered. Remaining qualification is general-mesh V&V (skew/non-orthogonal/tetrahedral/polyhedral), the steep-gradient/discontinuity matrix and final measured scheme-level gates |
| N5 temporal | PARTIAL | measured order matrix verified (Euler~1, CN~2, BDF2~2, RK2~2, RK3~3), variable-step BDF2 + history lifecycle + restart semantics implemented, and production transient MMS closure includes 3-D spatial order plus stiff temporal order via exact linear solves. Remaining limitation is production `advance_time` implicit-Picard robustness for stiff `dt >> h` cases, plus broader production transient qualification |
| N6 CFL/time-step | PARTIAL | deterministic global/local convective CFL, geometry-aware h=2V/A characteristic length, bounded adaptive controller, auditable decision history, production nonlinear rollback/retry, bounded stagnant-cell pseudo-time policy, dual-time stepping + adaptive field integration, and the transactional physical-time retry/RESTART lifecycle are implemented and unit-verified. Remaining work is broader production qualification of the time-step/CFL controls; Newton-based stiff implicit solving belongs to the solver-level work rather than an N6 timestep-control gap |
| N7 nonlinear convergence | PARTIAL | convergence monitor + stagnation/divergence detection + bounded deterministic adaptive relaxation + continuation/homotopy are implemented with independent residual/continuity/QoI acceptance gates, deterministic reports and a Re=100 cavity qualification corridor. Provided/uniform/restart initialization modes are explicitly validated. Remaining gaps are automatic initialization selection, a canonical continuation failure/recovery corridor, broader geometry/mesh/parameter qualification, and the fact that Ghia refinement is long-validation-only rather than merge-gating |
| N8 linear/preconditioners | PARTIAL | AMG/FieldSplit/Schur advanced; exact-Schur oracle + SIMPLE/SIMPLEC algebraic Schur (vs-oracle error qualified); LSC/BFBt/PCD, null-space and scaling remain |
| N9 pressure-velocity | PARTIAL | multiple algorithms are registered; common quantitative acceptance matrix is incomplete |
| N10 difficult meshes | PARTIAL | #423 supplies mesh robustness infrastructure; numerical degradation campaign is not fully integrated |
| N11 conservation/boundedness | PARTIAL | #479/#480 provide strong diagnostics; full benchmark campaign remains |
| N12 MMS | PARTIAL | scalar diffusion and analytical campaigns exist; reusable all-equation MMS laboratory is incomplete |
| N13 cross-backend | PARTIAL | equivalence/campaign infrastructure exists, including serial/MPI, assembled/matrix-free, restart and explicit CUDA-BLOCKED semantics; final status remains execution-dependent because numerical-equivalence evidence is not yet fully executed and retained |
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

N2 is **not yet closed**. The executable campaigns now establish useful method-specific domains and a second-order quadratic-LS result on the controlled non-affine tetrahedral family, but the package-level acceptance still requires independent face reconstruction and a complete controlled polyhedral/skew/non-orthogonal matrix. A successful subset of gradient tests must not be promoted to package qualification. The two-point cell-centred Green–Gauss operator is **not** a general polyhedral gradient: on tetrahedra the owner/neighbour centre-line interpolation is not face-centroid consistent and is therefore retained as an implemented, domain-restricted method rather than being artificially promoted. The point-linear Green–Gauss operator is the polyhedral Green–Gauss path used for skew/non-affine meshes.

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
- docs/validation/03-conservation/chapter.py
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


## 2026-10-04 — Source-of-truth consistency audit

PR #623 corrected the principal maturity-status contradictions. This follow-up keeps the three authoritative layers aligned: `NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` records requirement-level evidence, `NUMERICAL_METHOD_MATURITY_AUDIT.json` records package-level maturity, and this report records the human-readable audit. N2 does not own gradient-limiter qualification; limiters belong to N4. N5/N6 gaps describe current solver/qualification boundaries rather than already-closed temporal or dual-time work. N7 reflects the merged #617 state without claiming broader qualification than the evidence supports. N8 remains PARTIAL pending the open production-selection/AMG/Schur qualification work.
