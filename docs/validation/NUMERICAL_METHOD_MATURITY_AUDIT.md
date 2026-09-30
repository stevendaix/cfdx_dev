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
| N1 contracts | PARTIAL | numerical contracts, registry and explicit selection resolver exist; complete case/report traceability is incomplete |
| N2 gradients/reconstruction | PARTIAL | cell/vertex/point-linear Green-Gauss + linear/quadratic least-squares with orthogonal/skewed/stretched observed-order V&V, pinned zero-gradient boundary policy, and a tetrahedral campaign (constant exact all, linear exact LS/vertex/point; quadratic LS most accurate on tets ~0.02 L2 at n=16 but ~0.8-1.2 order; genuinely second-order polyhedral gradients still research-grade) |
| N3 diffusion | PARTIAL | orthogonal/corrected/limited/over-relaxed with explicit skewness-correction contract and quantified skew ladder (order ~2; uncorrected non-convergent under skew); tetrahedra: corrected operator linear-exact only with LS-family gradients (Linf ~5e-13) but smooth Laplacian remains O(1) — a discretisation gap, not just a gradient swap |
| N4 convection | PARTIAL | selectable upwind/TVD (MUSCL, 5 limiters) + central with 1D boundedness/linear-exactness/order/conservation V&V; QUICK, blended, multidimensional policy and production wiring incomplete |
| N5 temporal | PARTIAL | Euler/CN/BDF2 exist; history/restart/RK/variable-step/order campaign incomplete |
| N6 CFL/time-step | PARTIAL | adaptive CFL/pseudo-transient controllers exist; full characteristic-length/rollback/audit contract incomplete |
| N7 nonlinear convergence | PARTIAL | convergence infrastructure exists; unified conservation/QoI/continuation qualification incomplete |
| N8 linear/preconditioners | PARTIAL | AMG/FieldSplit/Schur are advanced; exact-Schur oracle, LSC/BFBt/PCD and scaling remain |
| N9 pressure-velocity | PARTIAL | multiple algorithms are registered; common quantitative acceptance matrix is incomplete |
| N10 difficult meshes | PARTIAL | #423 supplies mesh robustness infrastructure; numerical degradation campaign is not fully integrated |
| N11 conservation/boundedness | PARTIAL | #479/#480 provide strong diagnostics; full benchmark campaign remains |
| N12 MMS | PARTIAL | scalar diffusion and analytical campaigns exist; reusable all-equation MMS laboratory is incomplete |
| N13 cross-backend | PARTIAL | MPI/matrix-free/GPU foundations exist; equivalence campaign is incomplete |
| N14 reference matrix | PARTIAL | Fluent/VMFL and solver-comparison documents exist; method-level evidence matrix is incomplete |
| N15 performance | MISSING | scattered benchmarks exist, but no complete reproducible maturity campaign |
| N16 advanced methods | PLANNED | intentionally deferred by #461 |
| N17 turbulence | PARTIAL | 13-model registry/qualification work exists; physical qualification is incomplete |

## Critical findings

### Registry status is not qualification status

The existing numerical-method registry and capability matrix contain useful implementation metadata.
They must not be interpreted as proof of verification, validation or qualification.

The turbulence work in #473/#474 already uses the stronger separation:
equation verification, solver verification and physical validation.

The same discipline must be applied to the complete #461 matrix.

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
