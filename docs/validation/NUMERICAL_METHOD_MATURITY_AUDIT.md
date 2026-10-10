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
| N7 nonlinear convergence | PARTIAL | convergence monitor + stagnation/divergence detection + bounded deterministic adaptive relaxation + continuation/homotopy are implemented with independent residual/continuity/QoI acceptance gates, deterministic reports and a Re=100 cavity qualification corridor. Provided/uniform/restart initialization modes are explicitly validated. Remaining gaps are automatic initialization selection, a canonical continuation failure/recovery corridor, broader geometry/mesh/parameter qualification, spatial observed-order evidence, which is merge-gated by a 16/32/64 tripwire (test_ghia_cavity_order, order U_RMS=2.05, V_RMS=1.16) while the full 32/64/128 plus Re=400/Re=1000 campaign now runs weekly on schedule and establishes neither the spatial accuracy of the continuation or adaptive paths nor a general order on other schemes, meshes or geometries |
| N8 linear/preconditioners | PARTIAL | AMG/FieldSplit/Schur infrastructure is implemented; exact-Schur oracle, SIMPLE/SIMPLEC, LSC/BFBt and PCD have repository implementations and dedicated evidence. The common coupled production path now exposes BlockLocal, PCD, SIMPLE, SIMPLEC, LSC and BFBT with explicit requested/resolved model diagnostics and lifecycle/CSR-graph guards. Remaining N8 work is package-level qualification: representative CFD Schur acceptance envelopes, broader pressure null-space evidence, difficult-matrix AMG evidence and exact-HEAD qualification governance. Physical pressure-velocity qualification belongs to N9 |
| N9 pressure-velocity | PARTIAL | multiple algorithms are registered; common quantitative acceptance matrix is incomplete |
| N10 difficult meshes | PARTIAL | #423 public fixtures are now independently acquired/verified and exercised through production import/topology/geometry-quality validation; the N10 fixture acceptance matrix is reconciled by `scripts/audit_n10_fixture_evidence.py`. Numerical solver-support evidence is explicitly limited to the applicable `meshio-vtk-unstructured` execution/convergence smoke path; package promotion still requires retained quantitative evidence from the exact CI HEAD and broader fixture-specific accuracy/physical-validation evidence |
| N11 conservation/boundedness | PARTIAL | #479/#480 and merged PRs #767/#776/#780/#787/#814 provide conservation/boundedness diagnostics, momentum/scalar/energy qualification slices, quantitative campaign infrastructure and an independent MPI Poisson flux-balance gate. N11 remains partial: complete Couette/Poiseuille/scalar-MMS/energy/skew acceptance matrix, broad independent production-flux reconstruction and per-variable positivity/clipping telemetry are not yet qualified. PR #814 covers the active distributed Poisson path only, not production incompressible/thermal MPI |
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


## 2026-10-04 — N8 tracking re-baselining

This entry corrects N8 statements that the merged evidence had already overtaken, and records gaps that no layer stated. It changes no tolerance, threshold, gate or solver behaviour, and it claims no new qualification.

**Corrected as stale.**

- *"LSC/BFBt/PCD Schur approximations — missing"* in `NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` was wrong on two of its three subjects. LSC and BFBt are implemented behind the common `SchurApproximation` interface (`src/cfdx/core/linalg/lsc_bfbt_schur.h`, merged in #562) and are covered by four unit files. Only PCD is missing. The requirement is split accordingly, and LSC/BFBt is recorded as `implemented` — the same status the audit already gives SIMPLEC — not as qualified: the #568 harness keeps approximation quality diagnostic because no acceptance envelope exists yet.
- *"Resolved automatic linear plan is not yet recorded in the campaign report"* was overtaken by #620. `scripts/n8_physical_qualification.py` now carries `audit_linear_plan`, whose `explicit_linear_request_must_be_honored` policy is the one gate in the report that can fail the campaign.

**Recorded as new gaps.**

- The generic Schur layer is test-only. No file under `src/` or `apps/` includes `schur_approximation.h`, `exact_schur.h`, `lsc_bfbt_schur.h` or `simplerc_schur.h`, and no case-file or solver-options setting selects a Schur approximation. Production selects the separate `CoupledBlockSchurAMGPreconditioner`/`CoupledBlockSchurPreconditioner` pair, which approximates the velocity block cell-locally rather than by `Auu^-1`. The #481 definition of done — SIMPLEC/LSC/BFBt/PCD interchangeable in production — is therefore not met by the current tree, and this was not previously stated in any of the three layers.
- The Schur layer provides only the constant (mean-zero) pressure null-space projector. Pinned pressure in production is handled by reference-cell gauge elimination, a different mechanism the Schur layer does not use.

**Closed after this entry was written.** #634 gives `SimplerSchurApproximation` the same two contracts the exact and LSC/BFBt approximations already had: `update_values` rejects a CSR graph change, and an optional pressure `NullSpaceProjector` is honoured (opt-in, so existing callers keep their algebra). The two gaps this entry listed for SIMPLE/SIMPLEC are therefore addressed by #634 rather than outstanding; the pinned-pressure policy is not.

**Status unchanged.** N8 stays PARTIAL. Nothing here promotes a package or requirement to `qualified`; the outstanding work is production integration of the Schur layer, PCD, the null-space and lifecycle gaps above, dispatcher-fed matrix characteristics, the anisotropic AMG sweep, and the measured Schur acceptance envelope.


## 2026-10-05 — N9 Schur-family source-of-truth re-baselining

PR #660 is now merged as the N9.6 qualification-matrix contract. This entry re-baselines the Schur-family audit against the exact merged tree and supersedes the older N8 statements that pre-date the PCD implementation.

### Current repository evidence

- **Exact Schur** is implemented as the reference/oracle in `src/cfdx/core/linalg/exact_schur.h`, with dedicated N9 algebraic qualification evidence.
- **SIMPLE/SIMPLEC** are implemented through `src/cfdx/core/linalg/simplerc_schur.h` and its tests. The lifecycle/null-space contracts added by #634 are retained.
- **LSC/BFBt** are implemented through `src/cfdx/core/linalg/lsc_bfbt_schur.h`, with unit, null-space, oracle-comparison and quantitative tests.
- **PCD** is implemented through `src/cfdx/core/linalg/pcd_schur.h`, has explicit setup/update-values lifecycle checks, and is wired into `CoupledBlockSchurAMGPreconditioner` through `CoupledSchurApproximationModel::PCD`.
- **MGR** remains an N8 coupled-preconditioning path and is to be physically qualified by N9; its implementation is not duplicated here.
- **SIMPLER, Yosida, Cahouet–Chabard, Augmented-Lagrangian Schur and Inexact Uzawa** are not evidenced as distinct production implementations on the current tree. They remain planned/under assessment and must not be inferred from similarly named infrastructure.

### N8/N9 boundary

N8 owns reusable linear-algebra and preconditioner infrastructure. N9 owns pressure–velocity algorithm integration and common physical qualification. Therefore the presence of PCD/LSC/BFBt classes does **not** by itself close N8 or N9.

The current N9.6 matrix in `docs/validation/N9_PRESSURE_VELOCITY_QUALIFICATION_MATRIX.md` remains a qualification contract. It deliberately does not promote unexecuted cells.

### Required next work

1. Make the common `SchurApproximation` family production-selectable without silent model substitution.
2. Expose requested/resolved Schur model diagnostics.
3. Complete production PCD/LSC/BFBt selection and lifecycle coverage.
4. Build representative non-zero-flow/Oseen qualification evidence.
5. Keep Exact Schur as a controlled oracle/reference only.
6. Classify advanced methods before implementation; move methods that belong to N8 or a later advanced-solvers package out of N9.

No tolerance, iteration limit, acceptance gate or fallback policy is relaxed by this audit update.


## 2026-10-07 — N10 public-fixture acceptance-matrix synchronization

The N10 public-fixture chain is now reflected in the package-level maturity audit without promoting N10 to qualified. The four pinned public fixtures are independently verified and passed through the production importer with explicit 2-D/3-D topology and geometry-quality handling. The authoritative acceptance gate reconciles the independent verification, production qualification and persisted applicability fixture sets and rejects duplicate identifiers.

The current numerical-support boundary remains explicit: only `meshio-vtk-unstructured` has production execution/convergence smoke evidence and fixture-specific numerical reconstruction evidence. SU2, Gmsh and OpenFOAM fixtures remain solver applicability **NOT_CLAIMED** rather than being forced through an unsupported 3-D path. Physical validation and broad fixture-specific accuracy are also **NOT_CLAIMED**.

This synchronization is evidence governance only. It changes no numerical tolerance, convergence criterion, solver fallback or validation gate. N10 remains **PARTIAL** until the complete N10 quantitative evidence set is retained from the exact tested HEAD and the declared acceptance matrix is sufficient for package-level promotion.


## 2026-10-07 — N8 production Schur integration synchronization

The previous N8 entries stating that SIMPLE/SIMPLEC/LSC/BFBT were test-only are stale. Issue #723 has been completed after PR #804 centralized production coupled-block extraction and an audit confirmed explicit Schur-family selection, requested/resolved diagnostics, lifecycle/graph guards and production regression coverage. Exact Schur remains an oracle/reference and is not an automatic production fallback. MatrixCharacteristics was already completed by #790 and is not a remaining N8 item. This synchronization changes documentation only; it does not promote N8 from PARTIAL to QUALIFIED.


## N11 re-audit — 2026-10-10

Baseline: `master` at `e0026f4ae3934929290ef727d152415567fd284c` (merge #837). PRs #767, #776, #780, #787 and #814 are merged. The N11 package remains **PARTIAL**; implementation and selected executable evidence are not equivalent to package-level qualification.

- **Conservation evidence present:** momentum gate (#767), scalar conservation/boundedness (#776), energy/CHT conservation (#780), quantitative campaign evidence (#787), and independent post-solve flux reconstruction plus per-cell balance on the active distributed Poisson path (#814).
- **MPI boundary:** #814 does not qualify production incompressible/thermal MPI; those production physics paths are not distributed under MPI.
- **Remaining acceptance gaps:** complete quantitative Couette/Poiseuille/scalar-MMS/energy/skew campaign; independent flux reconstruction across all applicable production equations; complete global/local mass and boundary sign/orientation matrix; per-variable admissible bounds and positivity where required; explicit floor/clipping telemetry; and retained exact-HEAD campaign records.
- **Governance:** issue #480 was reopened and its checklist reconciled to the merged evidence and remaining gaps. No tolerances or acceptance gates were weakened.
