# CFDX Verification & Validation Governance

## Purpose
This document is the authoritative V&V governance layer for CFDX numerical-method maturity.

It defines what evidence is required, how evidence is classified, what constitutes a gate, and how evidence is traced to a CFDX requirement. It does not replace individual qualification cases, executable tests, or issue trackers.

The framework separates:
1. Code verification — is the numerical implementation correct?
2. Solution verification — is the numerical solution sufficiently converged and is discretisation error controlled?
3. Validation — does CFDX reproduce an independent physical/reference observation under matched conditions?
4. Uncertainty and limitations — what numerical, reference, and model-form limitations remain?

The methodology is aligned conceptually with ASME V&V 20, MMS-based code-verification practice, and Richardson/GCI-style grid-convergence analysis. These references provide methodology; they are not automatic CFDX acceptance oracles.

## 1. Evidence hierarchy
| Evidence | Question answered | Typical evidence |
|---|---|---|
| Implementation | Does the code path exist? | source, registry, unit test |
| Code verification | Is the numerical implementation consistent with its mathematical contract? | exact solution, MMS, polynomial reproduction, conservation identities |
| Solution verification | Is the computed solution numerically resolved? | residuals, conservation, mesh/time refinement, QoI stabilization |
| Validation | Does the physical prediction agree with an independent reference? | experiment, DNS, published benchmark, independent solver |
| Qualification | Is the declared use case supported by all applicable gates? | combined evidence record and reproducible CI campaign |

Implementation evidence cannot substitute for verification. Verification evidence cannot substitute for physical validation.

## 2. Status vocabulary
- IMPLEMENTED: executable code path exists.
- VERIFIED: applicable mathematical/code-verification evidence exists.
- VALIDATED: applicable independent physical/reference comparison exists.
- QUALIFIED: all mandatory gates for the declared scope are satisfied.
- PARTIAL: some required evidence exists but the scope is incomplete.
- DIAGNOSTIC: executable evidence exists but is intentionally not a qualification gate.
- READY: executable infrastructure exists and the case can be run; no PASS is implied.
- BLOCKED: evidence cannot currently be produced because of an identified prerequisite.
- PLANNED: intentionally not implemented or tested yet.

A registry entry must never infer QUALIFIED from IMPLEMENTED.

## 3. Verification classes
### 3.1 Code verification
Code verification targets the mathematical implementation. Preferred evidence includes constant-field preservation, linear-field exactness where applicable, polynomial reproduction, manufactured solutions, discrete conservation identities, independent operator checks, and limiting/positivity invariants where part of the contract.

For an observed-order claim, the exact or manufactured field, mesh family, refinement ratio, norm, and fitting procedure must be frozen.

A unit test such as a limiter-coefficient equality is component verification, not complete scheme verification.

### 3.2 Solution verification
Solution verification concerns a complete numerical solution. Applicable campaigns should report nonlinear residuals, true/recomputed residuals where available, global conservation imbalance, relevant local conservation diagnostics, QoI stabilization, spatial refinement, temporal refinement for transient problems, and solver-tolerance sensitivity where relevant.

A solution is not converged merely because a maximum iteration count was reached.

### 3.3 Validation
Validation requires an independent reference with documented provenance and matched geometry, physical parameters, boundary/initial conditions, dimensionality, nondimensional numbers, QoI definition, and sampling convention.

Reference uncertainty must be reported when available. Agreement with a reference cannot compensate for an unresolved numerical-verification failure.

## 4. Qualification gates
### V0 — provenance
Freeze the CFDX commit, build/compiler, numerical configuration, mesh generation inputs/version, reference dataset/version, and physical parameters.

### V1 — execution integrity
Demonstrate that CFDX actually executed the declared case. Retain executable/test identifier, exit status, configuration, solver log, and machine-readable result.

### V2 — convergence
Demonstrate numerical convergence using criteria appropriate to the case. Residual-only convergence is insufficient when physical conservation or QoI evidence is available.

### V3 — conservation
Independently check relevant mass, momentum, scalar, energy, CHT interface, and radiation balances. Acceptance thresholds must be declared by the case contract and must not be changed after observing a failure merely to obtain PASS.

### V4 — discretisation and temporal verification
For spatial or temporal accuracy claims, use a documented refinement sequence, fixed mathematical problem, declared refinement ratio, appropriate L1/L2/Linf norms, and observed order. Use GCI or another declared uncertainty estimate when appropriate.

### V5 — independent reference comparison
Freeze the reference, define the comparison metric before evaluating the result, compare matched QoIs, report absolute and relative differences, and distinguish reference uncertainty from CFDX numerical error.

### V6 — regression reproducibility
A promoted case must be reproducible from repository-controlled inputs. CI should retain machine-readable results, logs, plots/report, and exact test/build identifiers.

### V7 — auditability
The evidence must answer: what was run, with which numerical method, on which mesh, against which reference, using which acceptance criteria, and why was it promoted?

## 5. Quantitative accuracy and observed order
For errors E_i on systematically refined meshes with characteristic sizes h_i:

p_obs = ln(E_i / E_i+1) / ln(h_i / h_i+1)

For constant refinement ratio r:

p_obs = ln(E_coarse / E_fine) / ln(r)

Every order claim must identify the error norm, field/domain, boundary-cell treatment, refinement ratio, mesh family, and whether the solution is in the asymptotic range.

When GCI is used, document refinement ratio, observed order, safety factor, fine-grid error, extrapolated value where applicable, and whether the grids satisfy the method assumptions.

GCI is a solution-verification tool, not a validation metric.

Negative evidence is retained. CFDX must not weaken an order threshold, remove an inconvenient mesh, replace a difficult mesh with an easier one while keeping the same claim, or label polynomial exactness as second-order convergence.

## 6. Mesh and time-step studies
Refinement studies must distinguish mesh refinement, time-step refinement, iterative error, and model/reference discrepancy.

For transient studies, spatial and temporal errors must not be conflated. For mesh studies, retain topology, cell count, characteristic size, quality metrics, refinement ratio, and wall resolution where applicable.

A mesh labelled fine has no V&V meaning without a quantitative definition.

## 7. Reference-code comparisons
OpenFOAM, SU2, Fluent, CFX and STAR-CCM+ may be used as external technical references. They are not acceptance oracles merely because they produce a different numerical answer.

Cross-code campaigns must document code/version, mathematical formulation, gradient/reconstruction, convection/diffusion, pressure-velocity coupling, solver settings, mesh, boundary conditions, stopping criteria, QoIs, and differences.

Classify comparisons as: same mathematical method; similar method; different method with similar purpose; or not comparable.

Cross-code agreement is supporting evidence, not proof of correctness.

## 8. Qualification case record
Every promoted case should contain identity/scope, reference provenance, physical parameters, mesh family, boundary/initial conditions, numerical schemes, solver/preconditioner settings, convergence criteria, conservation criteria, QoI definitions, reference values, CFDX results, error metrics, observed order/GCI where applicable, runtime/hardware, regression evidence, verification result, validation result, uncertainty/limitations, and promotion decision.

Use docs/validation/CFDX_QUALIFICATION_CASE_TEMPLATE.md as the case-sheet template.

## 9. Numerical-method maturity traceability
Issue #461 follows:

requirement → implementation → executable test → verification evidence → validation evidence (if applicable) → qualification status

The chain may legitimately stop at verification when physical validation is not mathematically applicable.

Examples:
- Gradient operator: code and solution verification are primary; physical validation is normally indirect.
- Convection scheme: code/solution verification plus boundedness, conservation, and applicable order evidence are primary.
- Turbulence model: equation verification, solver verification, and physical validation are distinct.
- Linear preconditioner: algebraic verification, convergence/robustness evidence, and performance/scaling evidence are distinct.

## 10. CI policy
CI has two responsibilities.

### Diagnostic stage
Run tests and collect evidence even when an individual validation executable fails. Publish test output, numerical metrics, missing executables, solver failures, and generated reports.

### Gate stage
After evidence collection, enforce the declared acceptance gates. This keeps failures inspectable rather than hiding them behind an early CI exit.

Expensive total validation remains separate from ordinary PR CI and is requested explicitly by the repository validation workflow.

## 11. No false-positive qualification
The following are prohibited:
- tolerance inflation after a failed run;
- disabling a failing case to obtain green CI;
- fixed-iteration PASS;
- residual-only PASS when independent conservation/QoI evidence is required;
- hard-coded QoIs;
- using the reference oracle as the CFDX result;
- silently changing numerical schemes;
- silently changing mesh/reference configuration;
- claiming second-order accuracy from polynomial exactness alone;
- claiming GPU equivalence from GPU API availability;
- claiming reference-code parity from matching option names.

When evidence is incomplete, use PARTIAL, DIAGNOSTIC, READY, or BLOCKED.

## 12. Evidence retention and reproducibility
For each promoted campaign, retain CFDX revision, build/compiler/platform, case inputs, mesh provenance, numerical configuration, reference data/version, test command, machine-readable metrics, logs, figures, and final report.

Human-readable documents summarize evidence; machine-readable artifacts remain the regression source of truth.

## 13. Relationship between validation documents
- This document: governance, terminology, evidence, and gates.
- CFDX_QUALIFICATION_CAMPAIGN.md: qualification programme and mandatory case evidence.
- CFDX_QUALIFICATION_CASE_TEMPLATE.md: per-case calculation sheet.
- CFDX_QUALIFICATION_REGISTRY.json: machine-readable case status.
- NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json: requirement-to-evidence traceability for #461.
- NUMERICAL_METHOD_CAPABILITY_MATRIX.json: numerical-method capability metadata.
- NUMERICAL_METHOD_MATURITY_AUDIT.md: repository-level maturity interpretation.
- VALIDATION_REPORT.md: executable report-generation and CI evidence pipeline.
- Specialist documents: domain-specific evidence contracts.

No downstream document may promote a capability beyond the evidence supported by this governance layer.

## 14. Definition of qualification
A CFDX result is QUALIFIED only when the declared scope is explicit, implementation exists, applicable code verification is complete, applicable solution verification is complete, applicable physical validation is complete, conservation and convergence gates pass, quantitative error/order evidence is available where applicable, reference provenance is frozen, CI reproduces the evidence, and remaining limitations are documented.

Qualification is therefore a statement about a defined scope, not a claim that the entire CFDX solver is universally correct.

## 15. References and methodology
- ASME V&V 20, Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer.
- Celik et al., Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications / GCI methodology.
- Roache, Verification of Codes and Calculations.
- MMS-based code-verification methodology for manufactured solutions.

These references provide methodological guidance. CFDX acceptance thresholds remain case- and contract-specific and must be documented in repository evidence.

## 16. Documentation architecture

The V&V directory is a controlled evidence system, not a collection of independent notes.

Target structure:

docs/validation/
├── README.md                         # entry point and navigation
├── 00_governance/                    # terminology, gates, lifecycle, rules
├── 01_registries/                    # machine-readable source of truth
├── 02_method_verification/           # mathematical/operator verification
├── 03_solver_verification/           # complete solver/integration campaigns
├── 04_qualification/                 # promotion programme and templates
├── 05_cases/                         # frozen individual benchmark cases
├── 06_references/                    # provenance and reference datasets
└── 07_reporting/                     # report generation and publication

This is the target architecture. Migration must be incremental and must not break links or machine-readable consumers.

### 16.1 Single responsibility
Every document has exactly one primary role: governance, registry, method verification, solver verification, qualification, case, reference, or reporting. If a document accumulates another role, split it and link the two documents.

### 16.2 Source of truth
Each fact has one authoritative owner. V&V terminology/gates belong to governance; case status to the qualification registry; #461 requirement status to the requirements audit; capability catalogue to the capability matrix; reference provenance to reference records; case configuration/results to the case sheet plus retained machine-readable artifact; executable CI result to CI artifacts; human-readable summaries to reports. Other documents summarize and link; they do not create competing facts.

### 16.3 No duplicate matrices
Do not create a new matrix when an existing registry already owns the same dimension. A new matrix is allowed only when it answers a different question and explicitly states its scope. Capability, requirements, qualification, and campaign-result matrices are related but not interchangeable.

### 16.4 Document lifecycle
1. Proposal: define role and owner.
2. Scope: state what it does and does not claim.
3. Reference contract: identify mathematical/literature sources.
4. Implementation: add executable evidence where applicable.
5. Evidence: retain machine-readable outputs.
6. Review: audit against governance gates.
7. Promotion: update the relevant registry only after evidence exists.
8. Maintenance: update or explicitly supersede stale documents.
A document must never be created merely to make an issue appear complete.

### 16.5 Required evidence header
Every method, solver, qualification, and case document must state: Scope; Status; authoritative source/registry; evidence location; related issue/PR; last verified commit/date; known limitations. Numerical-result documents must additionally identify CFDX commit, build/compiler/platform, case ID, mesh/configuration ID, and test/CI job ID.

### 16.6 Naming rules
Use stable semantic names for new documents. Prefer N2_GRADIENT_VERIFICATION.md, N8_AMG_QUALIFICATION.md, and GHIA_RE100_CASE.md. Avoid new generic phase names such as PHASE9_ACCEPTANCE.md or NEW_FINAL_VALIDATION.md. Historical phase documents may remain while their authority is clarified.

### 16.7 Status synchronization
A result update must modify the authoritative registry and evidence document in the same PR. Never change Markdown status without updating the machine-readable registry; never mark PASS without retained executable evidence; never erase failed evidence. A regression must move the status back to the appropriate non-qualified state.

### 16.8 Promotion
Promotion is a state transition: PLANNED → IMPLEMENTED → VERIFIED → VALIDATED (when physical validation applies) → QUALIFIED. READY and DIAGNOSTIC are execution states, not maturity promotions.

### 16.9 Case implementation rule
A new qualification case starts with a stable case ID and registry entry, then frozen reference provenance, executable implementation, machine-readable result output, independent diagnostics, convergence/refinement studies where applicable, CI execution, case-sheet evidence, and only then registry promotion.

### 16.10 Failure handling
Failures are first-class V&V evidence. Retain failure class, failing metric, reproduction command, artifact, known root cause, and next action. Do not delete a case, disable a test, or rewrite acceptance criteria solely because it failed.

### 16.11 V&V PR review checklist
Before merge: correct document category; one authoritative owner per fact; explicit scope; implementation separated from verification/validation; frozen references; predeclared acceptance criteria; independent recomputation where appropriate; separated mesh/time/iterative error; retained failures; synchronized registry and narrative; retained CI artifacts; and sufficient evidence for any requested status transition.

### 16.12 Migration policy
The existing documentation contains valuable evidence but has grown organically. Do not perform a destructive mass rename. Migrate in small PRs: establish architecture and README; classify every document; identify overlaps; choose authoritative owners; add cross-links/deprecation notices; rename only after links/tooling are migrated; remove superseded duplicates only after replacement is proven. Specialist documents such as AMG, MGR, turbulence, conservation, and wall-distance remain separate until ownership is explicit.
