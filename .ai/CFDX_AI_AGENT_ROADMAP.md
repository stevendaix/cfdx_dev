# CFDX AI-Agent Roadmap

**Status:** roadmap / planning document  
**Parent epic:** [#557 — CFDX AI framework](https://github.com/stevendaix/cfdx_dev/issues/557)  
**Tracking issue:** [#968 — Make CFDX robust and effective for AI agents](https://github.com/stevendaix/cfdx_dev/issues/968)

## 1. Purpose

Make CFDX straightforward for AI agents to understand, operate, diagnose, validate, and improve while preserving numerical integrity, reproducibility, security, and the existing software architecture.

This roadmap is an incremental extension of the current `.ai/` framework. It is **not** a proposal to rebuild it from scratch. Existing skills, code intelligence, permission profiles, MCP implementations, execution artifacts, and case-persistence contracts must be audited and reused before new APIs are introduced.

The target operating model is:

- **Agents reason and orchestrate.**
- **CFDX's deterministic APIs and numerical core perform the work.**
- **Independent tests and validation criteria judge the evidence.**
- **Permissions and resource budgets constrain every operation.**
- **Every conclusion can be traced to source, inputs, commit, execution, and validation evidence.**

## 2. Non-negotiable engineering principles

1. **Repository-grounded:** current source, executable tests, current CI, and measured numerical results take precedence over model memory.
2. **Vendor-neutral:** canonical agent knowledge and contracts live in `.ai/`; provider-specific configuration references or derives from them.
3. **No duplicate domain semantics:** GUI, TUI, MCP adapters, and agent prompts reuse shared domain/application APIs rather than reimplementing solver rules.
4. **Numerical integrity:** never inflate tolerances, disable tests to get a green result, silently fall back to another method, or claim qualification from implementation or residual convergence alone.
5. **Evidence classification:** distinguish implementation, testing, verification, validation, and qualification. Preserve negative results and known limitations.
6. **Least privilege:** read-only by default; development, execution, write, and maintainer operations require distinct explicit permissions.
7. **Bounded execution:** expensive CFD runs require limits for time, iterations, memory/cores where enforceable, and output retention.
8. **Reproducibility:** meaningful experiments record commit SHA, effective configuration, input identities, environment, budget, command/tool, outputs, and stop reason.
9. **Small changes:** one coherent purpose per PR; tests and evidence must refer to the exact tested HEAD.
10. **Honest uncertainty:** missing evidence is not success. Use explicit outcomes such as `PASS`, `FAIL`, `NOT_RUN`, and `INCONCLUSIVE` where applicable.

## 3. Target architecture

Keep two logical MCP surfaces separate:

- **Development MCP:** repository navigation, code intelligence, documentation, tests, V&V metadata, CI and GitHub operations.
- **Runtime MCP:** canonical CFDX cases, mesh/case inspection, controlled configuration and execution, monitoring, checkpoints, fields, probes, metrics, and validation.

Both should consume stable, typed contracts and the existing permission model. Neither MCP is the source of truth for CFDX domain semantics.

The intended flow is:

`agent request → capability/permission check → typed domain operation → bounded execution → structured evidence → independent validation → evidence-based decision`

A multi-agent orchestration framework is optional. First make the CFDX APIs, contracts, and tools reliable; only add specialized agents when evaluation demonstrates a real benefit.

## 4. Phased implementation plan

The phases below are an ordered work plan, not claims that the work is absent. P0 must determine which items are already implemented and qualified. Do not create a new component when an existing one satisfies the requirement with adequate evidence.

### P0 — Audit current implementation and establish a gap register

**Work**
- Inventory `.ai/`, its rules, skills, agents, workflows, adapters, schemas, scripts, and current MCP implementations.
- Audit the code-intelligence index: deterministic build, incremental update, freshness, symbol/reference/dependency queries, parser provenance, and explicit test/V&V links.
- Audit Development MCP and Runtime MCP tool schemas, permissions, path confinement, timeout behavior, error handling, resource limits, and tests.
- Inspect current source, relevant tests, open PRs, and CI evidence on the current HEAD.
- Classify every capability as implemented, tested, verified, validated, qualified, incomplete, planned, or unknown, using only the evidence appropriate to that claim.
- Record gaps, dependencies, risks, and follow-up PRs in the existing epic/tracking issue.

**Deliverable:** evidence-backed gap register with source paths, tests, current status, risk, and acceptance criteria.

**Exit criteria**
- No new major API is proposed before the existing implementation has been checked.
- Planned behavior is clearly distinguished from shipped behavior.
- Every gap has a next action or a documented reason not to proceed.

### P1 — Contracts, schemas, and machine-readable evidence

**Work**
- Inventory existing contracts first; reuse the execution artifact and case persistence formats already established.
- Define only the missing versioned contracts needed for case summaries, solver capabilities, run manifests, validation reports, diagnostics, benchmarks, and patch/CI evidence.
- Represent status, units, missing values, stop reasons, and provenance explicitly.
- Record the exact commit, effective configuration, case/mesh identity, environment, resource budget, and artifact references for meaningful runs.
- Reject invalid or non-finite JSON values and inconsistent execution/convergence status.
- Add schema compatibility, malformed-input, missing-artifact, and failure-path tests.

**Exit criteria**
- A consumer can interpret reports without guessing from prose or default values.
- Execution success is not treated as scientific validity.
- Schema changes have documented compatibility/migration behavior.

### P2 — Agent knowledge and reusable skills

**Work**
- Audit and consolidate existing software skills (C++20, Python, CMake, testing, CI, review) and numerical skills (FVM, mesh topology/quality, gradients, boundary conditions, fluxes, pressure–velocity coupling, linear/nonlinear solvers, AMG/Schur, turbulence, multiphysics, conservation, V&V).
- Standardize each skill around applicability, prerequisites, procedure, actual available tools, tests, decision criteria, failure modes, references, and limitations.
- Keep skills small and composable; load task-relevant skills instead of one monolithic prompt.
- Use English source comments/docstrings and repository conventions where applicable.
- Add executable Python/Jupytext validation material where useful; avoid documents that contain only unexecutable claims.

**Exit criteria**
- Skills reference real paths and tools.
- Skills do not promote unsupported APIs or turn assumptions into facts.
- Representative tasks select the relevant skills without loading unrelated material.

### P3 — Qualify code intelligence

**Work**
- Test clean and incremental index builds, freshness checks, and deterministic output.
- Test symbol, reference, dependency, test, and validation queries against manually verified examples.
- Preserve parser backend, revision/content identity, coverage, failures, and diagnostics.
- Keep the conservative C/C++ scanner explicitly limited to coarse navigation; do not imply complete semantic coverage when parsing is partial.
- Create test-to-symbol or validation-to-symbol links only from explicit evidence; name similarity or directory proximity is not proof.
- Measure search precision and recall on a reviewed task set before considering vector search or external indexing services.

**Exit criteria**
- Stale/partial index results are visible to the agent.
- Results lead to verifiable repository locations.
- Absence of an index result is not interpreted as proof that a symbol or relationship does not exist.

### P4 — Harden Development MCP

**Work**
- Compare current tools with required workflows before adding or renaming tools.
- Provide stable typed capabilities for code/symbol/dependency search, documentation, test discovery/selection/execution, V&V evidence, CI inspection, and review preparation.
- Keep read-only, development, and maintainer permissions separate.
- Define tool timeouts, output limits, error types, and cancellation behavior.
- Test malformed input, stale indexes, path traversal, unavailable tools, subprocess failure, timeout, and large outputs.

**Exit criteria**
- Every advertised tool exists and has a tested schema and permission profile.
- Failures are surfaced rather than silently converted into empty or successful results.
- Read-only access cannot mutate the repository.

### P5 — Harden Runtime MCP and case lifecycle

**Work**
- Qualify case inspection, configuration, execution, monitoring, checkpoint/restart, result access, comparison, and validation operations.
- Preserve the canonical `.cfdx.h5` case-definition contract and the distinct state/checkpoint contract; do not invent parallel case schemas.
- Resolve and test the canonical case/session factory prerequisite before adding case creation.
- Keep write and execute gates independent; retain root confinement, symlink-escape protection, atomic persistence, clear no-overwrite behavior, and bounded timeouts.
- Make stop reasons and process outcomes explicit; document platform-specific termination and concurrency limits.
- Do not imply that synchronous Runtime MCP execution is a job scheduler.

**Exit criteria**
- Authorized operations are reproducible and bounded.
- Failure leaves inputs intact and does not leave misleading success artifacts.
- Runtime evidence remains descriptive; V&V claims come from the validation layer.

### P6 — Structured numerical diagnostics and performance

**Work**
- Expose structured histories for momentum/pressure residuals, velocity and pressure corrections, continuity/mass balance, linear iterations, preconditioner/setup/application timing, physical metrics, memory where available, and stop reasons.
- Distinguish linear convergence, nonlinear convergence, conservation, physical convergence, verification, validation, and qualification.
- Implement a three-stage diagnosis: deterministic symptom classification → evidence-backed hypotheses → controlled discriminating experiments.
- Benchmark solver and preconditioner variants with comparable case, mesh, matrix/RHS where applicable, initial state, stopping criteria, hardware, and build settings.
- Separate setup cost, solve/application cost, iteration counts, and end-to-end time.
- Optimize only subject to correctness, conservation, and relevant physical validation gates.

**Exit criteria**
- Known failures can be diagnosed without simply increasing iteration budgets or weakening tolerances.
- Benchmark reports explain whether comparisons are valid and preserve negative results.
- Recommendations state observed facts, measurements, interpretation, uncertainty, and the next discriminating test.

### P7 — Agent-assisted development workflow

**Work**
- Use isolated worktrees/branches and focused PRs.
- Standardize the workflow: requirement → repository search → impact analysis → plan → patch → focused tests → broader validation → independent review → CI → evidence report.
- Require an independent review/validation role for numerical changes.
- Bind test and CI evidence to the exact final HEAD.
- Make expensive tests explicit and keep required validation gates authoritative.
- Allow agents to prepare PRs, but do not let them bypass branch protection or declare a merge safe without required evidence.

**Exit criteria**
- A representative issue can be taken from investigation to a reviewable PR.
- The PR describes tests actually run, results, known limitations, and outstanding checks.
- Stale evidence is detected and rejected.

### P8 — GUI/TUI integration

**Work**
- Reuse the shared application/domain state and execution services.
- Expose agent proposals, permissions, effective parameters, budgets, campaign status, diagnostic evidence, and validation gaps.
- Keep domain rules out of GUI/TUI/MCP adapter-specific code.
- Add contract tests for equivalent operations across GUI, TUI, and MCP where applicable.

**Exit criteria**
- Interfaces do not disagree about case configuration or execution semantics.
- Important actions and their consequences are visible to the user.
- The user can inspect evidence and intervene before sensitive actions.

### P9 — Agent evaluation, hardening, and governance

**Work**
- Establish a fixed task suite covering code navigation, impact analysis, test selection, numerical diagnosis, conservation, V&V, performance comparison, CI analysis, and permission enforcement.
- Start with approximately 10 representative tasks, then expand to 30–50 as the suite matures.
- Measure diagnostic correctness, false positives/negatives, code-reference accuracy, regression rate, successful test selection, wall time, compute cost, tool calls, and unsupported claims.
- Establish a baseline and compare changes to models, prompts, skills, and orchestration against it.
- Add adversarial tests for prompt injection in repository content, unsafe paths, unauthorized operations, misleading logs, and incomplete evidence.
- Keep architecture, tool catalogue, schema compatibility, limitations, and evaluation results current.

**Exit criteria**
- The agent quality baseline is reproducible.
- Improvements are measured rather than inferred from anecdotal successes.
- Permission and numerical-integrity failures are treated as release blockers.

## 5. Recommended delivery order

Deliver small PRs in this order, adjusting only when P0 provides evidence of dependencies or existing completion:

1. Audit and gap register.
2. Missing contracts and evidence tests.
3. Code-index reliability and evaluation.
4. Skills and task fixtures.
5. Development MCP gaps.
6. Runtime MCP gaps.
7. Numerical diagnostics and benchmark comparability.
8. Controlled agent development workflow.
9. GUI/TUI integration.
10. Agent evaluation and hardening.

Parallel work is acceptable only where ownership and dependencies are clear. Do not open speculative PRs for capabilities that the audit shows are already implemented.

## 6. PR acceptance template

Every implementation PR should state:

- **Scope:** issue and acceptance criteria addressed.
- **Evidence inspected:** paths, APIs, tests, and prior behavior.
- **Changes:** smallest coherent implementation.
- **Tests:** commands and exact results, distinguishing local from CI.
- **Numerical evidence:** formulation/operator implications, conservation/stability considerations, reference and tolerances where relevant.
- **Provenance:** commit SHA, case/configuration, environment and campaign IDs where relevant.
- **Limitations:** unrun tests, unsupported configurations, known risks and remaining work.
- **Decision:** PASS, FAIL, NOT_RUN, or INCONCLUSIVE for each applicable criterion.

A green CI is evidence that configured checks passed, not a general qualification claim.

## 7. Tracking and status

The detailed implementation status belongs in [issue #968](https://github.com/stevendaix/cfdx_dev/issues/968), under the parent [epic #557](https://github.com/stevendaix/cfdx_dev/issues/557).

For each phase, record:
- status: planned / in progress / implemented / tested / verified / validated / qualified;
- relevant PR and exact merged commit;
- applicable tests and CI evidence;
- unresolved risks and the next action.

Do not mark a phase complete merely because documentation or code exists. Update this roadmap when repository evidence changes the plan; keep it consistent with the canonical `.ai/` rules and current implementation.
