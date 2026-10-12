# CFDX AI Global Operating Rules

These rules apply to every CFDX AI agent, workflow and MCP unless a more restrictive rule is explicitly defined.

## 1. Repository-grounded operation

- Inspect the current repository state before making or asserting changes.
- Prefer current executable code, tests, CI and repository documentation over model memory.
- Search broadly enough to identify callers, interfaces, tests and documentation affected by a change.
- Do not invent files, APIs, symbols, tests, CI results, benchmark results or repository history.
- If a tool fails, a result is unavailable, or evidence is incomplete, say so explicitly.

## 2. Scope and requirements

- Identify the issue, acceptance criteria and intended maturity level before implementation.
- Keep requested scope separate from opportunistic cleanup.
- Do not silently change requirements, acceptance criteria or test meaning.
- Preserve backward compatibility unless a breaking change is explicitly intended and documented.
- Identify assumptions and unresolved requirements before relying on them.

## 3. Engineering integrity

- Prefer the smallest coherent change that satisfies the requirement.
- Preserve existing architecture and ownership boundaries unless the task explicitly changes them.
- Do not duplicate domain semantics in GUI, TUI, adapters or agents when a shared core API/state model exists.
- Do not hide errors with broad exception handling, silent fallbacks, automatic degradation or implicit recovery.
- Fail clearly when a required capability is unavailable.

## 4. Numerical integrity

- Never alter tolerances, iteration limits, convergence criteria or reference values merely to make a test pass.
- Never disable, skip or weaken a test merely to obtain green CI.
- Do not claim numerical correctness, verification, validation or qualification from implementation existence or residual convergence alone.
- Preserve negative evidence, failed cases and known limitations.
- For numerical changes, state the continuous formulation, discrete operator and relevant conservation/stability/consistency implications.

## 5. Evidence and claims

Every substantive conclusion should be traceable to evidence.

Classify claims as applicable:
- implementation;
- test;
- verification;
- validation;
- qualification.

A green CI check proves only the configured checks passed. It does not prove broader numerical qualification.

When reporting status, distinguish:
- observed fact;
- measured result;
- interpretation;
- unresolved risk.

Do not turn an absence of evidence into evidence of success.

## 6. Change workflow

Use this default sequence:

requirement -> repository search -> impact analysis -> plan -> implementation -> focused tests -> broader tests/validation -> documentation -> diff/review -> CI -> evidence-based report.

After each material failure, investigate the cause before changing acceptance criteria.

For PR work:
- keep commits/PRs focused;
- rebase or recreate stale branches when necessary;
- inspect the final diff against the intended scope;
- do not merge until the required CI and review conditions are actually satisfied.

## 7. Safety and permissions

- Follow the active permission profile and least-privilege principle.
- Do not perform destructive operations without the required authorization.
- Bound expensive CFD workloads unless explicitly approved.
- Never expose secrets, credentials or private data in logs, generated files or reports.
- Treat generated scripts and external inputs as untrusted until inspected.

## 8. Reproducibility

- Record commands, configuration, commit SHA and relevant environment information for meaningful experiments.
- Prefer deterministic tests and controlled seeds where applicable.
- Distinguish local results from CI results.
- Make expensive or long-running validation explicitly identifiable.

## 9. Documentation and synchronization

- Documentation must describe the current implementation, not intended future behavior.
- Update affected user/developer/theory/V&V documentation when behavior or public interfaces change.
- Keep agent rules, skills and workflows consistent; do not create contradictory local instructions.

## 10. Completion criteria

An agent must not declare a task complete merely because code was changed.

Completion requires the applicable acceptance criteria, tests/evidence, documentation and CI/review conditions to be checked. Any remaining gap must be stated explicitly.
