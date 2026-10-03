# CFDX Development and Validation Rules

This document defines the repository rules that apply to implementation, testing, CI, numerical maturity, and qualification.

The objective is to keep CFDX development fast and reproducible **without weakening numerical assurance**.

## 1. Core principles

1. **Do not make the CI green by weakening the numerical contract.**
   - Do not relax tolerances.
   - Do not increase acceptance thresholds.
   - Do not disable failing tests to obtain a green run.
   - Do not remove qualification cases because they are expensive.
   - Do not shorten physically meaningful campaigns.

2. **Optimize infrastructure, not numerical evidence.**
   Safe optimization targets include build caching and reuse, CMake/CTest parallelism, artifact reuse, dependency caching, workflow scheduling, avoiding redundant compilation, and diagnostic collection.

3. **Implementation and qualification are different states.**
   A feature can be implemented and tested without being numerically qualified. Issue #461 must only mark a qualification item complete when the required evidence and acceptance criteria are actually demonstrated.

4. **A green test is evidence, not automatically qualification.**
   Qualification requires the appropriate campaign, reference solution or manufactured solution where applicable, convergence/order/conservation metrics, and explicit acceptance criteria.

## 2. Pull-request rules

Each PR should have a narrow, auditable scope.

A PR description should state:
- the issue or maturity level addressed;
- what is implemented;
- what is tested;
- what is actually qualified;
- what remains outside the scope;
- any known limitations.

Small follow-up PRs are preferred to large mixed-purpose changes.

When a branch becomes stale, rebase it onto the current target branch before the final validation when practical.

## 3. Merge-gate CI

The normal merge gate is designed to be fast while preserving its mandatory coverage.

### CFDX CI

The merge gate runs Release and DebugSanitizers.

The matrix uses `fail-fast: false` so an early failure does not hide independent failures.

The fast CTest selection excludes tests classified as `long` (and, for sanitizer execution, `benchmark`). This is an execution-time classification, not a statement that those tests are unimportant.

The merge gate must retain:
- the complete mandatory fast test selection;
- Python regression tests where required;
- diagnostics before the final decision;
- an explicit final pass/fail gate.

`--no-tests=error` is used where appropriate so an unexpectedly empty or missing test selection cannot silently pass.

## 4. Numerical maturity workflow

The numerical maturity workflow covers N1–N16.

Rules:
- Build the complete Release numerical suite once per workflow run when possible.
- Reuse the resulting build artifact across maturity levels rather than recompiling the same code for every level.
- Keep `fail-fast: false` so all maturity levels remain independently diagnosable.
- Preserve the test regexes and qualification scope of each N level.
- Keep diagnostics before the per-level final gate.
- N16 is an audit/coverage level and must not introduce an unnecessary C++ rebuild.
- Reusing artifacts must preserve executable permissions and the exact build output.

A CI optimization is acceptable only if it preserves the same numerical test content and acceptance criteria.

## 5. Long and full validation

Long or physically expensive campaigns are intentionally separated from the fast merge gate.

They may be activated explicitly when required for qualification or release evidence.

Being excluded from the fast merge gate does **not** make a long campaign optional when it is part of a qualification requirement.

For transient campaigns in particular:
- do not reduce the physical time horizon solely to save CI time;
- do not reduce mesh levels needed to establish an order of accuracy;
- do not replace a production campaign with a cheaper surrogate unless the qualification specification explicitly allows it.

## 6. Numerical integrity rules

For numerical-method changes, preserve independent evidence whenever possible.

Examples include conservation and flux-balance checks; boundedness and positivity; manufactured-solution convergence; spatial and temporal order; Richardson extrapolation; restart equivalence; nonlinear convergence and stagnation/divergence evidence; reference-solution comparisons; and solver/preconditioner energy or contraction metrics.

Negative tests are evidence too. A documented rejection, stagnation, divergence, or invalid-input path must not be converted into a pass merely to simplify the gate.

## 7. Diagnostics before gates

Tests should expose useful diagnostics before the final pass/fail decision.

A failure should identify, where applicable:
- the test/case;
- mesh or resolution level;
- measured metric;
- expected acceptance criterion;
- convergence/order calculation;
- relevant iteration/time information;
- whether the failure is an implementation defect, numerical regression, infrastructure failure, or missing qualification evidence.

Do not hide a failure behind a later aggregate gate.

## 8. Qualification and issue #461

Issue #461 is the numerical maturity/qualification roadmap.

For each item:
- `[x]` means the required implementation **and qualification evidence** is complete;
- an implemented feature without sufficient campaign evidence remains unchecked or explicitly marked partial;
- a green CI run alone does not justify checking a qualification item;
- related PRs should be referenced so the evidence remains traceable.

When a PR closes only part of a maturity level, update the issue to reflect the remaining boundary rather than declaring the whole level complete.

## 9. Safe CI optimization checklist

Before merging a CI optimization, verify:

- [ ] No numerical tolerance changed.
- [ ] No acceptance threshold changed.
- [ ] No required test was disabled.
- [ ] No qualification case was removed.
- [ ] No physically meaningful campaign was shortened.
- [ ] No mesh/resolution level required for an order study was removed.
- [ ] The same mandatory test selection remains.
- [ ] Artifact reuse preserves executable permissions and build contents.
- [ ] Diagnostics remain available.
- [ ] Final gates still fail on genuine numerical regressions.

## 10. Rule of thumb

**Make the computation and CI infrastructure faster; never make the numerical proof weaker just to make CI faster.**

When speed and qualification evidence conflict, optimize the build, scheduling, caching, artifact transport, or test orchestration first.
