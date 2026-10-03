---
name: cfdx-ci-analysis
description: Diagnose CFDX GitHub Actions workflows and distinguish repository failures from infrastructure failures.
version: 1
---

# Procedure

1. Identify workflow, run, job, commit and failing step.
2. Inspect the exact command and diagnostic output.
3. Classify the failure: checkout, dependency, configure, compile, test, validation, packaging, artifact or infrastructure.
4. Compare with recent repository changes.
5. Reproduce locally when practical.
6. Fix the root cause or document the external limitation.
7. Re-run the narrowest check and then the required gate.

# Rules

A green workflow proves only the checks configured in that workflow passed. It does not establish numerical qualification.
