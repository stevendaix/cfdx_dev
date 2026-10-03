# CFDX AI Workflows

Workflows describe repeatable multi-step procedures. They orchestrate skills and tools; they do not replace them.

Core workflows:

- implement-feature
- fix-bug
- implement-numerical-method
- audit-pr
- audit-issue
- investigate-ci-failure
- prepare-validation
- update-documentation
- develop-gui-feature
- develop-tui-feature
- run-cfdx-case

Standard sequence: requirement -> repository/documentation search -> dependency/impact analysis -> plan -> change -> tests -> validation -> documentation -> diff review -> CI -> evidence-based report.

Numerical changes additionally require the mathematical/discrete formulation and independent verification/validation evidence.
## Mandatory cross-cutting rules

All workflows inherit `.ai/AGENTS.md`. A workflow may add stricter requirements but must not weaken them. In particular, workflows must preserve negative evidence, distinguish required from optional checks, avoid silent fallback/error swallowing, and report tool or evidence failures explicitly.

For PR/CI workflows, the agent must inspect the final repository state and actual CI results before declaring completion or merge readiness. A planned command is not evidence that the command succeeded.
