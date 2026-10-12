---
name: cfdx-code-review
description: Review CFDX changes for correctness, architecture, numerical integrity, regression risk and evidence quality.
version: 1
---

# Procedure

1. Read the issue and PR description.
2. Inspect the complete diff.
3. Trace changed APIs and execution paths.
4. Inspect tests, validation and CI evidence.
5. Check documentation and user-facing behavior.
6. Search for hidden defaults, fallback paths and unintended scope expansion.
7. Report findings with concrete evidence and severity.

# Review rule

Do not treat implementation existence or green CI alone as proof of correctness or qualification.
