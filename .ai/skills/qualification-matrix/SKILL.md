---
name: cfdx-qualification-matrix
description: Map CFDX tests and V&V evidence to explicit maturity and qualification criteria.
version: 1
---

# Purpose

Maintain an auditable mapping between requirements, cases, tests, evidence and qualification status.

# Procedure

For each criterion record:

- requirement;
- case/test;
- required or optional status;
- implementation exercised;
- evidence type;
- measured result;
- acceptance criterion;
- current status;
- evidence location;
- known limitation.

Use precise states such as implemented, tested, verified, validated, qualified, partial, blocked or unknown.

# Rules

Do not infer qualification from the existence of a test or from green CI alone. Preserve failed and deferred cases. If a criterion is not actually exercised, mark it as not demonstrated rather than passing it by inference.
