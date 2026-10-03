---
name: cfdx-vv-workflow
description: Plan and execute reproducible CFDX verification and validation workflows.
version: 1
---

# Purpose

Provide a common procedure for CFDX V&V work from requirement to evidence-backed conclusion.

# Scope

Analytical solutions, MMS, benchmark/reference cases, mesh/time refinement, convergence, conservation, boundedness, uncertainty and qualification evidence.

# Procedure

1. Identify the requirement and maturity/qualification level.
2. Define the quantity of interest and acceptance criterion before running the case.
3. Identify the independent reference or analytical basis.
4. Define mesh/time/refinement and convergence protocols.
5. Locate the existing CFDX case, implementation and tests.
6. Run the smallest diagnostic case first.
7. Record configuration, commit, environment and commands.
8. Measure the predefined quantities.
9. Separate solver convergence from discretisation accuracy and model validation.
10. Preserve failures, anomalies and non-monotone trends.
11. Compare results against the predefined criterion.
12. Report implementation/test/verification/validation/qualification status separately.

# Failure handling

Do not change tolerances, reference values or acceptance criteria after seeing an unfavorable result unless the change is independently justified and explicitly documented. Do not suppress a failing case.

# Expected output

Case definition, protocol, commands, measurements, plots/tables where applicable, acceptance result, limitations and reproducibility metadata.
