# N14 — Reference-code comparison matrix

Issue #461 N14 is a traceability and methodology-comparison package. It documents how CFDX numerical methods compare with publicly documented methods in OpenFOAM, SU2, Ansys Fluent, Ansys CFX and Siemens Simcenter STAR-CCM+.

## Rules

1. CFDX is the starting point: the matrix describes what CFDX actually implements.
2. A shared method name is never evidence of mathematical equivalence.
3. External codes are references for formulation/capability documentation, not acceptance oracles.
4. Implementation, verification, validation and qualification remain independent statuses.
5. Every comparison is classified as one of:
   - `same_mathematical_method`
   - `similar_method`
   - `different_method_same_purpose`
   - `unavailable_in_cfdx`
6. Proprietary or undocumented implementation details are recorded as `not_publicly_documented`, not inferred.
7. Benchmark compatibility is reported only when the CFDX case and reference definition are demonstrably comparable.

## Artifacts

- `N14_REFERENCE_CODE_MATRIX.json`: machine-readable schema and comparison records.
- `N14_REFERENCE_CODE_MATRIX.md`: human-readable matrix and audit rules.
- `openfoam.md`, `su2.md`, `fluent.md`, `cfx.md`, `starccm.md`: source dossiers and scope limits.
- `N14_BENCHMARK_CROSS_REFERENCE.md`: links between N14 comparisons and existing N12/N13/N4/N8/N9 evidence.

## Evidence vocabulary

- `official_documentation`: vendor/project documentation describing the numerical method.
- `official_source`: inspectable source code from an open-source reference.
- `publication`: peer-reviewed or institutional technical publication.
- `cfdx_executable`: executable CFDX verification/validation evidence.
- `not_available`: no sufficiently authoritative public formulation was located.

N14 does not create PASS/FAIL gates against a reference code. Quantitative acceptance remains governed by the relevant CFDX V&V package.
