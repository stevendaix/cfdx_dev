# CFDX V&V Documentation Map

This map is the migration inventory for docs/validation. It is deliberately separate from the qualification registry: it governs documentation ownership, not numerical status.

| Current document | Role | Target area | Action |
|---|---|---|---|
| CFDX_VV_GOVERNANCE.md | Governance | 00_governance | authoritative |
| CFDX_QUALIFICATION_CAMPAIGN.md | Qualification | 04_qualification | authoritative |
| CFDX_QUALIFICATION_CASE_TEMPLATE.md | Case template | 04_qualification | authoritative |
| CFDX_QUALIFICATION_REGISTRY.json | Registry | 01_registries | authoritative |
| NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json | Registry | 01_registries | authoritative for #461 |
| NUMERICAL_METHOD_CAPABILITY_MATRIX.json | Registry | 01_registries | authoritative capability catalogue |
| NUMERICAL_METHOD_MATURITY_AUDIT.md/.json | Audit | 01_registries | consolidate/clarify paired human+machine roles |
| NUMERICAL_MODEL_VERIFICATION_MATRIX.md | Method verification | 02_method_verification | authoritative method matrix |
| CONSERVATION_BOUNDEDNESS.md | Method verification | 02_method_verification | retain as specialist contract |
| AMG_PRECONDITIONER_QUALIFICATION.md | Solver verification | 03_solver_verification | retain as N8 specialist campaign |
| MGR_PRECONDITIONER_QUALIFICATION.md | Solver verification | 03_solver_verification | retain as N8 specialist campaign |
| SOLVER_LEVEL_BENCHMARKS.md | Solver verification | 03_solver_verification | authoritative solver benchmark overview |
| PHASE3_6_SKEW_CAMPAIGN.md | Solver verification | 03_solver_verification | rename after link audit |
| PHASE9_ACCEPTANCE.md | Solver verification | 03_solver_verification | rename after link audit |
| M1_M4_VERIFICATION_PLAN.md | Method verification | 02_method_verification | split if it mixes plan and evidence |
| M1_M4_ANALYTICAL_BENCHMARKS.md | Method verification | 02_method_verification | retain as analytical benchmark catalogue |
| EXTENDED_BENCHMARK_MATRIX.md | Benchmark catalogue | 03_solver_verification | reconcile overlap with solver benchmarks |
| PERFORMANCE_AND_PHYSICS_WAVE.md | Campaign plan | 03_solver_verification | classify explicitly; avoid duplicate status |
| VALIDATION_REPORT.md | Reporting | 07_reporting | authoritative pipeline documentation |
| LEVEL_B_REFERENCES.md | Reference | 06_references | authoritative provenance index |
| FLUENT_VMFL_MATRIX.md | Reference/compatibility | 06_references | retain; not qualification authority |
| SU2_NACA0012_READER.md | Reference | 06_references | retain as source-data reader documentation |
| *_QUALIFICATION_CASE.md | Case | 05_cases | migrate by case family after registry audit |
| VMFL036_* | Case | 05_cases | migrate as VMFL036 case evidence |
| TURBULENCE_*.md/json | Specialist method/model | 02_method_verification + 05_cases | keep model equations separate from physical cases |
| THERMOPHYSICAL_MODELS.md | Method/model verification | 02_method_verification | classify against numerical-model matrix |
| THERMAL_RADIATION_REGRESSION_MATRIX.md | Solver verification | 03_solver_verification | retain as regression campaign |
| WALL_DISTANCE_TURBULENCE_INTEGRATION.md | Specialist dependency | 02_method_verification | retain; link from turbulence ownership |

## Important overlaps to resolve

### 1. Capability vs requirements vs qualification
NUMERICAL_METHOD_CAPABILITY_MATRIX, NUMERICAL_METHOD_REQUIREMENTS_AUDIT, NUMERICAL_METHOD_MATURITY_AUDIT, and CFDX_QUALIFICATION_REGISTRY answer different questions. They must not be collapsed into one giant matrix.

### 2. Benchmark catalogues
EXTENDED_BENCHMARK_MATRIX, SOLVER_LEVEL_BENCHMARKS, M1_M4_ANALYTICAL_BENCHMARKS and CFDX_QUALIFICATION_CAMPAIGN overlap in case inventories. The qualification registry should own case status; benchmark documents should describe scope and methodology only.

### 3. Phase documents
PHASE3_6, PHASE9 and PHASE13 documents encode historical implementation waves. They are useful engineering history but should not become parallel authorities for current qualification status.

### 4. Turbulence
The turbulence JSON is the machine-readable model registry. TURBULENCE_QUALIFICATION_MATRIX.md explains model-level qualification semantics. TURBULENCE_CLOSURE_EQUATIONS.md owns equation/reference contracts. Physical case documents own case evidence. These roles should remain separate.

## Migration rule
Do not move all files in one PR. First establish ownership and links. Then migrate one family at a time, with CI and repository-wide link checks. Delete or archive a duplicate only after its facts have been absorbed by the authoritative document.