# CFDX Verification, Validation & Qualification

This directory is the controlled evidence system for CFDX numerical verification, solution verification, validation and qualification.

## Start here
1. Governance: CFDX_VV_GOVERNANCE.md
2. Qualification programme: CFDX_QUALIFICATION_CAMPAIGN.md
3. Case template: CFDX_QUALIFICATION_CASE_TEMPLATE.md
4. Case registry: CFDX_QUALIFICATION_REGISTRY.json
5. #461 requirements audit: NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
6. Numerical maturity audit: NUMERICAL_METHOD_MATURITY_AUDIT.md
7. Executable reporting: VALIDATION_REPORT.md

## Document roles
| Role | Purpose | Source of truth |
|---|---|---|
| Governance | V&V definitions, gates, lifecycle and documentation rules | V&V governance |
| Registry | Machine-readable status and traceability | JSON registry |
| Method verification | Mathematical/operator evidence | Method evidence |
| Solver verification | Complete solver/integration evidence | Solver campaigns |
| Qualification | Promotion programme and case requirements | Qualification campaign |
| Case | One frozen benchmark configuration and its evidence | Case sheet + artifacts |
| Reference | Provenance and reference-data definition | Reference record |
| Reporting | Generation/publication of campaign evidence | Reporting pipeline |

## Current directory migration
The directory currently contains legacy/organic documents as well as authoritative registries and specialist campaigns. Do not create another document until its role and authoritative owner are identified.

The target architecture and migration rules are defined in CFDX_VV_GOVERNANCE.md, section 16.

## Mandatory rules for new V&V work
- Define the scope before implementation.
- Assign a stable case or requirement ID.
- Identify the authoritative registry.
- Freeze reference provenance before comparing results.
- Build an executable test before declaring solver evidence.
- Retain machine-readable evidence.
- Separate code verification, solution verification and validation.
- Define acceptance criteria before observing the result.
- Keep failures visible.
- Update registry and evidence documentation in the same PR.
- Never infer qualification from implementation presence.

## Important distinction
An analytical/reference oracle is not a CFDX solver result. A component test is not a complete scheme verification. A solver convergence result is not physical validation. A cross-code match is not proof of correctness.

Qualification is always scoped to the declared mathematical model, numerical method, physical configuration and evidence gates.