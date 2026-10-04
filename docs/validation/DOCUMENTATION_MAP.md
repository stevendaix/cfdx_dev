# CFDX Validation Documentation Map

This map defines **ownership and migration**, not numerical status.

## Canonical architecture

| Area | Primary source | Human-facing layer | Machine/evidence layer |
|---|---|---|---|
| V&V methodology | `docs/vv/*/chapter.py` | `docs/vv/README.md` | executable checks / CI |
| Validation campaigns | `docs/validation/**/chapter.py` | directory `README.md` | JSON + retained artifacts |
| Benchmark cases | `14-benchmarks/*/chapter.py` | `14-benchmarks/README.md` | campaign artifacts |
| Qualification | executable qualification campaign | `15-qualification/README.md` | `CFDX_QUALIFICATION_REGISTRY.json` |
| Capability catalogue | — | concise generated/summary view | `NUMERICAL_METHOD_CAPABILITY_MATRIX.json` |
| #461 maturity audit | executable evidence + registry | audit summary | `NUMERICAL_METHOD_MATURITY_AUDIT.json` |
| Reporting | `scripts/validation_report.py` | `VALIDATION_REPORT.md` | `results.json` / CI artifacts |

## Markdown policy

Markdown remains appropriate for navigation and directory READMEs, governance/policy that is not an executable scientific campaign, provenance/reference descriptions, and migration/history documents that have not yet been absorbed.

Markdown is **not** the preferred source for new quantitative campaign results, observed-order claims, benchmark measurements, or qualification status.

## Current migration

Validation directories 00–11 currently contain short READMEs but no campaign `chapter.py`. That is acceptable until executable evidence exists. When a campaign is implemented, add a Jupytext `chapter.py` and keep the README as navigation.

The benchmark family already follows the desired model: each case has an executable `chapter.py`, while the family README provides scope and navigation.

## Legacy flat documents

The flat specialist Markdown files under `docs/validation/` are **not deleted by this cleanup**. They are migration candidates. Before removal, classify their facts as V&V methodology (link to `docs/vv/`), executable campaign (migrate to `chapter.py`), machine-readable status (move into JSON), provenance/reference material (retain), or historical plan (retain explicitly until no longer useful).

No duplicate is deleted solely because a newer file exists.

## Status ownership

PASS/FAIL/READY/QUALIFIED state must come from executable evidence and the machine-readable registry. Human-facing Markdown may explain the state and scope, but must not become a manually maintained competing status matrix.
