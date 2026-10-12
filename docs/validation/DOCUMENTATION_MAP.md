# Validation Evidence Map

This file is retained as a migration map. It is **not** a numerical status matrix.

## Canonical ownership

- V&V methodology: `docs/vv/*/chapter.py`.
- Executed validation campaigns: `docs/validation/**/chapter.py`.
- Benchmark campaigns: `docs/validation/14-benchmarks/*/chapter.py`.
- Qualification status: `CFDX_QUALIFICATION_REGISTRY.json`, backed by retained executable evidence.
- Capability/maturity catalogues: their JSON registries.
- Report mechanics: `scripts/validation_report.py`.

## Migration decisions

| Existing material | Decision |
|---|---|
| Numbered directory README files | Keep as navigation/scope. |
| Benchmark `chapter.py` files | Keep as canonical executable campaign sources. |
| Qualification registry JSON | Keep as machine-readable status source. |
| V&V governance Markdown | Keep only as governance/provenance while methodology remains in `docs/vv/`. |
| Specialist qualification/campaign Markdown | Audit and migrate facts incrementally; do not delete evidence blindly. |
| Hand-maintained quantitative matrices | Treat as migration/reporting views, not new sources of truth. |
| Reference-code descriptions | Keep as provenance/reference material. |
| Report-generation documentation | Keep as user/developer documentation for the reporting pipeline, not campaign evidence. |

## Conversion rule

When a specialist Markdown file contains executable mathematics, quantitative acceptance criteria, measured results, or reproducible data reduction, its **new canonical form** should be a Jupytext percent-format `chapter.py`. The Markdown may then become a short pointer or be removed after all inbound links and ownership are audited.

Do not fabricate a campaign notebook from prose alone. Migration creates executable source only when the underlying evidence or computation can actually be reproduced.

## Current priority

1. Benchmark campaigns — already Python-first; standardize their common evidence contract.
2. Lower-level verification campaigns 00–11 — add executable chapters when corresponding tests/evidence exist.
3. Qualification — make the executable campaign/report and JSON registry authoritative; keep prose as a view.
4. Legacy specialist Markdown — classify, migrate or retire one document at a time.
5. CI — publish the executable sources and machine-readable evidence, while excluding standalone support scripts from the Sphinx source set.
