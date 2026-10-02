# CFDX Documentation Tooling

This directory is intentionally small. It contains documentation-only orchestration: execution of designated Theory sources, collection of machine-readable results, figure generation, reference checks and publication helpers.

It must not contain CFDX numerical-domain logic. Numerical algorithms remain in CFDX itself or in explicitly labelled independent reference experiments under `docs/theory`.

## Source policy

- Python `py:percent` is the executable Theory source.
- Markdown is navigation, contracts and prose where execution is not useful.
- Generated HTML, notebooks, figures and reports are build products.
- No custom web application is introduced.

## Reproducibility

Executable studies must record the CFDX revision when they call CFDX, numerical parameters, mesh/resolution, solver tolerances and relevant platform/dependency information.
