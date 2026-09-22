# CFDX validation report generation

The validation campaign produces a reproducible PDF from the executable CFDX results.

## Design

The pipeline deliberately separates four layers:

1. **Reference definition** — the Fluent VMFL matrix and published analytical/reference values.
2. **CFDX execution** — CTest executables run the numerical cases.
3. **Data reduction** — `scripts/validation_report.py` parses the executable output, computes report metrics and creates figures.
4. **Document rendering** — Python writes the LaTeX document and `latexmk` produces the final PDF.

This is preferable to generating a PDF directly from Python because the numerical data remain machine-readable while LaTeX controls typography, tables, references and pagination. `latexmk` also automatically resolves LaTeX dependencies and repeated compilation passes.

## Local use

After building CFDX:

```bash
python3 scripts/validation_report.py \
  --build-dir build \
  --output-dir build/validation-report
```

The output directory contains:

- `cfdx_validation_report.pdf` — human-readable report;
- `cfdx_validation_report.tex` — exact LaTeX source used to create the PDF;
- `results.json` — machine-readable extracted results;
- executable logs;
- generated figures.

## Interpretation

The report intentionally distinguishes:

- **reference oracle**: analytical/reference calculation only;
- **solver result**: CFDX actually solved the physical configuration;
- **PASS/FAIL**: quantitative solver-level assessment;
- **READY/PARTIAL/BLOCKED**: capability state from the VMFL matrix.

A reference formula is therefore never promoted to a CFDX solver PASS.

The current report contains the full VMFL001--VMFL078 matrix and executable results for the validation tests available in the build. As additional Fluent-aligned cases are implemented, the same reporting interface should consume their numerical result records rather than adding ad-hoc report code.

## GitHub on-demand campaign

The complete validation campaign can be launched manually from GitHub Actions with the
**CFDX VMFL validation campaign** workflow. It builds the current validation suite, runs
`scripts/validation_report.py`, and publishes the complete evidence as workflow artifacts.

The artifact contains:

- `cfdx_validation_report.pdf` — human-readable validation report;
- `results.json` — machine-readable solver/reference status and validation gate;
- executable stdout/stderr logs;
- generated figures and LaTeX source;
- CTest evidence from the build.

A workflow run is intentionally diagnostic as well as gating: missing validation
executables and non-zero solver exits remain visible in the report, while the workflow
uploads the evidence before enforcing the validation gate.

This workflow is **manual (`workflow_dispatch`)** so the full campaign can be run on
demand without consuming CI resources on every commit. It should be used after a
validation implementation, numerical fix, or solver change to inspect the complete
campaign and identify regressions.
