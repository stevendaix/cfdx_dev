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
