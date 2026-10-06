# Generic CFDX validation tools

These utilities are intentionally solver-agnostic **data-quality and comparison
primitives**. They are not a qualification engine.

## Data flow

1. CFDX (or another solver) produces standardized CSV/JSON evidence.
2. `aggregate_validation_runs.py` combines runs deterministically and records
   source-file provenance.
3. `compare_validation_results.py` compares explicitly requested QoIs.
4. A campaign-specific Jupytext chapter performs convergence analysis, plots,
   scientific interpretation and qualification decisions.

No NACA-specific reference value or acceptance tolerance belongs in these tools.

## Aggregation

```bash
python scripts/aggregate_validation_runs.py \
  --input-dir results/n9-s8 \
  --pattern '*.csv' \
  --key case \
  --output results/n9-s8/all.csv
```

The aggregator rejects missing/invalid headers, malformed rows and duplicate
keys. It creates the output directory and excludes the output file itself from
the input glob. The `source_file` column preserves provenance.

## Comparison

```bash
python scripts/compare_validation_results.py \
  --computed results/n9-s8/cfdx.csv \
  --reference references/n9-s8.csv \
  --qoi Cd Cl \
  --output results/n9-s8/comparison.json
```

The comparator rejects missing cases, missing QoIs, duplicate cases,
non-numeric values and non-finite values. A zero reference has no relative-error
denominator and is reported with `relative_error: null`.

These utilities deliberately do **not** assign PASS/FAIL. Campaign contracts
own tolerances and scientific decisions.
