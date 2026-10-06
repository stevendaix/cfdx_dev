# Generic CFDX validation tools

These tools are deliberately solver-agnostic. They separate data production from
comparison so the same workflow can be reused for numerical maturity cases,
N9 external-flow cases, and later physical qualification.

## 1. Aggregate runs

```
python scripts/aggregate_validation_runs.py --input-dir results/n9-s8 --pattern '*.csv' --output results/n9-s8/all.csv
```

## 2. Compare CFDX with an external reference

```
python scripts/compare_validation_results.py --computed results/n9-s8/cfdx.csv --reference references/n9-s8.csv --qoi Cd Cl
```

The output reports absolute and relative discrepancies. A zero reference uses
absolute error only; the tool never invents a relative error denominator.

## 3. V&V discipline

The tools do not define pass/fail tolerances. Acceptance criteria remain part
of each validation campaign. This prevents a generic utility from silently
weakening a qualification contract.

The same artifacts can therefore be reused for N9-S8, future N9 studies,
thermal/radiation validation, and physical qualification without embedding
NACA-specific assumptions.
