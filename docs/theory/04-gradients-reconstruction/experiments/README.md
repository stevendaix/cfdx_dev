# Executable Experiments

These experiments are part of the Theory chapter. They are intentionally small and independently testable.

## Reproducibility contract

Every numerical experiment should record source revision, Python version, dependencies, mesh definition, field definition, numerical method, boundary policy, limiter state, convergence criteria, error norm and raw results.

Generated plots and tables must be reproducible from these sources.

## Current experiments

- `linear_field.py`: analytic linear field and exact constant gradient.
- `quadratic_field.py`: analytic quadratic field and exact spatially varying gradient.
- `convergence.py`: observed-order calculation.
- `metrics.py`: volume-weighted relative L2 metric.
- `make_report.py`: machine-readable convergence report.
