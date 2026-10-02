# Executable N2 reference experiments

These scripts are deliberately **reference experiments**, not CFDX solver code.

## Reproducibility contract

Each experiment must state:

- the analytic field and exact solution;
- the mesh family and refinement parameter;
- the reconstruction algorithm;
- boundary/stencil policy;
- error norm and weighting;
- observed-order calculation;
- conditioning diagnostics where relevant;
- whether the result is demonstration data or qualification evidence.

The current reference study uses a smooth skewed quadrilateral mesh and evaluates
interior cells with complete four-neighbour stencils. This isolates spatial
reconstruction from boundary closure.

## Experiments

- `linear_field.py` — manufactured linear field and exact constant gradient.
- `quadratic_field.py` — manufactured quadratic field and exact gradient.
- `mesh.py` — deterministic structured-but-skewed quadrilateral mesh generator.
- `reconstruction.py` — independent Green–Gauss, LS and WLS reference algorithms.
- `metrics.py` — volume-weighted relative L2 error.
- `convergence.py` — observed-order calculation.
- `run_gradient_study.py` — complete refinement study producing JSON.
- `self_check.py` — deterministic linear-field exactness and conditioning checks.

The study is intentionally independent from the CFDX C++ implementation. It is
therefore suitable as an external mathematical reference, but it is not by
itself a qualification oracle.

Run locally with:

    python self_check.py
    python run_gradient_study.py

The expected scientific result is determined by the computed data. No
hard-coded convergence values are embedded in the scripts.
