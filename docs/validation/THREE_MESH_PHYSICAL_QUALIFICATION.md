# CFDX — three-mesh physical qualification (new opt-in tool)

**Scope:** `scripts/qualify_mesh_series.py` is a fail-closed **evidence reader**, not
an alternative CFD solver, a numerical-model repair, or a source of physical
reference values. Its synthetic unit tests must never be used to promote any
real solver to `PASS`. Existing CI, qualification registries, solver tolerances,
and PR acceptance gates are unchanged.

## Purpose and limitations

- Evaluate three **distinct, completed and converged** runs generated at the
  **same full Git commit**, with SHA-256-verified mesh and machine-readable
  result artifacts. A ZIP without Git provenance **cannot be promoted**.
- Enforce residual and local/global conservation separately. Solving a linear
  system is not equivalent to conservation of mass or energy.
- Compute an observed order and a normalized fine-mesh GCI with the **actual**
  two mesh-refinement ratios (including small non-uniformities), and compare
  each QoI against an **independently sourced reference**.
- In `naca_laminar` / `naca_sst` mode, require distinct oracles for pressure
  drag, viscous drag, total drag and lift, and check drag decomposition on
  **every** mesh. A converged Cd_total alone cannot qualify wall shear.
- In `thermal` mode, require local and global energy balances (normalized as
  defined by the simulation exporter) and physical heat/temperature QoIs.
- In `sst` mode, the caller chooses independently checked quantities such as
  `Cf`, velocity profile errors, or displacement thickness. **This script by
  itself does not qualify SST**: the accompanying campaigns must still gate
  positivity of k/omega, wall distance/y+, turbulence production/dissipation,
  domain extent, inflow conditions and reference uncertainty.

**Not proved by the gate:** mesh connectivity quality, dimensional
consistency of quantities, validity of the exported measured residual,
credibility of an external literature reference, RANS closure or a final CFDX
production certification. These are independent campaign responsibilities.
The gate never invents physical data or changes a failed status to `PASS`.

## Input format (schema 1.0)

All numbers must be **finite** JSON numbers (`NaN` and `Infinity` are rejected).
`source_sha` and `--expected-source-sha` must be identical full 40-hex-character
Git commit hashes, not branch names or the temporary local commit from an
extracted ZIP. Use a genuine, clean source checkout for promotion.

```json
{
  "schema_version": "1.0",
  "case_id": "naca0012-re1000-alpha0",
  "physics": "naca_laminar",
  "source_sha": "<40 lowercase hex characters>",
  "criteria": {
    "max_true_residual": 1e-8,
    "max_mass_global": 1e-7,
    "max_mass_local": 1e-7,
    "max_drag_decomposition_error": 1e-10,
    "min_refinement_ratio": 1.5,
    "max_refinement_ratio_difference": 0.01,
    "min_observed_order": 1.0,
    "max_gci": 0.02
  },
  "references": {
    "cd_pressure": {"value": 0.01, "scale": 0.01, "max_relative_error": 0.05},
    "cd_viscous": {"value": 0.04, "scale": 0.04, "max_relative_error": 0.05},
    "cd_total": {"value": 0.05, "scale": 0.05, "max_relative_error": 0.05},
    "cl": {"value": 0.0, "scale": 0.01, "max_relative_error": 0.05}
  },
  "levels": [
    {"mesh": {"path": "meshes/coarse.msh", "sha256": "<64 lowercase hex characters>"},
     "result": {"path": "results/coarse.json", "sha256": "<64 lowercase hex characters>"}},
    {"mesh": {"path": "meshes/medium.msh", "sha256": "<64 lowercase hex characters>"},
     "result": {"path": "results/medium.json", "sha256": "<64 lowercase hex characters>"}},
    {"mesh": {"path": "meshes/fine.msh", "sha256": "<64 lowercase hex characters>"},
     "result": {"path": "results/fine.json", "sha256": "<64 lowercase hex characters>"}}
  ]
}
```

**Important:** the numerical reference values and thresholds shown above are
*format illustrations only, not validated NACA reference values*. They must be
replaced by independently justified values and a documented error budget.
No `PASS` may be inferred from this example. `scale` is a strictly positive,
explicit **normalizing magnitude**; do not divide by a zero lift reference.
For `thermal`, use `max_energy_global` and `max_energy_local` instead of mass
thresholds; for generic `laminar` and `sst` select relevant QoI names.

Each `result` points to a completed, **retained JSON export of a production
solver run**. The exporter must generate these values from the solver and its
faces/cells, not copy manual guesses from a configuration file:

```json
{
  "case_id": "naca0012-re1000-alpha0",
  "physics": "naca_laminar",
  "source_sha": "<same full source commit>",
  "mesh_sha256": "<SHA-256 of actual mesh input>",
  "status": "COMPLETED",
  "converged": true,
  "cells": 12336,
  "h": 0.008,
  "true_residual": 1e-9,
  "balances": {"mass_global": 1e-9, "mass_local": 5e-9},
  "qoi": {"cd_pressure": 0.01, "cd_viscous": 0.04, "cd_total": 0.05, "cl": 0.0}
}
```

All numeric outputs in this second block are **illustrative**, not measured
CFDX results. The `h` values must be physical characteristic mesh spacings
computed consistently for a fixed geometry/dimension (for example,
`h=(domain_volume/cell_count)^(1/d)` for a genuinely d-dimensional refinement
family), not hand-assigned arbitrary labels. Independent mesh-quality,
wall-gradient, far-field and domain-size checks are mandatory for real cases.

Balances must be dimensionless **absolute values** normalized against a fixed
physical flux scale, with no cancellation hiding large local imbalances:

- `mass_global`: magnitude of net boundary mass flux divided by a documented
  nonzero reference throughput. `mass_local`: largest absolute cell flux
  imbalance divided by the same reference scale (or another stated, fixed
  physically meaningful one).
- `energy_global`: magnitude of source/input minus sink/output (including
  storage for transient cases), normalized by a documented nonzero heat-flow
  scale. `energy_local`: largest cell energy imbalance, using the same fixed
  scaling convention.
- `true_residual`: norm of the **real assembled** linear residual `b-Ax`
  divided by the run's declared residual normalization. The exporter must
  retain the normalization choice and linear/nonlinear iteration history.
  This is not a substitute for nonlinear momentum, pressure, k/omega or energy
  residuals, which need separate per-physics campaign checks.

## Run

From the root of an actual Git checkout containing the retained mesh and result
artifacts, with fixed, reviewed campaign criteria:

```bash
python3 scripts/qualify_mesh_series.py \
  --input <evidence-directory>/case.json \
  --expected-source-sha "$(git rev-parse HEAD)" \
  --report <evidence-directory>/three_mesh_verdict.json
```

Return code **0** indicates that all supplied data passed the declared gates;
return code **1** indicates `FAIL` and writes the diagnostic reason to JSON.
This CLI must be run *after* the original solver campaigns and independent
reference verification. No new production qualification is claimed until
actual CFDX runs produce and retain the required files.

### Numerical convergence rule

Three sizes are ordered coarse → medium → fine. For a QoI `f` and ratios
`r0 = h0/h1`, `r1 = h1/h2`, the positive observed order `p` solves

`abs((f0-f1)/(f1-f2)) = (r0^p - 1)/(1 - r1^(-p))`.

The monitored fine-grid GCI is
`1.25 * abs(f2-f1) / ((r1^p - 1) * reference_scale)`.
The gate rejects non-monotone sequences, implausible order, inconsistent or
insufficient refinement, and unqualified references. Mesh-independent exact
values can pass if they satisfy the physical reference error gate, with order
reported symbolically as `mesh_independent`.

### Integration and PR separation

The tool is **opt-in**. Only its **synthetic self-tests** are registered under
`test_three_mesh_physical_gate_selftest` in CTest. It does **not** override any
existing readiness/physical-validation status or any required CI check.
Use independently reviewed solver exports, then associate the retained verdict
with the appropriate laminar, SST or thermal case in the qualification registry
in a separate evidence-promotion review. Do not promote NACA based on the
currently failing n128 run.
