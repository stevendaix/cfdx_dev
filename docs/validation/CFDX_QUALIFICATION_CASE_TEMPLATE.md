# CFDX Qualification Case Sheet — Template

## Case identity

- ID:
- Domain:
- Case:
- Status:
- Implementation PR:
- Validation issue:
- Date / revision:

## Reference

- Reference type: analytical / DNS / experiment / literature / independent code
- Citation:
- Dataset/version:
- Provenance:
- Reference configuration equivalence notes:

## Parameters

| Quantity | Symbol | Value | Unit |
|---|---|---:|---|
| | | | |

Dimensionless parameters:

| Parameter | Definition | Value |
|---|---|---:|
| | | |

## Mesh

| Level | Cells/nodes | Topology | Quality | Wall resolution | Refinement ratio |
|---|---:|---|---|---|---:|
| Coarse | | | | | |
| Medium | | | | | |
| Fine | | | | | |

Mesh generation command and exact version:

## Boundary conditions

| Boundary | Physics | Value / law | Units |
|---|---|---|---|
| | | | |

Pressure reference / gauge:

## Numerical schemes

- Gradient:
- Interpolation:
- Convection:
- Diffusion:
- Pressure equation:
- Temporal:
- Turbulence:
- Thermal:
- Radiation angular/spatial:
- Other:

## Solver settings

- Algorithm:
- Linear solver:
- Preconditioner:
- Linear tolerance:
- Linear max iterations:
- Nonlinear max iterations:
- Relaxation:
- Initialization:
- Parallelism:
- Precision:

## Convergence criteria

- Momentum:
- Continuity:
- Pressure:
- Energy:
- Turbulence:
- Radiation:
- QoI stabilization:
- Statistical stationarity:

## Conservation criteria

- Mass imbalance:
- Momentum imbalance:
- Energy imbalance:
- CHT interface heat-flux mismatch:
- Radiation closure:
- Other:

## Quantity of interest

For every QoI define the mathematical extraction formula, sampling location and units.

| QoI | Definition | Units | Extraction |
|---|---|---|---|
| | | | |

## Reference value

| QoI | Reference | Source |
|---|---:|---|
| | | |

## CFDX result

| Mesh | QoI | CFDX | Absolute error | Relative error |
|---|---|---:|---:|---:|
| | | | | |

## Observed order

For error values E_i on systematically refined meshes with ratio r:

p_obs = ln(E_coarse / E_fine) / ln(r)

For three or more levels, report the full sequence and use GCI where appropriate.

| Pair / level | Error | Refinement ratio | Observed order | GCI |
|---|---:|---:|---:|---:|
| | | | | |

## Runtime

| Mesh | Mesh generation | Solver | Iterations | CPU/runner | Parallelism |
|---|---:|---:|---:|---|---:|
| | | | | | |

## Regression status

- Executable test:
- CTest name:
- CI workflow:
- Artifact/log:
- Determinism:
- Regression acceptance:
- Current status: PLANNED / READY / RUNNING / DIAGNOSTIC / PASS / BLOCKED

## Audit conclusion

State separately:

1. code verification result;
2. solution verification result;
3. physical validation result;
4. remaining uncertainty;
5. exact reason for promotion or non-promotion.
