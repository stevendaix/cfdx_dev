# SU2 NACA0012 validation case — reader qualification

## Purpose

This case is the external-solver I/O qualification case for CFDX.

It replaces the axisymmetric sphere case as the **format-reader validation
target** for the unified solver-I/O path. It does **not** remove or invalidate
the native CFDX axisymmetric VMFL036 capability; that remains a separate
physics/solver qualification.

## Source case

The case is the SU2 NACA0012 Euler validation/tutorial family:

- 2-D NACA0012 external flow;
- Euler solver;
- Mach number 0.8;
- angle of attack 1.25 degrees;
- freestream pressure 101325 Pa;
- freestream temperature 288.15 K;
- SU2 mesh format;
- two boundary markers: `airfoil` and `farfield`.

The repository fixture is the existing SU2 NACA0012 mesh/config pair under
`tests/data/su2/`. The mesh contains 5233 vertices and 10216 triangular
cells.

SU2's current regression suite also contains a NACA0012 Euler regression
case, with documented iteration/QoI reference values. The CFDX PR intentionally
uses the SU2 input files as an **input-format/readability qualification**, not
as evidence that CFDX reproduces SU2's solver result.

## What the PR verifies

The test reconstructs a conventional SU2 case directory:

```text
naca0012/
├── mesh.su2
└── config.cfg
```

and calls the same directory-level `Su2Adapter.convert()` entry point used
by the unified importer.

It verifies:

1. SU2 mesh dimensionality;
2. vertex and cell cardinalities;
3. triangular cell connectivity;
4. boundary-marker discovery;
5. SU2 configuration parsing;
6. Euler physics identification;
7. freestream pressure/temperature;
8. iteration count;
9. conversion without blocking Gap Analysis findings.

## Results are deliberately not used

The repository's small `mismatched_solution.csv` fixture contains only 10 rows
while the NACA0012 mesh contains 10216 cells. It is therefore not a valid
solution for this mesh, and it is deliberately kept out of the mesh directory
so the adapter does not auto-discover it and block every conversion of this
case.

The test does **not** weaken the cardinality check and does not use that
fixture as validation evidence. Missing solution output is recorded as a
non-blocking diagnostic because this PR validates case reading, not result
import.

A future result-import qualification should use an actual SU2 solution whose
field localisation/cardinality matches the mesh.

## Acceptance

The reader qualification passes only when:

- the real SU2 mesh is parsed;
- the real SU2 configuration is parsed;
- mesh/config metadata are internally coherent;
- no blocking Gap Analysis finding is produced.

No CFDX numerical result is compared here. Numerical validation remains a
separate solver-level campaign.
