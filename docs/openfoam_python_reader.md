# OpenFOAM Python reader: native parser and integration limits

## Decision

CFDX uses a small native Python parser for the common ASCII OpenFOAM syntax
needed by case inspection, setup dictionaries, scalar/vector fields, and
polyMesh lists. It has no runtime dependency on an external OpenFOAM
Python package. This avoids a GPL-3.0-only dependency and preserves support for
the project's Python 3.10 baseline.

The parser is deliberately scoped: it handles ordinary dictionaries, nested
dictionary blocks, quoted strings, comments, numeric values, parenthesized
lists, dimension bracket lists, uniform fields, and nonuniform scalar/vector
volume fields. It is not a complete implementation of the OpenFOAM file
grammar.

## Supported import path

- The adapter reads ASCII constant/polyMesh/points, faces, owner, neighbour,
  and boundary files.
- Native face/owner/neighbour topology is preserved as CFDX CSR. Surface
  triangles exported to mesh.vtk are for visualization only and are not
  treated as volume cells.
- Numeric time directories are sorted numerically; --results imports
  cell-centred scalar/vector volume fields from the latest time.
- Common setup dictionaries are parsed and preserved in case metadata.
  Supported controls are mapped explicitly; fvSchemes and fvSolution
  remain source metadata until equivalent CFDX numerical registry mappings
  are qualified.
- Unparsed setup dictionaries are preserved as raw text and reported as a
  non-blocking gap. Requested result imports fail clearly when field syntax
  cannot be decoded; values are never silently coerced.

## Explicit limitations

- Binary polyMesh and binary field files are not supported by this native
  Python parser. They are detected and reported rather than being treated as
  ASCII. The original case remains unchanged.
- Tensor, surface, point, and decomposed processor fields are not imported as
  cell-centred scalar/vector arrays.
- Boundary-condition values, material properties, turbulence/radiation
  closure equivalence, and all numerical scheme/solver controls are not yet
  translated comprehensively into CFDX runtime configuration.
- The parser is not intended to execute OpenFOAM directives such as arbitrary
  #include/macro expansion or reproduce every release-specific grammar
  extension.

## Qualification contract

Keep source metadata separate from the CFDX canonical representation. Never
claim numerical equivalence merely because a case was parsed. Qualify
representative ASCII cases against source mesh/fields and the resulting
case.cfdx.h5, including cell counts, face ownership, boundary patch ranges,
field locations, dimensions, and numeric values. Binary support requires a
separate implementation and fixtures.

## Running CFDX from an imported result field

The production incompressible solver supports an explicit, opt-in velocity
initialization from the cell-centred OpenFOAM field stored at
`/fields/vector/U`:

```bash
cfdx-io convert /path/to/openfoam-case --results --output /path/to/converted
cfdx_production_solver --mesh /path/to/converted/case.cfdx.h5 \
  --output-dir /path/to/run --iterations 500 --initial-velocity-field
```

The option is intentionally explicit; without it, existing CFDX initialization
behavior is unchanged. It requires a CFDX case HDF5 file containing a vector
field named `U`, with exactly one finite three-component value per cell.
Missing fields, wrong cell counts, non-finite values, or a non-case mesh input
fail clearly rather than falling back to a uniform velocity.

Only velocity is seeded from the imported result in this first integration
step. Pressure remains the case-declared initial pressure (or zero when none is
declared). OpenFOAM pressure is not imported as the production pressure initial
guess because common incompressible cases store kinematic pressure while CFDX's
production field is declared in pascals; the current interchange metadata does
not yet provide a qualified dimensional conversion. The imported field is an
initial guess, not a claim that OpenFOAM boundary conditions, turbulence,
thermal/radiation models, or discretization schemes are numerically equivalent.

## First real-case qualification: icoFoam lid-driven cavity

Use the official OpenFOAM **icoFoam/cavity** tutorial as the first real-case
qualification. This is a transient, incompressible, laminar, Newtonian,
single-phase flow solved with PISO. It exercises standard polyMesh topology
and cell-centred `U` and kinematic pressure `p`. It is not a VOF case.

The OpenCFD tutorial guide lists the path as
`$FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity`; layouts may vary by
distribution/version. Record the OpenFOAM distribution/version and actual
source path in the qualification evidence. See the [official lid-driven cavity
tutorial](https://www.openfoam.com/documentation/tutorial-guide/2-incompressible-flow/2.1-lid-driven-cavity-flow).

### Reproducible end-to-end reader procedure

With the OpenFOAM environment sourced and CFDX Python I/O installed, run from
the CFDX repository root:

```bash
bash scripts/qualify_openfoam_icofoam_cavity.sh /tmp/cfdx-icofoam-qualification
```

The script locates the installed tutorial instead of assuming a version-specific
path, copies it to a clean working directory, runs the tutorial's `Allrun`
when available (otherwise `blockMesh` and `icoFoam`), selects the latest
numeric time containing both `U` and `p`, converts it through `cfdx-io`,
and checks mesh counts, topology array sizes/index ranges, field shapes, and
finite values. It writes provenance to `provenance.txt` and machine-readable
reader evidence to `evidence/reader_import_evidence.json`. Use a fresh output
directory for each run.

This script is a reproducible structural/import check, not the whole solver
qualification: it does not yet compare every imported value bit-for-bit against
the OpenFOAM source, validate patch-range semantics independently, or execute
the CFDX production solver. Do not mark the issue complete based on this script
alone. Keep generated cases and evidence outside the repository unless
redistribution/licensing has been checked.

### Acceptance checks

For full qualification, retain evidence for:

- OpenFOAM distribution/version, tutorial provenance, and selected time;
- points, faces, cells, owner/neighbour topology, and boundary patch ranges;
- imported `U` dimensions `[0 1 -1 0 0 0 0]`, cell-centred location,
  shape `(n_cells, 3)`, finite values, and exact comparison with source values;
- imported `p` dimensions `[0 2 -2 0 0 0 0]` if it is imported as a scalar;
  do not map it to CFDX pressure in Pa without qualified dimensional conversion;
- exact preservation of imported `U` after production-solver initialization
  and before the first residual.

The tutorial's mesh and time-directory contents depend on the installed
OpenFOAM distribution/version and run configuration. Record actual values;
never hard-code unverified counts or claim a run passed without retained output.

## Required-vs-optional production solver check

Ordinary Python runs may skip solver-backed tests when
`CFDX_PRODUCTION_SOLVER` is absent. A dedicated qualification must fail rather
than skip in that situation. Configure both variables and run:

```bash
export CFDX_PRODUCTION_SOLVER=/absolute/path/to/cfdx_production_solver
export CFDX_REQUIRE_OPENFOAM_SOLVER_QUALIFICATION=1
pytest -q tests/python/test_openfoam_results_cli.py
```

This validates the fail-closed configuration and, when the solver is available,
the post-initialization velocity evidence. It does not establish equivalence
with OpenFOAM. The current integration initializes only `U`; pressure,
boundary-condition mapping, discretization equivalence, turbulence,
thermal/radiation, and VOF remain outside this qualification.
