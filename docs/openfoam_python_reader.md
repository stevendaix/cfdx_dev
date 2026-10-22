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
claim numerical equivalence merely because a case was parsed. Before relying on
the importer for production cases, qualify representative ASCII cases against
the OpenFOAM source mesh/fields and the resulting case.cfdx.h5, including
cell counts, face ownership, boundary patch ranges, field locations, and
numeric values. Binary support requires a separate implementation and fixtures.


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


## First real-case qualification: VOF dam-break

For VOF qualification, use the official OpenFOAM **damBreak** tutorial with
`interFoam`, not the `icoFoam/cavity` tutorial. `icoFoam` is a transient,
incompressible, laminar, single-phase solver; it has no transported phase
fraction. `interFoam` solves two incompressible, isothermal immiscible fluids
using VOF and therefore exercises the phase-fraction field `alpha.water`
(or the phase name used by that OpenFOAM distribution), in addition to
velocity and pressure fields.

The canonical tutorial is commonly located at
`$FOAM_TUTORIALS/multiphase/interFoam/laminar/damBreak/damBreak` in OpenCFD
distributions. Tutorial directory layouts vary by distribution/version; locate
the installed tutorial instead of assuming the path exists. Record the exact
OpenFOAM distribution/version and source path in the evidence.

### Reproducible acquisition and conversion

Run in a shell with the OpenFOAM environment sourced:

```bash
printf 'OpenFOAM version: '; foamVersion
printf 'Tutorial root: %s\\n' "$FOAM_TUTORIALS"
find "$FOAM_TUTORIALS" -type d -path '*/interFoam/*/damBreak' -print
```

Select the installed laminar dam-break tutorial path found above, copy it to a
working directory, and run its documented meshing and solver commands
(e.g. `blockMesh` and `setFields`, then `interFoam`, where required by that
tutorial's README/scripts). Do not assume the exact command sequence across
OpenFOAM distributions. Then run:

```bash
cfdx-io convert ./openfoam-dambreak --results --output ./converted
```

This is a **VOF reader qualification**, so the acceptance checks must explicitly
include the phase-fraction field and not only `U`: confirm that the importer
recognizes the scalar volume field `alpha.water` (or the actual phase field
name), preserves its cell-centred values and dimensions, and checks values are
finite and lie within the expected physical interval [0, 1] up to the source
solver's documented numerical tolerances. Also verify mesh counts, owner/
neighbour consistency, boundary patch ranges, and selected time. The current
reader's generic scalar/vector field support does not itself establish VOF
solver equivalence or validate alpha transport.

Keep the original tutorial's provenance and license information. Do not commit
generated solver output by default; only add a small redistributable fixture
after checking upstream licensing and preserving attribution.

### Required-vs-optional production solver check

The production-solver integration tests skip when
`CFDX_PRODUCTION_SOLVER` is absent during ordinary Python test runs. A
dedicated qualification run must not count that skip as a pass. Set both
environment variables in that job:

```bash
export CFDX_PRODUCTION_SOLVER=/absolute/path/to/cfdx_production_solver
export CFDX_REQUIRE_OPENFOAM_SOLVER_QUALIFICATION=1
pytest -q tests/python/test_openfoam_results_cli.py
```

With the strict flag enabled, a missing solver executable configuration raises
a test failure instead of skipping. This fail-closed check prevents a false
green qualification status; it does not itself qualify VOF. The current
production integration flag `--initial-velocity-field` initializes only `U`.
It does not initialize or solve `alpha.water`, and must not be presented as
end-to-end VOF support. Issue #1048 should therefore remain open until CFDX has
a qualified alpha-field import path and the intended solver-side VOF workflow
is explicitly tested.
