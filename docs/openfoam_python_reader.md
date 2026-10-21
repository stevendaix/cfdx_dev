# OpenFOAM Python reader: native parser and integration limits

## Decision

CFDX uses a small native Python parser for the common ASCII OpenFOAM syntax
needed by case inspection, setup dictionaries, scalar/vector fields, and
polyMesh lists. It has no runtime dependency on foamlib or another OpenFOAM
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
