# OpenFOAM Python reader: dependency and integration decision

## Finding

CFDX already has two OpenFOAM mesh-import paths: the native C++ reader for
`constant/polyMesh` and a Python `OpenFOAMAdapter`. The Python adapter also
contains hand-written list/field parsing. It is therefore not necessary to
start another parser from scratch, but the current adapter should not be
mistaken for a complete case/results reader.

The preferred parser candidate is [foamlib](https://github.com/gerlero/foamlib).
It provides a standalone OpenFOAM file parser, dict-like access to dictionaries
and fields, ASCII/binary field support, and time-directory traversal. Its
documented API includes `FoamFile`, `FoamFieldFile`, and `FoamCase`.
[fluidfoam](https://github.com/fluiddyn/fluidfoam) is a useful alternative for
post-processing, while pyOFM is primarily a polyMesh reader.

## Licensing and runtime constraints

Do not make foamlib a mandatory dependency of the CFDX core. Keep it optional
until the repository's Python-version and licensing policies have been checked.
The current latest foamlib release declares Python >=3.12; older releases may
support earlier interpreters, but the chosen version must be resolved and tested
in CFDX's supported Python matrix. foamlib is GPL-3.0-only, so distribution and
dependency policy needs explicit review before bundling it as a required
dependency.

## Integration contract

- Keep the existing native C++ polyMesh reader as the topology path until
  parity tests show the Python path is better for a specific use case.
- Use foamlib for OpenFOAM dictionary and field parsing; do not copy its parser.
- Keep the source representation and the CFDX canonical representation separate.
- Preserve unknown dictionary entries and report unsupported conversions.
- Keep field data usable for post-processing even if CFDX cannot reproduce the
  source physics model.
- Never silently claim numerical equivalence after converting a boundary
  condition, turbulence model, discretization scheme, radiation model, or
  thermal model.

## First implementation slice

`cfdx.io.openfoam_case_reader` provides optional, read-only case inspection and
single-field reading. It is a parser integration boundary, not yet a complete
case-to-HDF5 converter. Next steps are to connect its result metadata to the
existing `ConversionResult` / gap-analysis pipeline, add versioned fixture
cases (including binary fields), and qualify the data mapping against the
native C++ mesh and the HDF5 writer.
