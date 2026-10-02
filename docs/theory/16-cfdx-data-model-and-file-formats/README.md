> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 16 — CFDX Data Model, Mesh Files and Case Formats

**Status: CODE-AUDITED against the current PR branch.**

## 1. File roles

The intended distinction is:

    case.cfdx.h5  = problem definition / mesh / setup / metadata
    case.dat.h5   = numerical checkpoint / restart state
    case_<time>.vtu = visualisation artifact

The invariant is

$$
\boxed{\text{case definition}\neq\text{numerical state}}.
$$

This chapter documents the implemented case-level HDF5 path. It does not claim that every future DAT/restart feature is already implemented.

## 2. Current HDF5 schema v1

The implementation defines:

    CFDX_HDF5_FORMAT_VERSION = 1
    CFDX_HDF5_SCHEMA_VERSION = 1
    CFDX_VERSION = "0.7"

The format and schema versions are stored as root string attributes. They are distinct from the application version.

### Required root datasets

| Dataset | HDF5 type | Shape | Meaning |
|---|---|---|---|
| points | float64 | [n_points, 3] | point coordinates |
| face_vertices | uint64 | [n_face_vertex_refs] | flattened polygon connectivity |
| face_offsets | uint64 | [n_faces+1] | CSR face ranges |
| owner | uint64 | [n_faces] | owner cell per face |
| neighbour | int64 | [n_faces] | neighbour cell, -1 for boundary |
| cell_faces | uint64 | [n_cell_face_refs] | flattened cell-face connectivity |
| cell_offsets | uint64 | [n_cells+1] | CSR cell ranges |

The reader checks rank, integer type/sign, terminal offsets, index ranges and owner/neighbour dimensions.

### Field groups

    /fields/scalar/<name>  float64 [n_cells]
    /fields/vector/<name>  float64 [n_cells, 3]

The case-level writer creates these datasets directly. They represent cell fields in the current implementation.

## 3. CSR invariants

For face connectivity:

$$
0=o_0\le o_1\le\cdots\le o_{n_f},
$$

with

$$
o_{n_f}=|face\_vertices|.
$$

Face f occupies [o_f,o_{f+1}).

For cells the equivalent invariants apply to cell_offsets and cell_faces.

The reader rejects non-monotone offsets and offsets beyond the flattened storage.

## 4. Owner/neighbour semantics

Internal faces have owner >= 0 and neighbour >= 0.

Boundary faces use neighbour = -1.

The reader rejects neighbour values below -1 and indices outside the cell range.

## 5. Boundary patch extension

Boundary metadata is optional. When present, the current writer stores:

    boundary_patches      root string attribute
    patch_face_ids        uint64 [n_patch_face_refs]
    patch_face_offsets    uint64 [n_patches+1]

The attribute encodes patch name, starting face, count and numeric patch type. The reader checks metadata count/ranges against the CSR representation.

If boundary metadata is absent, the current reader reconstructs a conservative generic boundary patch from neighbour = -1. This compatibility behaviour must not be confused with preservation of original patch names/types.

## 6. Integrity metadata

The writer computes FNV-1a hashes for:
- topology;
- geometry;
- complete mesh.

The reader validates these when the complete integrity metadata set is present. Older files without all hashes remain readable under the compatibility path.

This is integrity checking, not cryptographic provenance.

## 7. Provenance and setup

Current root attributes include format/schema/application version, dates, dimension, precision, endian, source solver/format/version, source case path/name, case_setup_json and gap_report_json. The implementation uses a small JSON parser for embedded setup/gap objects.

The case setup therefore bridges:

    setup model → case_setup_json → HDF5 → C++ CaseSetup → solver configuration

## 8. Reader acceptance rules

A case is rejected for missing required topology datasets, invalid ranks/types, inconsistent dimensions, invalid CSR offsets, out-of-range indices, invalid owner/neighbour values or inconsistent boundary metadata.

These checks are part of the format contract and should have explicit regression tests.

## 9. Import pipeline

Current native families include OpenFOAM and Gmsh, with generic mesh import and VTU output paths. Import must preserve or explicitly report:
- topology;
- orientation;
- boundary groups;
- dimensional interpretation;
- unsupported constructs;
- source provenance.

## 10. Compatibility rules still to close

The code-audit identifies future schema work rather than claiming it complete:
- machine-readable schema validation;
- explicit units/dimension metadata per stored field;
- deterministic serialisation policy;
- documented chunking/compression policy;
- complete DAT/restart schema;
- N-to-M restart compatibility;
- formal schema fixtures for every reader/writer dataset.

## 11. Source of truth

The authoritative implementation files are:

    src/cfdx/io/hdf5/schema.h
    src/cfdx/io/hdf5/case_hdf5_io.h
    src/cfdx/io/hdf5/case_hdf5_io.cpp
    src/cfdx/io/hdf5/hdf5_reader.cpp
    src/cfdx/io/hdf5/hdf5_writer.cpp
    src/cfdx/io/hdf5/mini_json.h

The documentation must be updated whenever the serialized contract changes.


## Scientific explanation standard

For every important equation, algorithm or data-model rule, document the motivation; definitions/units/sign conventions; derivation or formal rationale; discrete/FVM representation where applicable; actual CFDX execution path; example; error/limitation/sensitivity analysis; exact implementation files; executable verification; benchmark/V&V evidence; and stable bibliography/reference identifiers. Tables summarize explanations and do not replace them. Keep **Implemented / Verified / Validated / Qualified** distinct; documentation maturity never promotes numerical maturity.