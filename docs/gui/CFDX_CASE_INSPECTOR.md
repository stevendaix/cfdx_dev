# CFDX HDF5 case inspection in the Workbench

## Purpose

The Workbench can inspect a CFDX HDF5 artifact before deciding whether to load
its mesh into a renderer. Inspection is metadata-first: it reads root
attributes, group and dataset names, dataset shapes, dtypes, chunking and
compression metadata. It does not read coordinate, connectivity, cell-ID or
field-array payloads.

This distinction matters for production cases whose topology contains tens of
millions of entities. Opening a file must not implicitly allocate a second copy
of the mesh in Python or GPU memory.

## Supported layouts

The inspector recognises both layouts currently written by CFDX:

- Canonical project files: root attribute `format=CFDX`, schema metadata and
  the `/case/config` JSON dataset written by `case_io.save_case`.
- CFDX mesh interchange files: `format=CFDX-HDF5-mesh-v1`, root-level CSR
  topology datasets, optional `/fields/scalar` and `/fields/vector` groups,
  source metadata and `case_setup_json` provenance.

It also reports region-grouped layouts when the topology is stored under
`/regions/<region>`. Imported case metadata is not confused with the
canonical editable project configuration: a mesh interchange file carrying
`case_setup_json` is inspectable but is not silently loaded as a project.

## Workbench behavior

The Open action first inspects the HDF5 metadata. It then:

1. Loads the editable `CFDXSession` only when a valid canonical
   `/case/config` is present.
2. Displays mesh counts, regions, scalar/vector field names, case setup keys,
   HDF5 dataset metadata and file attributes in the Case File dock.
3. Loads the existing mesh browser and PyVista preview only for a supported
   single-region root layout with no more than 100,000 points, faces and cells.
4. Leaves larger or unsupported layouts in metadata-only mode and explains
   why the full in-memory preview was not attempted.

The 100,000-entity threshold is a conservative guard for the current renderer,
not a CFDX file-format limit. It must not be raised until measured memory and
latency benchmarks justify doing so.

The numerical checkpoint remains a separate `.dat.h5` artifact. Inspecting a
case does not load or activate a restart checkpoint.

## Performance contract

The metadata inspector must not call HDF5 dataset slicing or conversion on
mesh/field datasets. The dataset browser is bounded (300 entries by default)
so a pathological file cannot create an unbounded Qt tree.

The current PyVista mesh loader still materializes full topology and creates
Python-side face objects. It is therefore deliberately not used for large
cases. A production large-mesh viewport requires a separate renderer workstream:
chunked HDF5 reads, C++/VTK-owned topology, surface extraction or decimation,
spatially accelerated picking, and on-demand field transfer. The correct
acceptance benchmark is at least 1M, 10M and 50M cells, recording load latency,
peak RSS, GPU memory, frame rate and selection/cut latency.

## Validation

`tests/python/test_case_inspector.py` covers canonical and interchange
metadata, bounded dataset discovery, malformed files, and a sparse HDF5
fixture whose logical coordinate/connectivity arrays represent a 50-million
point case. The fixture verifies metadata inspection without materializing
those arrays.

The inspector is read-only. It does not repair files, rewrite attributes,
change case configuration or establish numerical validation/qualification.
