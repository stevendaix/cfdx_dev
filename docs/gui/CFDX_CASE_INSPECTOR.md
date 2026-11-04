# CFDX HDF5 case inspection in the Workbench

## Purpose

The Workbench can inspect a CFDX HDF5 artifact before deciding whether to load
its mesh into a renderer. Inspection is metadata-first: it reads root
attributes, group and dataset names, dataset shapes, dtypes, chunking and
compression metadata. It does not read coordinate, connectivity, cell-ID or
field-array payloads merely to open the file.

Selecting a dataset in the Case File dock reads a bounded sample (32 values
by default, hard maximum 256). The sample is read-only, uses HDF5 slicing,
preserves trailing vector components where possible, and never materializes
the full dataset. Numeric and fixed-width string datasets stored through internal HDF5 hard
links are supported; soft/external links, external raw storage, virtual,
variable-length, object and compound datasets are rejected. This is a
diagnostic sample, not a replacement for a full field reader or
post-processing pipeline.

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
not a CFDX file-format limit. Both PyVista mesh and DAT loading entry points
repeat the metadata preflight and reject oversized inputs before slicing
topology arrays, so direct renderer calls cannot bypass the Workbench guard.
The inspector records the source file device, inode, size and nanosecond
modification time; renderer entry points reject a summary if the file changed
since inspection, closing the stale-summary/time-of-check-to-time-of-use gap.
Preview eligibility also checks the topology
dataset shapes, caps eager source topology at 32 MiB and caps face connectivity
at two million entries. DAT preview additionally preflights the cell-to-face
CSR datasets, their integer dtypes, internal-link status, cell count, and the
same 32 MiB/2M-entry budgets before reading the checkpoint or connectivity.
After reading the bounded arrays, it also estimates the fully expanded VTK
polyhedron cell array in bounded chunks and refuses construction above two
million flattened entries, preventing Python-list expansion from multiplying
the source-array footprint. Immediately before any viewport mutation or VTK
object construction, the adapter validates finite
coordinates, face-offset and cell-offset CSR invariants, face vertex ranges,
patch face ranges, and cell-to-face ranges. Degenerate faces and empty cell
face lists are rejected with explicit errors. These checks defend the renderer
against malformed but superficially well-shaped HDF5 files. The limits are
deliberately conservative because the current PyVista adapter creates
additional Python-side connectivity lists. They must not be raised until
measured peak-RSS and latency benchmarks justify doing so.

The numerical checkpoint remains a separate `.dat.h5` artifact. Inspecting a
case does not load or activate a restart checkpoint.
The Case File dock displays datasets as a hierarchical HDF5 tree. Its filter
retains matching ancestor groups; dataset samples can be opened by mouse or
keyboard selection. The Copy dataset path action copies the selected absolute
HDF5 path without reading its payload. Expand all and Collapse all controls
support navigation of large metadata trees.


## Performance contract

Opening a file does not slice mesh/field datasets. Only an explicit dataset
selection reads a bounded sample. Chunked datasets whose individual chunks
exceed the 1 MiB decompression budget, virtual datasets, external storage and
non-internal HDF5 links are not sampled.
Dataset discovery stops at 300 entries by default (5,000 hard maximum) and
the HDF5 object-tree walk is capped independently at 20,000 objects. This also
bounds files containing very large numbers of empty groups. Scalar/vector field
catalogs and region names are each capped at 300 entries, and root attribute
display is bounded to 100 entries and 4 KiB per value. Truncated region
catalogs display lower-bound counts rather than presenting partial sums as
exact mesh totals. Variable-length
dataset storage is reported as unknown rather than estimated from pointer
size. Known small
format/provenance attributes and the case setup JSON are read for
identification; decoded setup JSON larger than 1 MiB is rejected before
parsing.

The current PyVista mesh loader still materializes full topology and creates
Python-side face objects. It is therefore deliberately not used for large
cases. A production large-mesh viewport requires a separate renderer workstream:
chunked HDF5 reads, C++/VTK-owned topology, surface extraction or decimation,
spatially accelerated picking, and on-demand field transfer. The correct
acceptance benchmark is at least 1M, 10M and 50M cells, recording load latency,
peak RSS, GPU memory, frame rate and selection/cut latency.

The GUI CI also runs `scripts/benchmark_cfdx_hdf5_inspector.py` at logical
sizes of 1M, 10M and 50M cells. The benchmark resolves the source package
explicitly because `scripts/cfdx.py` otherwise shadows the importable `cfdx`
package when Python executes scripts directly. It records metadata-inspection
latency, bounded-sample latency, peak RSS, logical dataset bytes, physical
sparse-file size and allocated bytes for the sampled points chunk. That one
chunk is deliberately materialized and checked against known values so the
sample measurement exercises actual HDF5 chunk I/O rather than only returning
implicit fill values. The JSON report is published as the
`cfdx-hdf5-inspector-benchmark` workflow artifact. These sparse fixtures
measure the file-browser path only; they do not simulate geometry rendering,
GPU upload or interactive picking.

## Validation

`tests/python/test_case_inspector.py` covers canonical and interchange
metadata, bounded dataset/field discovery, malformed count attributes,
bounded numeric/fixed-string samples and unsafe sample paths, plus a sparse HDF5 fixture
whose logical coordinate/connectivity arrays represent a 50-million point
case. `tests/python/test_mesh_topology_validation.py` exercises malformed CSR
offsets, out-of-range point/face references, degenerate faces and non-finite
coordinates before VTK construction. The fixture verifies metadata inspection
without materializing large coordinate/connectivity arrays.

The inspector is read-only. It does not repair files, rewrite attributes,
change case configuration or establish numerical validation/qualification.
