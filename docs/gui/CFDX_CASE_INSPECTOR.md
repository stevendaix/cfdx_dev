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

The Case File dock also provides an explicit **Compute statistics** action for
the selected numeric dataset. It scans in bounded blocks and reports minimum,
maximum, arithmetic mean, RMS, finite/non-finite counts and bytes scanned.
The default logical-read budget is 64 MiB and the default block/decompression
budget is 1 MiB; larger datasets or chunks are refused with a reason instead
of silently allocating more memory. Multidimensional fields are traversed by
first-axis slabs, and a single row wider than the block budget is refused. This quick-look treats every stored value
equally: it is not volume-, mass- or area-weighted and must not be interpreted
as a conservative CFD integral. Opening a case and selecting a dataset still
do not trigger this full-dataset scan; the user must explicitly request it.

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
The dataset tree also exposes a read-only storage audit for every discovered dataset:
logical payload size (shape × dtype), allocated HDF5 storage bytes, allocated/logical
ratio, chunk shape and compression filter. This distinguishes dense datasets from
chunked sparse datasets whose logical dimensions can be enormous while few chunks
are physically allocated. The ratio is not a compression benchmark: it reflects
allocated HDF5 storage for that dataset and can be zero for an unallocated chunked
dataset. Variable-length/object datasets retain unknown logical/allocated sizes
rather than reporting misleading pointer-based estimates. These values are read
from HDF5 metadata and do not materialize field or mesh arrays.

The Case File dock displays datasets as a hierarchical HDF5 tree. Its filter
retains matching ancestor groups; dataset samples can be opened by mouse or
keyboard selection. The Copy dataset path action copies the selected absolute
HDF5 path without reading its payload. Expand all and Collapse all controls
support navigation of large metadata trees.

The Mesh dock creates face rows in pages of 250 when a patch is expanded; a
large patch therefore does not create one Qt item per face on the first click.
Picking a rendered boundary patch selects its canonical patch ID and synchronizes
the tree, application selection and 3D highlight (right-click a boundary
patch in the viewport). Picking individual VTK cells or CFDX faces is
intentionally not inferred from renderer cell order: a stable
VTK-cell-to-CFDX-face mapping is still required before that interaction can be
offered safely.


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
`cfdx-hdf5-inspector-benchmark` workflow artifact. The JSON also includes
streaming statistics for a separate materialized 4,096-value scalar probe,
including extrema, finite count, bytes scanned and elapsed time. This exercises
the statistics path without scanning the sparse production-scale field or
confusing a synthetic probe with a full-case result. These sparse fixtures
measure the file-browser path only; they do not simulate geometry rendering,
GPU upload or interactive picking.

The generic `scripts/n15_benchmark.py` runner measures the direct child
process's peak RSS independently on POSIX systems using `wait4`; this avoids
reporting the cumulative high-water mark of earlier benchmark children. On
platforms without `wait4`, `peak_rss_bytes` is `null` with source
`unavailable` unless a metrics file supplies a measurement. This process-level
measurement does not include descendant processes and is not a GPU-memory
measurement.

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


## Mesh Intelligence: surface geometry diagnostics

The Case File dock exposes **Analyze surface quality** only when the metadata
preflight confirms the supported, bounded single-region surface topology. The
action checks the inspected file identity again and reads only the root
`points`, `face_vertices` and `face_offsets` datasets after validating
internal hard links, numeric dtypes, dimensions, entity/connectivity limits,
the 32 MiB source-data budget and the 1 MiB per-chunk decompression budget.

The current report contains:
- Polygon area minimum and maximum.
- Face edge-length ratio (longest / shortest edge), with min, p95 and max.
- Relative warpage: maximum vertex distance to the centroid plane normal to the
  polygon area vector, divided by the longest edge; min/p95/max summary and
  worst face IDs.
- Degenerate faces and faces with numerically zero-length edges; a zero-length edge makes a face degenerate even when its polygon area is otherwise nonzero.

These are **surface-face geometric indicators**, not volume-cell quality
metrics. A surface-only boundary mesh cannot establish cell skewness,
non-orthogonality, Jacobian, volume, cell aspect ratio or interior-cell
validity. Those require the complete cell polyhedra, cell centres and
neighbour relationships; the UI explicitly says so rather than inferring them
from boundary faces. The current implementation is an on-demand diagnostic,
not a mesh repair operation and it never modifies the source file.


The GUI workflow also runs a deterministic 128 × 128 warped-quad probe and
publishes `cfdx-surface-quality-report` as a JSON artifact. It verifies that
all 16,384 faces are processed and that the known warped surface produces a
nonzero warpage indicator. Its elapsed time is a CI smoke measurement, not a
cross-machine performance claim.
