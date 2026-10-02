# 16 — CFDX Data Model, Mesh Files and Case Formats

**Status: IN PROGRESS — repository-grounded format chapter.**

The numerical method cannot be understood independently of the discrete data on which it operates. This chapter explains topology, geometry, import, HDF5 storage and the distinction between a case and a restart.

## 1. Three file roles

~~~text
case.cfdx.h5
   mesh + setup + physics + numerics + metadata
             |
             v
       solver initial state
             |
             v
case.dat.h5
   checkpoint / restart state
             |
             v
case_<time>.vtu
   visualization / post-processing
~~~

The fundamental invariant is:

$$
\boxed{\text{case definition} \neq \text{numerical state}.}
$$

A CFDX case is the source of truth for the problem definition. A DAT is a numerical checkpoint. A VTU is an output artifact.

## 2. Topology versus geometry

Topology answers which points form a face, which cell owns a face, which cell neighbours a face, and which faces belong to a patch.

Geometry answers point coordinates, face centroid, face area, face area vector, cell centroid, cell volume and mesh-quality indicators.

~~~text
points + connectivity
        |
        v
     topology
        |
        v
 face/cell geometry
        |
        v
     FVM operators
~~~

This distinction is fundamental for reproducibility.

## 3. Actual CFDX-HDF5 v1 layout

The implementation does **not** currently use the previously proposed hierarchical `/mesh/...`, `/physics/...` and `/numerics/...` tree as its on-disk representation. The current case-level implementation is deliberately flatter.

The authoritative root-level datasets and groups are documented in `src/cfdx/io/hdf5/case_hdf5_io.h`:

~~~text
/
├── points                         float64 [n_points, 3]
├── face_vertices                  uint64  [n_face_vertex_refs]
├── face_offsets                   uint64  [n_faces + 1]
├── owner                          uint64  [n_faces]
├── neighbour                      int64   [n_faces]
├── cell_faces                     uint64  [n_cell_face_refs]
├── cell_offsets                   uint64  [n_cells + 1]
│
├── /fields/
│   ├── /scalar/<name>             float64 [n_cells]
│   └── /vector/<name>             float64 [n_cells, dim]
│
└── root attributes
    ├── format_version
    ├── schema_version
    ├── cfdx_version
    ├── topology_hash
    ├── geometry_hash
    ├── mesh_hash
    ├── creation_date
    ├── modification_date
    ├── dimension
    ├── precision
    ├── endian
    ├── source_solver
    ├── source_format
    ├── source_version
    ├── source_case_path
    ├── source_case_name
    ├── case_setup_json
    ├── gap_report_json
    ├── mesh_topology
    └── boundary_patches (optional)
~~~

This distinction is important: the Theory documentation must never present a **target architecture** as though it were the current file format.

The current C++ implementation defines:

~~~text
CFDX_HDF5_FORMAT_VERSION = 1
CFDX_HDF5_SCHEMA_VERSION = 1
CFDX_VERSION = "0.7"
~~~

The format therefore has three different versioning concepts:

1. **format version** — file-level interchange contract;
2. **schema version** — interpretation of the stored datasets/attributes;
3. **CFDX application version** — software version producing/consuming the file.

A change to the serialized layout must update the appropriate compatibility contract and its tests.

## 4. Schema metadata

The implementation currently defines:

~~~text
CFDX_HDF5_FORMAT_VERSION = 1
CFDX_HDF5_SCHEMA_VERSION = 1
CFDX_VERSION = "0.7"
~~~

The schema also carries provenance such as source solver/format/version, source case information, topology identifiers and hashes.

The Theory documentation must distinguish format version, schema version, application version and physical case version.

## 5. Variable-length connectivity

CFDX uses a CSR-like representation for polygonal connectivity:

~~~text
face_vertices
face_vertices_offsets
~~~

For example:

~~~text
face_vertices =
[0,1,2,3,4,5,6]

face_vertices_offsets =
[0,4,7]
~~~

Face i occupies

$$
[offset_i,offset_{i+1}).
$$

This permits arbitrary polygonal faces without fixed-size padding.

## 6. Owner/neighbour orientation

An internal face has owner >= 0 and neighbour >= 0.

A boundary face has neighbour = -1.

The orientation defines the sign of the face area vector. For an internal face:

$$
\mathbf S_{f,P}=-\mathbf S_{f,N}.
$$

This is a mathematical conservation invariant, not merely an implementation convention.

## 7. Boundary patches

A patch identifies a geometric set of boundary faces.

Current semantic roles include inlet, outlet, wall, symmetry, periodic and interface.

A patch role is not itself a scalar boundary condition.

~~~text
geometric patch
      |
      v
boundary role
      |
      v
field condition
      |
      v
mathematical constraint
      |
      v
FVM boundary contribution
~~~

The Theory chapter must derive how Dirichlet, Neumann, mixed and coupled conditions enter the control-volume equation.

## 8. Case setup serialization

The current HDF5 implementation also carries a serialized case_setup_json attribute.

The bridge is:

~~~text
Python setup model
      |
      v
case_setup_json
      |
      v
HDF5 case
      |
      v
C++ CaseSetup
      |
      v
solver configuration
~~~

The implementation contains a lightweight JSON parser because the setup is embedded as an HDF5 attribute.

## 9. Import pipeline

Current import families include OpenFOAM and Gmsh native paths plus Python/meshio adapters.

~~~text
OpenFOAM polyMesh ----Gmsh .msh -------------+--> import / mapping
meshio adapters -------/
                         |
                         v
                   CFDX topology
                         |
                         v
                   geometry build
                         |
                         v
                  mesh validation
                         |
                         v
                 hashes/provenance
                         |
                         v
                   case.cfdx.h5
~~~

Import is a data-model transformation. Each importer must preserve connectivity, orientation, boundary groups and provenance, or report what is lost.

## 10. Mesh validation chain

~~~text
read file
   |
parse topology
   |
construct Mesh
   |
validate indices/connectivity
   |
construct boundary patches
   |
compute geometry
   |
check volumes/orientation
   |
mesh quality
   |
hash/provenance
   |
accept numerical case
~~~

The solver must not silently repair an invalid mesh by changing the mathematical discretisation.

## 11. Restart semantics

A DAT may contain physical time, nonlinear iteration, solution fields, temporal history, global cell identifiers and restart metadata.

It must not silently redefine the physical model, boundary conditions, numerical method selection or mesh topology.

Case/DAT compatibility must be checked explicitly.

## 12. Implementation traceability

Principal current files include:

- src/cfdx/io/hdf5/schema.h
- src/cfdx/io/hdf5/case_hdf5_io.h
- src/cfdx/io/hdf5/case_hdf5_io.cpp
- src/cfdx/io/hdf5/hdf5_reader.h/.cpp
- src/cfdx/io/hdf5/hdf5_writer.h/.cpp
- src/cfdx/io/hdf5/mini_json.h
- src/cfdx/python/cfdx/io/hdf5_writer.py
- src/cfdx/io/openfoam/openfoam_importer.cpp
- src/cfdx/io/gmsh/gmsh_importer.cpp
- src/cfdx/io/mesh/mesh_importer.cpp
- src/cfdx/io/vtu/vtu_writer.cpp

Tests of the reader/writer and conversion pipeline are part of the evidence for the format.

## 13. Format chapter standard

For every supported input format the final chapter will document grammar, node representation, element representation, cell mapping, face extraction, boundary semantics, orientation, dimensionality, units, unsupported constructs, conversion loss, post-import validation, CFDX destination datasets, exact importer implementation and tests.

## 14. Future improvements

Possible improvements are tracked separately from current support:

- schema machine validation;
- explicit dtype/shape contracts;
- deterministic serialization;
- stronger topology/geometry hashes;
- richer units metadata;
- provenance graph;
- HDF5 chunking/compression policy;
- parallel HDF5;
- N-to-M restart guarantees.

No proposed improvement is described as implemented until code and evidence exist.
