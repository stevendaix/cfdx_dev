# CFDX mesh import architecture

## Design

CFDX now has one public entry point:

```cpp
#include "cfdx/io/mesh/mesh_importer.h"

cfdx::core::Mesh mesh;
if (!cfdx::io::mesh::import_mesh(path, mesh)) {
    // import failed
}
```

The dispatcher uses two backends:

- **OpenFOAM case**: native reader for `constant/polyMesh`.
- **File formats**: meshio Python bridge, followed by the CFDX HDF5 interchange reader.

This keeps the CFDX topology representation independent of individual file
formats.

## OpenFOAM

The native reader consumes:

- `points`
- `faces`
- `owner`
- `neighbour`
- `boundary`

It reconstructs cell-to-face CSR from owner/neighbour and preserves patch
names, patch types and patch face IDs. OpenFOAM stores faces first as internal
faces and then boundary faces; CFDX does not rely on this ordering and uses the
explicit owner/neighbour information.

The reader also accepts an empty `neighbour` list for a mesh with no internal
faces.

## meshio formats

The bridge uses `meshio.read()`, so supported input formats follow the
installed meshio version rather than a hard-coded CFDX extension list. This
covers, among others, Gmsh MSH, ANSYS MSH, CGNS, Exodus, MED/Salome, Nastran,
Netgen, STL, SU2, VTK/VTU, XDMF and UGRID.

Gmsh MSH 2.2/4.0/4.1, ASCII and binary are therefore handled by the same
backend, including sparse node tags and Gmsh physical-group metadata.

## Topology normalization

The bridge converts supported meshio cell blocks into CFDX face/cell CSR.

Supported volume topologies include:

- tetrahedral
- hexahedral
- wedge/prism
- pyramid
- voxel
- polyhedron cell blocks when supplied by meshio

Common second-order variants are reduced to their corner-node topology:
tetra10, hexahedron20/27, wedge15/18, pyramid13/14, triangle6 and quad8/9.

For surface-only meshes, CFDX builds edge faces and creates a deterministic
`boundary` patch when the source format does not provide physical groups.

Unsupported cell types are reported explicitly rather than silently converted.

## Validation

Every imported mesh is checked with `Mesh::topo_validate()` before being
accepted by the universal C++ importer.

The CI contains regression tests for:

- native OpenFOAM polyMesh import;
- meshio -> CFDX HDF5 conversion;
- CMake compilation;
- CTest;
- Release and AddressSanitizer/UndefinedBehaviorSanitizer builds.

## Important limitation

"All mesh formats supported by meshio" does not mean that every possible
element topology can be represented without loss in a finite-volume solver.
Format support and topology support are deliberately separated. A format can
be read by meshio while an unsupported cell topology is reported and rejected
or skipped instead of producing an invalid CFDX mesh.
