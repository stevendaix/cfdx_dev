# 03 — Meshes

**Status: REPOSITORY-GROUNDED. This chapter documents the current mesh model; it does not claim mesh qualification.**

## 1. Why the mesh is part of the numerical method

A finite-volume operator is defined on control volumes and oriented faces:

$$
\int_{V_P}\nabla\cdot\mathbf F\,dV
=\sum_{f\in\partial P}\mathbf F_f\cdot\mathbf S_{f,P}.
$$

Incorrect topology, orientation or geometry changes the discrete operator even when the PDE is unchanged.

## 2. CFDX topology

The current C++ Mesh contains point coordinates, face connectivity, owner/neighbour ownership, cell-to-face connectivity and boundary patches. Face and cell connectivity use CSR-like flat arrays.

For a face, the range is [face_offsets[f], face_offsets[f+1]). For a cell, the range is [cell_offsets[c], cell_offsets[c+1]).

Primary implementation families are core/mesh/mesh.h, face.h, cell.h, point.h, ownership.h, boundary.h and solver_mesh.h.

The Mesh validator checks offset consistency, index ranges, owner/neighbour consistency, patch consistency and cell-face references.

## 3. Face orientation

An internal face has one owner and one neighbour. Its oriented area vector changes sign when viewed from the neighbouring cell:

$$
\mathbf S_{f,N}=-\mathbf S_{f,P}.
$$

Boundary faces use neighbour = -1. This is a conservation invariant.

## 4. Geometry

The geometry layer derives face centres, face area vectors, cell centres and signed/absolute volumes. For a closed 3-D cell the current implementation uses

$$
V=\frac13\sum_f(\mathbf C_f-\mathbf C_P)\cdot\mathbf S_{f,P}.
$$

The signed volume is retained so inverted cells can be rejected.

Relevant files are core/geometry/face_geometry.h, cell_geometry.h, geometry_cache.h, mesh_quality.h and mesh_validator.h.

## 5. Mesh quality

The current quality layer computes skewness and non-orthogonality. For an internal face,

$$
\mathbf d=\mathbf C_N-\mathbf C_P
$$

is compared with the face area vector. Non-orthogonality is an angle; skewness measures transverse displacement of the face centre from the centre line.

The warning thresholds in mesh_validator.h are diagnostics, not universal physical limits.

## 6. Validation chain

$$
\text{input}\rightarrow\text{topology}\rightarrow\text{geometry}\rightarrow
\text{quality}\rightarrow\text{closure}\rightarrow\text{FVM}.
$$

The current validator checks finite coordinates, positive face areas, positive signed cell volumes, surface closure and quality diagnostics.

## 7. Import

Current import families include OpenFOAM and Gmsh native paths plus generic mesh import support. Import is a mathematical data-model transformation: element mapping, orientation, boundary groups and unsupported constructs must be recorded.

Relevant files are io/mesh/mesh_importer.*, io/openfoam/openfoam_importer.*, io/gmsh/gmsh_importer.* and io/vtu/vtu_writer.*.

## 8. Required V&V

The mesh campaign must include topology closure, analytical geometry checks, positive volume, controlled skew/non-orthogonality, stretched cells, polyhedra and valid near-degenerate cells. Invalid meshes must be rejected rather than repaired silently.
