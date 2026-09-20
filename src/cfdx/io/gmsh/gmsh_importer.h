#pragma once

#include <string>
#include <vector>
#include <meshio/meshio.h>

namespace cfdx::io::gmsh {

//! M0.10-T02 — Gmsh importer
//! Import Gmsh .msh files using meshio adapter per CFDX spec §23.
//!
//! Supported:
//!   - Mesh topology (points, lines, triangles, quads, tets, hexes, prisms)
//!   - Physical groups (boundary patches mapped to CFDX patches)
//!   - Cell/point fields (scalars, vectors)

//! Import Gmsh mesh file into CFDX Mesh structure
bool import_gmsh_mesh(const std::string& mshPath, Mesh& mesh);
//! Import Gmsh field file into CFDX scalar field
bool import_gmsh_scalar(const std::string& mshPath, const std::string& fieldName, ScalarCellField& field);
//! Import Gmsh vector field into CFDX vector field
bool import_gmsh_vector(const std::string& mshPath, const std::string& fieldName, VectorCellField& field);

} // namespace cfdx::io::gmsh
