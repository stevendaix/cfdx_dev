#pragma once

#include <string>
#include <vector>
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"

namespace cfdx::io::gmsh {

//! M0.10-T02 — Gmsh importer
//! Import Gmsh .msh files using meshio adapter per CFDX spec §23.
//!
//! Supported:
//!   - Mesh topology (points, lines, triangles, quads, tets, hexes, prisms)
//!   - Physical groups (boundary patches mapped to CFDX patches)
//!   - Cell/point fields (scalars, vectors)

//! Import Gmsh mesh file into CFDX Mesh structure
bool import_gmsh_mesh(const std::string& mshPath, cfdx::core::Mesh& mesh);
//! Import Gmsh field file into CFDX scalar field
bool import_gmsh_scalar(const std::string& mshPath, const std::string& fieldName, cfdx::core::ScalarCellField& field);
//! Import Gmsh vector field into CFDX vector field (stub)
// bool import_gmsh_vector(const std::string& mshPath, const std::string& fieldName, cfdx::core::VectorCellField& field);

} // namespace cfdx::io::gmsh
