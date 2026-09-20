#pragma once

#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/face_geometry.h"
#include <string>

// Local type definitions for OpenFOAM import skeleton
namespace cfdx::io::openfoam {
    using real3 = double[3];
    struct Face {
        int n;
        uint32_t* indices;
    };
    struct PatchDef {
        std::string name;
        int type;
        std::vector<uint32_t> face_ids;
    };
}

// Meshio adapter: optional; if meshio/meshio.h unavailable, basic parser used instead
#ifndef HAS_MESHIO
#define HAS_MESHIO 0
#endif

namespace cfdx::io::openfoam {

// M0.10-T01 — OpenFOAM importer
// Import OpenFOAM meshes and case data using meshio adapter per spec §23.
//
// Supported:
//   - Polyhedral mesh topology (points, faces, cells, owner/neighbour)
//   - Boundary patches with names and types
//   - Field data conversion (scalar, vector, tensor)
//   - Geometry export (pts/faces for CFDX native format)
//
// Not yet supported (future):
//   - Compressed/serialized case files
//   - Dynamic mesh motion
//   - AMR adaptive refinement

//! Import OpenFOAM case directory into CFDX Mesh
/*!
  \param ofCasePath Path to OpenFOAM case directory (contains system/constant)
  \param mesh Output CFDX Mesh containing imported topology and patches
  \return true on success, false on failure (invalid case format, unsupported features)
  \note Uses meshio as the underlying adapter for meshio-based format parsing.
  \note Polyhedral cell support: stores as unstructured cells with face connectivity.
  \note Boundary patches: union of all patches in constant/polyMesh/ directory.
*/
bool import_openfoam_case(const std::string& ofCasePath, Mesh& mesh);

//! Convert OpenFOAM field to CFDX field
/*!
  \param ofFieldName Name of the OpenFOAM field (e.g. "U", "p", "k", "epsilon")
  \param ofCasePath Path to OpenFOAM case directory
  \param cfdxField Output CFDX field (will be resized/reallocated)
  \return true on success
  \note Handles: volScalarField, volVectorField, volTensorField
  \note Supported types: scalar, vector, symmetric tensor (6 components)
  \note Unit metadata preserved when available in field metadata.
*/
bool import_openfoam_field(const std::string& ofFieldName,
                           const std::string& ofCasePath,
                           ScalarCellField& cfdxField);

//! Convert OpenFOAM field to CFDX field (vector version)
/*!
  \param ofFieldName Name of the OpenFOAM field (e.g. "U")
  \param ofCasePath Path to OpenFOAM case directory
  \param cfdxField Output CFDX vector field
  \return true on success
*/
bool import_openfoam_field(const std::string& ofFieldName,
                           const std::string& ofCasePath,
                           VectorCellField& cfdxField);

//! Read OpenFOAM polyhedral mesh using meshio adapter
/*!
  \param ofCasePath Path to OpenFOAM case directory
  \param points Output vector of point positions (3*n)
  \param faces Output face indices (n faces, each with face point count + indices)
  \param patches Output boundary patch definitions (name, type, face indices)
  \return true on success
  \note Uses meshio::read_mesh as the primary adapter.
  \note Polyhedral faces: triangulated on import (each poly face → 3 triangle faces).
  \note Conservation: face areas/Sf computed from vertex positions.
*/
bool read_openfoam_mesh_meshio(const std::string& ofCasePath,
                               std::vector<real3>& points,
                               std::vector<Face>& faces,
                               std::vector<PatchDef>& patches);

} // namespace cfdx::io::openfoam