#pragma once

#include "cfdx/core/mesh/mesh.h"
#include <string>

namespace cfdx::io::mesh {

enum class MeshFormat {
    OPENFOAM,
    MESHIO,
    UNKNOWN
};

MeshFormat detect_format(const std::string& path);

// Universal mesh entry point.
// OpenFOAM cases use the native polyMesh reader; all other files are sent
// through meshio and the CFDX-HDF5 interchange reader.
bool import_mesh(const std::string& path, cfdx::core::Mesh& mesh);

} // namespace cfdx::io::mesh
