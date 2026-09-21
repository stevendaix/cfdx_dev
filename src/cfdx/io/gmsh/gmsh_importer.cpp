#include "gmsh_importer.h"
#include "cfdx/io/mesh/mesh_importer.h"

namespace cfdx::io::gmsh {

bool import_gmsh_mesh(const std::string& path, cfdx::core::Mesh& mesh) {
    // meshio handles Gmsh MSH 2.2/4.0/4.1, ASCII and binary, including
    // sparse node tags and physical groups.
    return cfdx::io::mesh::import_mesh(path, mesh);
}

bool import_gmsh_scalar(const std::string&, const std::string&,
                       cfdx::core::ScalarCellField&) {
    return false;
}

} // namespace cfdx::io::gmsh
