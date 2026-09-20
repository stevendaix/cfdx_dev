#include "gmsh_importer.h"
#include <meshio/meshio.h>
#include <fstream>

using namespace cfdx::io::gmsh;

bool import_gmsh_mesh(const std::string&, Mesh&) {
    // Skeleton — meshio adapter for Gmsh .msh format
    return false; // Not fully implemented
}
bool import_gmsh_scalar(const std::string&, const std::string&, ScalarCellField&) {
    return false;
}
bool import_gmsh_vector(const std::string&, const std::string&, VectorCellField&) {
    return false;
}
