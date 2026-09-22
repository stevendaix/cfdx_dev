#include "cfdx/core/mesh/mesh.h"
#include "cfdx/io/mesh/mesh_importer.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using cfdx::core::Mesh;
using cfdx::core::Vec3;

void write_file(const std::filesystem::path& path, const std::string& text) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot create " + path.string());
    out << text;
}

std::string point_key(const Vec3& p) {
    std::ostringstream os;
    os.precision(17);
    os << p.x << "," << p.y << "," << p.z;
    return os.str();
}

std::string face_key(const Mesh& mesh, std::size_t face_id) {
    const auto& faces = mesh.faces();
    const auto begin = faces.offsets_data()[face_id];
    const auto end = faces.offsets_data()[face_id + 1];

    std::vector<std::string> points;
    points.reserve(static_cast<std::size_t>(end - begin));
    for (std::uint64_t i = begin; i < end; ++i) {
        points.push_back(point_key(mesh.points().get(faces.vertices_data()[i])));
    }
    std::sort(points.begin(), points.end());
    std::ostringstream os;
    for (const auto& p : points) os << p << ";";
    return os.str();
}

std::multiset<std::string> face_signature(const Mesh& mesh) {
    std::multiset<std::string> result;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f)
        result.insert(face_key(mesh, f));
    return result;
}

bool equivalent(const Mesh& a, const Mesh& b) {
    if (a.n_points() != b.n_points() ||
        a.n_faces() != b.n_faces() ||
        a.n_cells() != b.n_cells())
        return false;

    if (face_signature(a) != face_signature(b))
        return false;

    const auto sa = a.stats();
    const auto sb = b.stats();
    return sa.n_boundary_faces == sb.n_boundary_faces &&
           sa.n_internal_faces == sb.n_internal_faces &&
           a.topo_validate().ok && b.topo_validate().ok;
}

} // namespace

int main() {
    namespace fs = std::filesystem;
    const fs::path root = fs::temp_directory_path() / "cfdx_cross_import_equivalence";
    fs::remove_all(root);
    fs::create_directories(root / "constant" / "polyMesh");

    try {
        // Same unit tetrahedron represented by OpenFOAM native polyMesh and
        // Gmsh/meshio. The test compares semantic topology, not file ordering.
        write_file(root / "constant" / "polyMesh" / "points",
R"(FoamFile { version 2.0; format ascii; class vectorField; object points; }
4
(
(0 0 0)
(1 0 0)
(0 1 0)
(0 0 1)
)
)");

        write_file(root / "constant" / "polyMesh" / "faces",
R"(FoamFile { version 2.0; format ascii; class faceList; object faces; }
4
(
3(0 2 1)
3(0 1 3)
3(1 3 2)
3(2 3 0)
)
)");

        write_file(root / "constant" / "polyMesh" / "owner",
R"(FoamFile { version 2.0; format ascii; class labelList; object owner; }
4
(
0
0
0
0
)
)");

        write_file(root / "constant" / "polyMesh" / "neighbour",
R"(FoamFile { version 2.0; format ascii; class labelList; object neighbour; }
0
(
)
)");

        write_file(root / "constant" / "polyMesh" / "boundary",
R"(FoamFile { version 2.0; format ascii; class polyBoundaryMesh; object boundary; }
1
(
wall
{
    type wall;
    nFaces 4;
    startFace 0;
}
)
)");

        const fs::path gmsh = root / "tetra.msh";
        write_file(gmsh,
R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
1
1 4 2 1 1 1 2 3 4
$EndElements
)");

        Mesh openfoam_mesh;
        Mesh gmsh_mesh;
        if (!cfdx::io::mesh::import_mesh(root.string(), openfoam_mesh)) {
            std::cerr << "OpenFOAM import failed
";
            fs::remove_all(root);
            return 1;
        }
        if (!cfdx::io::mesh::import_mesh(gmsh.string(), gmsh_mesh)) {
            std::cerr << "Gmsh/meshio import failed
";
            fs::remove_all(root);
            return 1;
        }

        if (!equivalent(openfoam_mesh, gmsh_mesh)) {
            std::cerr << "OpenFOAM and Gmsh/meshio imports are not semantically equivalent
";
            fs::remove_all(root);
            return 1;
        }

        fs::remove_all(root);
        std::cout << "cross-import equivalence passed
";
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "
";
        fs::remove_all(root);
        return 1;
    }
}
