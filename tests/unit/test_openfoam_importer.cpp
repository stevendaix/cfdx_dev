#include "cfdx/io/openfoam/openfoam_importer.h"
#include "cfdx/core/mesh/mesh.h"

#include <filesystem>
#include <fstream>
#include <iostream>

int main() {
    namespace fs = std::filesystem;
    const fs::path root = fs::temp_directory_path() / "cfdx_openfoam_reader_test";
    const fs::path poly = root / "constant" / "polyMesh";
    fs::create_directories(poly);

    auto write = [](const fs::path& p, const char* text) {
        std::ofstream out(p);
        if (!out) return false;
        out << text;
        return true;
    };

    const bool ok =
        write(poly / "points", R"(FoamFile { version 2.0; format ascii; class vectorField; object points; }
4
(
(0 0 0)
(1 0 0)
(0 1 0)
(0 0 1)
)
)") &&
        write(poly / "faces", R"(FoamFile { version 2.0; format ascii; class faceList; object faces; }
4
(
3(0 2 1)
3(0 1 3)
3(1 2 3)
3(2 0 3)
)
)") &&
        write(poly / "owner", R"(FoamFile { version 2.0; format ascii; class labelList; object owner; }
4
(
0
0
0
0
)
)") &&
        write(poly / "neighbour", R"(FoamFile { version 2.0; format ascii; class labelList; object neighbour; }
0
(
)
)") &&
        write(poly / "boundary", R"(FoamFile { version 2.0; format ascii; class polyBoundaryMesh; object boundary; }
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

    if (!ok) return 1;

    cfdx::core::Mesh mesh;
    const bool imported = cfdx::io::openfoam::import_openfoam_case(root.string(), mesh);
    const auto stats = mesh.stats();
    const bool valid = imported && mesh.topo_validate().ok &&
                       stats.n_points == 4 && stats.n_faces == 4 &&
                       stats.n_cells == 1 && stats.n_boundary_faces == 4 &&
                       mesh.boundary().n_patches() == 1;

    fs::remove_all(root);
    if (!valid) {
        std::cerr << "OpenFOAM reader test failed\n";
        return 1;
    }
    return 0;
}
