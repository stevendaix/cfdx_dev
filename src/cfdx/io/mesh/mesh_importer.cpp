#include "mesh_importer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/openfoam/openfoam_importer.h"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>

namespace cfdx::io::mesh {
namespace {
std::string quote(const std::string& value) {
    std::string out = "'";
    for (char c : value) {
        if (c == static_cast<char>(39)) out += "'\\''";
        else out += c;
    }
    return out + "'";
}

bool is_openfoam_case(const std::filesystem::path& path) {
    const auto poly = path / "constant" / "polyMesh";
    return std::filesystem::is_directory(poly) &&
           std::filesystem::exists(poly / "points") &&
           std::filesystem::exists(poly / "faces") &&
           std::filesystem::exists(poly / "owner") &&
           std::filesystem::exists(poly / "neighbour");
}

bool import_with_meshio(const std::filesystem::path& input, cfdx::core::Mesh& mesh) {
#ifndef CFDX_SOURCE_DIR
    (void)input;
    (void)mesh;
    return false;
#else
    namespace fs = std::filesystem;
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const fs::path output = fs::temp_directory_path() /
        ("cfdx_meshio_" + std::to_string(stamp) + ".h5");
    const fs::path script = fs::path(CFDX_SOURCE_DIR) / "scripts" / "meshio_import.py";

    std::ostringstream command;
    command <<
#ifdef CFDX_PYTHON_EXECUTABLE
        quote(CFDX_PYTHON_EXECUTABLE)
#else
        "python3"
#endif
            << " " << quote(script.string()) << " "
            << quote(input.string()) << " " << quote(output.string());

    const int rc = std::system(command.str().c_str());
    if (rc != 0 || !fs::exists(output)) {
        std::error_code ec;
        fs::remove(output, ec);
        return false;
    }

    const bool ok = cfdx::io::read_mesh_hdf5(output.string(), mesh);
    std::error_code ec;
    fs::remove(output, ec);
    return ok && mesh.topo_validate().ok;
#endif
}
} // namespace

MeshFormat detect_format(const std::string& path) {
    if (is_openfoam_case(std::filesystem::path(path))) return MeshFormat::OPENFOAM;
    if (!std::filesystem::path(path).extension().empty()) return MeshFormat::MESHIO;
    return MeshFormat::UNKNOWN;
}

bool import_mesh(const std::string& path, cfdx::core::Mesh& mesh) {
    const std::filesystem::path input(path);
    if (is_openfoam_case(input))
        return cfdx::io::openfoam::import_openfoam_case(path, mesh);
    if (!std::filesystem::is_regular_file(input)) return false;
    return import_with_meshio(input, mesh);
}
} // namespace cfdx::io::mesh
