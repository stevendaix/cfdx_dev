#include "openfoam_importer.h"

#include "cfdx/core/mesh/boundary.h"
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <limits>
#include <algorithm>

namespace cfdx::io::openfoam {
namespace {

std::string strip_comments(std::string text) {
    text = std::regex_replace(text, std::regex(R"(//[^\n]*)"), "");
    text = std::regex_replace(text, std::regex(R"(/\*[\s\S]*?\*/)"), "");
    return text;
}

bool read_text(const std::filesystem::path& path, std::string& text) {
    std::ifstream in(path);
    if (!in) return false;
    std::ostringstream buffer;
    buffer << in.rdbuf();
    text = strip_comments(buffer.str());
    return true;
}

bool read_points(const std::filesystem::path& path, std::vector<cfdx::core::Vec3>& points) {
    std::string text;
    if (!read_text(path, text)) return false;
    const auto begin = text.find('(');
    const auto end = text.rfind(')');
    if (begin == std::string::npos || end == std::string::npos || end <= begin) return false;

    const std::string body = text.substr(begin + 1, end - begin - 1);
    std::regex point_re(R"(\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\))");
    for (std::sregex_iterator it(body.begin(), body.end(), point_re), e; it != e; ++it) {
        points.emplace_back(
            std::stod((*it)[1].str()),
            std::stod((*it)[2].str()),
            std::stod((*it)[3].str()));
    }
    return !points.empty();
}

bool read_label_list(const std::filesystem::path& path, std::vector<std::int64_t>& values) {
    std::string text;
    if (!read_text(path, text)) return false;
    const auto begin = text.find('(');
    const auto end = text.rfind(')');
    if (begin == std::string::npos || end == std::string::npos || end <= begin) return false;
    std::istringstream in(text.substr(begin + 1, end - begin - 1));
    std::int64_t value = 0;
    while (in >> value) values.push_back(value);
    return true;
}

bool read_faces(const std::filesystem::path& path,
                std::vector<std::vector<std::uint32_t>>& faces) {
    std::string text;
    if (!read_text(path, text)) return false;
    const auto begin = text.find('(');
    const auto end = text.rfind(')');
    if (begin == std::string::npos || end == std::string::npos || end <= begin) return false;
    const std::string body = text.substr(begin + 1, end - begin - 1);

    std::regex face_re(R"(\b(\d+)\s*\(([^()]*)\))");
    for (std::sregex_iterator it(body.begin(), body.end(), face_re), e; it != e; ++it) {
        const std::size_t count = static_cast<std::size_t>(std::stoull((*it)[1].str()));
        std::istringstream in((*it)[2].str());
        std::vector<std::uint32_t> face;
        std::uint64_t index = 0;
        while (in >> index) {
            if (index > std::numeric_limits<std::uint32_t>::max()) return false;
            face.push_back(static_cast<std::uint32_t>(index));
        }
        if (face.size() != count || face.size() < 3) return false;
        faces.push_back(std::move(face));
    }
    return !faces.empty();
}

cfdx::core::PatchType patch_type(const std::string& type) {
    if (type == "symmetryPlane" || type == "symmetry") return cfdx::core::PatchType::SYMMETRY;
    if (type == "cyclic" || type == "cyclicAMI" || type == "cyclicACMI") return cfdx::core::PatchType::PERIODIC;
    if (type == "processor") return cfdx::core::PatchType::INTERFACE;
    if (type == "empty") return cfdx::core::PatchType::EMPTY;
    if (type == "wall") return cfdx::core::PatchType::WALL;
    if (type == "patch") return cfdx::core::PatchType::UNKNOWN;
    return cfdx::core::patch_type_from_string(type);
}

bool read_boundary(const std::filesystem::path& path,
                   cfdx::core::BoundaryPatches& boundary) {
    std::string text;
    if (!read_text(path, text)) return false;

    std::regex patch_re(
        R"(([A-Za-z_][A-Za-z0-9_.-]*)\s*\{([\s\S]*?)\})");
    std::regex type_re(R"(\btype\s+([^;]+);)");
    std::regex nfaces_re(R"(\bnFaces\s+(\d+)\s*;)");
    std::regex start_re(R"(\bstartFace\s+(\d+)\s*;)");

    bool found = false;
    for (std::sregex_iterator it(text.begin(), text.end(), patch_re), e; it != e; ++it) {
        const std::string name = (*it)[1].str();
        const std::string body = (*it)[2].str();
        std::smatch type_match, count_match, start_match;
        if (!std::regex_search(body, type_match, type_re)) continue;
        if (!std::regex_search(body, count_match, nfaces_re)) continue;
        if (!std::regex_search(body, start_match, start_re)) continue;
        const std::size_t nfaces = static_cast<std::size_t>(std::stoull(count_match[1].str()));
        const std::size_t start = static_cast<std::size_t>(std::stoull(start_match[1].str()));

        cfdx::core::Patch patch;
        patch.name = name;
        patch.type = patch_type(type_match[1].str());
        patch.face_ids.reserve(nfaces);
        for (std::size_t i = 0; i < nfaces; ++i) {
            patch.face_ids.push_back(static_cast<cfdx::core::FaceIndex>(start + i));
        }
        boundary.add_patch(patch);
        found = true;
    }
    return found;
}

} // namespace

bool import_openfoam_case(const std::string& case_path, cfdx::core::Mesh& mesh) {
    namespace fs = std::filesystem;
    const fs::path root(case_path);
    const fs::path poly = root / "constant" / "polyMesh";
    if (!fs::is_directory(poly)) return false;

    std::vector<cfdx::core::Vec3> points;
    std::vector<std::vector<std::uint32_t>> faces;
    std::vector<std::int64_t> owner;
    std::vector<std::int64_t> neighbour;
    cfdx::core::BoundaryPatches boundary;

    if (!read_points(poly / "points", points) ||
        !read_faces(poly / "faces", faces) ||
        !read_label_list(poly / "owner", owner) ||
        !read_label_list(poly / "neighbour", neighbour)) {
        return false;
    }
    if (owner.size() != faces.size() || neighbour.size() > owner.size()) return false;

    mesh.clear();
    mesh.points().resize(points.size());
    for (std::size_t i = 0; i < points.size(); ++i)
        mesh.points().set(i, points[i].x, points[i].y, points[i].z);

    mesh.faces().build_from_scratch(
        std::vector<std::vector<cfdx::core::FaceIndex>>(faces.begin(), faces.end()));

    mesh.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        if (owner[f] < 0) return false;
        mesh.ownership().set_owner(f, static_cast<cfdx::core::CellIndex>(owner[f]));
        mesh.ownership().set_neighbour(
            f, f < neighbour.size() ? neighbour[f] : cfdx::core::FaceOwnership::BOUNDARY);
    }

    // Reconstruct cell -> face CSR directly from owner/neighbour.
    std::size_t n_cells = 0;
    for (std::int64_t value : owner) {
        if (value < 0) return false;
        n_cells = std::max(n_cells, static_cast<std::size_t>(value) + 1);
    }
    for (std::int64_t value : neighbour) {
        if (value >= 0)
            n_cells = std::max(n_cells, static_cast<std::size_t>(value) + 1);
    }

    std::vector<std::vector<cfdx::core::FaceIndex>> cell_faces(n_cells);
    for (std::size_t f = 0; f < faces.size(); ++f) {
        cell_faces[static_cast<std::size_t>(owner[f])].push_back(
            static_cast<cfdx::core::FaceIndex>(f));
        if (f < neighbour.size() && neighbour[f] >= 0)
            cell_faces[static_cast<std::size_t>(neighbour[f])].push_back(
                static_cast<cfdx::core::FaceIndex>(f));
    }
    for (const auto& cf : cell_faces) mesh.cells().push_cell(cf);

    const fs::path boundary_file = poly / "boundary";
    if (fs::exists(boundary_file)) {
        if (!read_boundary(boundary_file, boundary)) return false;
        mesh.set_boundary(boundary);
    }

    const auto validation = mesh.topo_validate();
    return validation.ok;
}

bool read_openfoam_mesh_meshio(const std::string& case_path,
                               std::vector<cfdx::core::Vec3>& points,
                               std::vector<OpenFOAMFace>& faces,
                               std::vector<PatchDef>& patches) {
    cfdx::core::Mesh mesh;
    if (!import_openfoam_case(case_path, mesh)) return false;

    points.resize(mesh.n_points());
    for (std::size_t i = 0; i < mesh.n_points(); ++i)
        points[i] = {mesh.points().x(i), mesh.points().y(i), mesh.points().z(i)};

    faces.resize(mesh.n_faces());
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto off = mesh.faces().offsets_data()[f];
        const auto end = mesh.faces().offsets_data()[f + 1];
        faces[f].n = static_cast<int>(end - off);
        faces[f].indices.assign(
            mesh.faces().vertices_data() + off,
            mesh.faces().vertices_data() + end);
    }
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        const auto& patch = mesh.boundary().patch(p);
        patches.emplace_back(
            patch.name, static_cast<int>(patch.type), static_cast<int>(patch.size()),
            std::vector<std::uint32_t>(patch.face_ids.begin(), patch.face_ids.end()));
    }
    return true;
}

} // namespace cfdx::io::openfoam
