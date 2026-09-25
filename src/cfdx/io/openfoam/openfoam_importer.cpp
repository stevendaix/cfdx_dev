#include "openfoam_importer.h"

#include "cfdx/core/mesh/boundary.h"
#include <fstream>
#include <sstream>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <regex>

namespace cfdx::io::openfoam {
namespace {

std::string strip_comments(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    bool line=false, block=false;
    for(std::size_t i=0;i<text.size();++i) {
        const char c=text[i];
        const char n=(i+1<text.size()?text[i+1]:'\0');
        if(line) { if(c=='\n') { line=false; out.push_back(c); } continue; }
        if(block) { if(c=='*'&&n=='/') { block=false; ++i; } else if(c=='\n') out.push_back('\n'); continue; }
        if(c=='/'&&n=='/') { line=true; ++i; continue; }
        if(c=='/'&&n=='*') { block=true; ++i; continue; }
        out.push_back(c);
    }
    return out;
}
bool read_text(const std::filesystem::path& path, std::string& text) {
    std::ifstream in(path);
    if(!in) return false;
    std::ostringstream buffer; buffer << in.rdbuf();
    text=strip_comments(buffer.str());
    return true;
}
struct Cursor {
    const char* p; const char* end;
    void ws() { while(p<end && std::isspace(static_cast<unsigned char>(*p))) ++p; }
    bool expect(char c) { ws(); if(p<end && *p==c){++p;return true;} return false; }
    bool size(std::size_t& v) {
        ws(); if(p>=end || *p<'0' || *p>'9') return false;
        std::uint64_t x=0;
        while(p<end && std::isdigit(static_cast<unsigned char>(*p))) {
            x=x*10+static_cast<unsigned>(*p-'0'); ++p;
            if(x>std::numeric_limits<std::size_t>::max()) return false;
        }
        v=static_cast<std::size_t>(x); return true;
    }
    bool i64(std::int64_t& v) {
        ws(); char* q=nullptr; errno=0;
        const double probe=0.0; (void)probe;
        long long x=std::strtoll(p,&q,10);
        if(q==p || errno==ERANGE || q>end) return false;
        p=q; v=static_cast<std::int64_t>(x); return true;
    }
    bool real(double& v) {
        ws(); char* q=nullptr; errno=0;
        double x=std::strtod(p,&q);
        if(q==p || errno==ERANGE || q>end || !std::isfinite(x)) return false;
        p=q; v=x; return true;
    }
};
bool read_declared_count(const std::string& text, std::size_t& count) {
    const auto pos=text.find('\n');
    const char* b=text.data(); const char* e=b+text.size();
    Cursor c{b,e};
    while(c.p<c.end) {
        c.ws(); std::size_t v=0;
        if(c.size(v)) { c.ws(); if(c.p<c.end && *c.p=='(') { count=v; return true; } }
        while(c.p<c.end && *c.p!='\n') ++c.p;
    }
    (void)pos; return false;
}
bool read_points(const std::filesystem::path& path, std::vector<cfdx::core::Vec3>& points) {
    std::string text; if(!read_text(path,text)) return false;
    std::size_t declared=0; if(!read_declared_count(text,declared) || declared==0) return false;
    const auto pos=text.find('(', text.find(std::to_string(declared)));
    if(pos==std::string::npos) return false;
    Cursor c{text.data()+pos,text.data()+text.size()};
    if(!c.expect('(')) return false;
    points.clear(); points.reserve(declared);
    for(std::size_t i=0;i<declared;++i) {
        if(!c.expect('(')) return false;
        double x,y,z; if(!c.real(x)||!c.real(y)||!c.real(z)||!c.expect(')')) return false;
        points.emplace_back(x,y,z);
    }
    c.ws(); return c.expect(')') && c.p==c.end;
}
bool read_label_list(const std::filesystem::path& path, std::vector<std::int64_t>& values) {
    std::string text; if(!read_text(path,text)) return false;
    std::size_t declared=0; if(!read_declared_count(text,declared)) return false;
    const auto pos=text.find('(', text.find(std::to_string(declared)));
    if(pos==std::string::npos) return false;
    Cursor c{text.data()+pos,text.data()+text.size()}; if(!c.expect('(')) return false;
    values.clear(); values.reserve(declared);
    for(std::size_t i=0;i<declared;++i) { std::int64_t v; if(!c.i64(v)) return false; values.push_back(v); }
    c.ws(); return c.expect(')') && c.p==c.end && values.size()==declared;
}
bool read_faces(const std::filesystem::path& path,std::vector<std::vector<std::uint32_t>>& faces) {
    std::string text; if(!read_text(path,text)) return false;
    std::size_t declared=0; if(!read_declared_count(text,declared)||declared==0) return false;
    const auto pos=text.find('(', text.find(std::to_string(declared)));
    if(pos==std::string::npos) return false;
    Cursor c{text.data()+pos,text.data()+text.size()}; if(!c.expect('(')) return false;
    faces.clear(); faces.reserve(declared);
    for(std::size_t i=0;i<declared;++i) {
        std::size_t n=0; if(!c.size(n)||n<3||!c.expect('(')) return false;
        std::vector<std::uint32_t> face; face.reserve(n);
        for(std::size_t j=0;j<n;++j) { std::size_t v; if(!c.size(v)||v>std::numeric_limits<std::uint32_t>::max()) return false; face.push_back(static_cast<std::uint32_t>(v)); }
        if(!c.expect(')')) return false;
        faces.push_back(std::move(face));
    }
    c.ws(); return c.expect(')') && c.p==c.end && faces.size()==declared;
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
        for (std::size_t i = 0; i < nfaces; ++i)
            patch.face_ids.push_back(static_cast<cfdx::core::FaceIndex>(start + i));
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

    for (const auto& face : faces)
        for (const auto vertex : face)
            if (vertex >= points.size()) return false;

    std::vector<std::vector<cfdx::core::FaceIndex>> face_vertices;
    face_vertices.reserve(faces.size());
    for (const auto& face : faces) {
        face_vertices.emplace_back();
        face_vertices.back().reserve(face.size());
        for (const std::uint32_t vertex : face) {
            face_vertices.back().push_back(
                static_cast<cfdx::core::FaceIndex>(vertex));
        }
    }
    mesh.faces().build_from_scratch(std::move(face_vertices));

    mesh.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        if (owner[f] < 0) return false;
        if (f < neighbour.size() && neighbour[f] < -1) return false;
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
    if (!fs::exists(boundary_file) || !read_boundary(boundary_file, boundary))
        return false;
    for (std::size_t p = 0; p < boundary.n_patches(); ++p) {
        const auto& patch = boundary.patch(p);
        for (const auto face_id : patch.face_ids)
            if (face_id >= faces.size()) return false;
    }
    mesh.set_boundary(boundary);

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
