#include "cfdx/core/field/field.h"
#include "gmsh_importer.h"
#include "cfdx/core/mesh/mesh.h"
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace cfdx {
namespace io {
namespace gmsh {

using namespace cfdx::core;

// Helper: read all nodes from $Nodes ... $EndNodes section
static bool parseNodes(const std::string& content,
                        std::map<uint32_t, cfdx::core::Vec3>& points_map) {
    std::istringstream stream(content);
    std::string line;
    bool in_nodes = false;
    while (std::getline(stream, line)) {
        std::istringstream ls(line);
        std::string token;
        if (line.find("$Nodes") != std::string::npos) { in_nodes = true; continue; }
        if (line.find("$EndNodes") != std::string::npos) { in_nodes = false; break; }
        if (!in_nodes) continue;
        std::string first_token;
        uint32_t id;
        double x, y, z;
        ls >> id >> x >> y >> z;
        points_map[id] = cfdx::core::Vec3(x, y, z);
    }
    return !points_map.empty();
}

// Helper: read elements and build face topology from $Elements ... $EndElements
// Simplified: triangulate elements, build owner/neighbour for each face
static bool parseElements(const std::string& content,
                            Mesh& mesh,
                            const std::map<uint32_t, cfdx::core::Vec3>& points_map) {
    std::istringstream stream(content);
    std::string line;
    bool in_elements = false;
    std::map<uint32_t, std::vector<uint32_t>> cell_faces;
    // For simplicity: build a basic cell with points from first triangle
    // Real implementation would parse all element types and build full CSR
    while (std::getline(stream, line)) {
        std::istringstream ls(line);
        std::string token;
        if (line.find("$Elements") != std::string::npos) { in_elements = true; continue; }
        if (line.find("$EndElements") != std::string::npos) { in_elements = false; break; }
        if (!in_elements) continue;
        uint32_t id, type, tag;
        ls >> id >> type >> tag;
        // Read node IDs based on element type
        // Simplified: assume triangles (type 2 in Gmsh)
        std::vector<uint32_t> node_ids;
        uint32_t nid;
        while (ls >> nid) node_ids.push_back(nid);
        if (node_ids.size() >= 3 && id > 0) {
            mesh.cells().push_cell({node_ids.begin(), node_ids.end()});
        }
    }
    return true;
}

// Helper: parse $PhysicalNames section
static std::map<uint32_t, std::string> parsePhysicalNames(const std::string& content) {
    std::map<uint32_t, std::string> physical_names;
    std::istringstream stream(content);
    std::string line;
    bool in_section = false;
    while (std::getline(stream, line)) {
        if (line.find("$PhysicalNames") != std::string::npos) { in_section = true; continue; }
        if (line.find("$EndPhysicalNames") != std::string::npos) { in_section = false; break; }
        if (!in_section) continue;
        std::istringstream ls(line);
        uint32_t dim, tag;
        std::string name_quoted;
        ls >> dim >> tag;
        std::getline(ls, name_quoted);  // rest of line including quotes
        // Trim quotes and whitespace
        size_t start = name_quoted.find('"');
        size_t end = name_quoted.rfind('"');
        std::string name = (start != std::string::npos && end != std::string::npos && end > start)
            ? name_quoted.substr(start + 1, end - start - 1) : "";
        physical_names[tag] = name;
    }
    return physical_names;
}

// Helper: parse $PhysicalGroups section
static std::map<uint32_t, std::vector<uint32_t>> parsePhysicalGroups(const std::string& content) {
    std::map<uint32_t, std::vector<uint32_t>> groups;
    std::istringstream stream(content);
    std::string line;
    bool in_section = false;
    while (std::getline(stream, line)) {
        if (line.find("$PhysicalGroups") != std::string::npos) { in_section = true; continue; }
        if (line.find("$EndPhysicalGroups") != std::string::npos) { in_section = false; break; }
        if (!in_section) continue;
        std::istringstream ls(line);
        uint32_t dim, tag, num_entities;
        ls >> dim >> tag >> num_entities;
        uint32_t entity;
        std::vector<uint32_t> entities;
        while (ls >> entity) {
            entities.push_back(entity);
        }
        groups[tag] = entities;
    }
    return groups;
}

bool import_gmsh_mesh(const std::string& filename, Mesh& mesh) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << "\n";
        return false;
    }

    // ============================================================
    // 1. LECTURE DES POINTS
    // ============================================================
    std::string line;
    std::vector<std::array<double, 3>> points_raw;
    bool in_nodes = false;
    std::size_t expected_nodes = 0;
    std::size_t node_read = 0;

    while (std::getline(file, line)) {
        if (line == "$Nodes") {
            std::getline(file, line);  // nombre de noeuds
            expected_nodes = std::stoul(line);
            points_raw.resize(expected_nodes);
            in_nodes = true;
            continue;
        }
        if (line == "$EndNodes") {
            in_nodes = false;
            break;
        }
        if (!in_nodes) continue;
        std::istringstream iss(line);
        int id;
        iss >> id;
        double x, y, z;
        iss >> x >> y >> z;
        // Gmsh utilise des indices 1-based, CFDX 0-based
        if (id > 0 && static_cast<std::size_t>(id) <= expected_nodes) {
            points_raw[id - 1][0] = x;
            points_raw[id - 1][1] = y;
            points_raw[id - 1][2] = z;
        }
        node_read++;
    }

    // Remplir mesh.points()
    cfdx::core::PointCloud points(points_raw.size());
    for (std::size_t i = 0; i < points_raw.size(); ++i) {
        points.set(i, points_raw[i][0], points_raw[i][1], points_raw[i][2]);
    }
    mesh.points() = points;

    std::cout << "[GMSH_DEBUG] Points lus: " << mesh.n_points() << " (attendu=" << expected_nodes << ")\n";

    // ============================================================
    // 2. LECTURE DES ÉLÉMENTS (cellules 2D uniquement : triangles type 2, quads type 3)
    // ============================================================
    file.clear();
    file.seekg(0);
    std::vector<std::vector<cfdx::core::PointIndex>> raw_cells;
    bool in_elements = false;
    std::size_t expected_elements = 0;
    std::size_t element_read = 0;

    while (std::getline(file, line)) {
        if (line == "$Elements") {
            std::getline(file, line);  // nombre d'éléments
            expected_elements = std::stoul(line);
            in_elements = true;
            continue;
        }
        if (line == "$EndElements") {
            in_elements = false;
            break;
        }
        if (!in_elements) continue;

        std::istringstream iss(line);
        int id, type, tag;
        iss >> id >> type;
        // Ignorer le tag physique (lecture simplifiée)
        int phys_tag;
        if (!(iss >> phys_tag)) phys_tag = 0;

        std::vector<cfdx::core::PointIndex> nodes;
        int nid;
        // Lire les IDs des noeuds en fonction du type
        // Type 2 = triangle (3 noeuds), Type 3 = quad (4 noeuds)
        int expected_nodes_in_element = (type == 2) ? 3 : (type == 3) ? 4 : 0;
        if (expected_nodes_in_element > 0) {
            for (int k = 0; k < expected_nodes_in_element; ++k) {
                if (iss >> nid) {
                    // Conversion 1-based (Gmsh) -> 0-based (CFDX)
                    nodes.push_back(static_cast<cfdx::core::PointIndex>(nid - 1));
                }
            }
            if (!nodes.empty() && nodes.size() == static_cast<std::size_t>(expected_nodes_in_element)) {
                raw_cells.push_back(nodes);
            }
        }
        element_read++;
    }

    std::cout << "[GMSH_DEBUG] Éléments 2D lus: " << raw_cells.size()
              << " (attendu=" << expected_elements << ")\n";

    // ============================================================
    // 3. CONSTRUCTION DE LA TOPOLOGIE FVM
    // ============================================================
    std::cout << "[GMSH_DEBUG] Construction topologie FVM...\n";

    // ============================================================
    // 3. CONSTRUCTION DE LA TOPOLOGIE FVM (version accumulation)
    // ============================================================
    std::cout << "[GMSH_DEBUG] Construction topologie FVM (accumulation)...\n";

    std::map<std::vector<cfdx::core::PointIndex>, std::pair<std::size_t, std::size_t>> face_map;
    std::vector<std::vector<cfdx::core::PointIndex>> all_faces_acc;

    mesh.faces().clear();
    mesh.ownership().clear();
    mesh.cells().clear();

    for (std::size_t c = 0; c < raw_cells.size(); ++c) {
        const auto& nodes = raw_cells[c];
        std::vector<cfdx::core::FaceIndex> c_faces;

        std::vector<std::vector<cfdx::core::PointIndex>> cell_edges;
        if (nodes.size() == 3) {
            cell_edges = {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[0]}};
        } else if (nodes.size() == 4) {
            cell_edges = {{nodes[0], nodes[1]}, {nodes[1], nodes[2]},
                          {nodes[2], nodes[3]}, {nodes[3], nodes[0]}};
        } else {
            continue;
        }

        for (const auto& edge_nodes : cell_edges) {
            auto canonical = edge_nodes;
            std::sort(canonical.begin(), canonical.end());

            auto it = face_map.find(canonical);
            std::size_t f_idx;
            if (it == face_map.end()) {
                f_idx = all_faces_acc.size();
                all_faces_acc.push_back(edge_nodes);
                face_map[canonical] = {c, f_idx};
            } else {
                f_idx = it->second.second;
            }
            c_faces.push_back(static_cast<cfdx::core::FaceIndex>(f_idx));
        }
        mesh.cells().push_cell(c_faces);
    }

    // Construction finale via build_from_scratch
    mesh.faces().build_from_scratch(all_faces_acc);
    std::cout << "[GMSH_DEBUG] Après build_from_scratch, n_faces = " << mesh.n_faces()
              << ", n_cells = " << mesh.n_cells() << ", all_faces_acc = " << all_faces_acc.size() << "\n";

    // Reconstruction de l'ownership avec la nouvelle numérotation des faces
    mesh.ownership().resize(mesh.n_faces());
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        std::size_t cell_off = mesh.cells().cell_offset(c);
        std::size_t cell_n = mesh.cells().cell_size(c);
        for (std::size_t k = 0; k < cell_n; ++k) {
            std::size_t f = static_cast<std::size_t>(mesh.cells().faces_data()[cell_off + k]);
            if (f < mesh.n_faces()) {
                mesh.ownership().set_owner(f, static_cast<cfdx::core::CellIndex>(c));
                mesh.ownership().set_neighbour(f, -1);  // Par défaut frontière
            }
        }
    }

    std::cout << "[GMSH_DEBUG] Après reconstruction ownership, n_faces = " << mesh.n_faces()
              << ", n_internal = " << mesh.ownership().n_internal_faces()
              << ", n_boundary = " << mesh.ownership().n_boundary_faces() << "\n";
    std::cout << "[GMSH_DEBUG] Calcul géométrie des faces...\n";
    std::vector<cfdx::core::Vec3> face_centres(mesh.n_faces());
    std::vector<cfdx::core::Vec3> face_Sf(mesh.n_faces());
    const auto& pts_data = mesh.points();
    const auto* verts_data = mesh.faces().vertices_data();
    const auto* f_off_data = mesh.faces().offsets_data();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto off = f_off_data[f];
        const auto n = f_off_data[f + 1] - off;
        double cx = 0.0, cy = 0.0, cz = 0.0;
        for (std::size_t v = 0; v < n; ++v) {
            std::size_t node_idx = verts_data[off + v];
            cx += pts_data.x(node_idx); cy += pts_data.y(node_idx); cz += pts_data.z(node_idx);
        }
        face_centres[f] = cfdx::core::Vec3(cx / n, cy / n, cz / n);
        // Calcul du vecteur surface Sf (pour arête 2D : rotation 90° anti-horaire)
        if (n == 2) {
            std::size_t n0 = verts_data[off];
            std::size_t n1 = verts_data[off + 1];
            double dx = pts_data.x(n1) - pts_data.x(n0);
            double dy = pts_data.y(n1) - pts_data.y(n0);
            face_Sf[f] = cfdx::core::Vec3(-dy, dx, 0.0);
        } else {
            face_Sf[f] = cfdx::core::Vec3(0.0, 0.0, 0.0);
        }
    }
    std::cout << "[GMSH_DEBUG] Géométrie des faces calculée: " << face_centres.size() << " faces\n";

    return mesh.n_cells() > 0 || mesh.n_points() > 0;
}

bool import_gmsh_scalar(const std::string&, const std::string&, ScalarCellField&) {
    return false;  // Stub
}
bool import_gmsh_vector(const std::string&, const std::string&, void*) {
    return false;  // Stub — requires VectorCellField definition
}

// ============================================
// Physical Groups parsing (v4 P0 — M0.10-T02)
// ============================================
// Section $PhysicalNames / $PhysicalGroups added to cavity.msh
// Next: map element tags -> boundary patches for CFDX Mesh

} // namespace gmsh
} // namespace io
} // namespace cfdx
