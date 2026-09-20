#include "cfdx/core/field/field.h"
#include "gmsh_importer.h"
#include "cfdx/core/mesh/mesh.h"
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

using namespace cfdx::io::gmsh;
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
    if (!file.is_open()) return false;

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string full = buffer.str();

    std::map<uint32_t, cfdx::core::Vec3> points_map;
    // Find Nodes section in full file
    size_t nodes_start = full.find("$Nodes");
    size_t nodes_end = full.find("$EndNodes");
    size_t elements_start = full.find("$Elements");
    size_t elements_end = full.find("$EndElements");

    if (nodes_start != std::string::npos && nodes_end != std::string::npos) {
        std::string nodes_content = full.substr(nodes_start, nodes_end - nodes_start);
        parseNodes(nodes_content, points_map);
    }

    mesh.points().resize(points_map.size());
    std::map<uint32_t, uint32_t> point_index_map;  // Gmsh node ID -> CFDX index
    uint32_t idx = 0;
    for (const auto& kv : points_map) {
        mesh.points().set(idx, kv.second.x, kv.second.y, kv.second.z);
        point_index_map[kv.first] = idx;
        idx++;
    }

    if (elements_start != std::string::npos && elements_end != std::string::npos) {
        std::string elements_content = full.substr(elements_start, elements_end - elements_start);
        parseElements(elements_content, mesh, points_map);
    }

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
