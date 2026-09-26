// M0.10-T03: SU2 adapter implementation
#include "su2_adapter.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace cfdx {
namespace io {
namespace su2 {

// ---------------------------------------------------------------------------
// Helper: trim string
// ---------------------------------------------------------------------------
static std::string trim_s(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static std::string to_lower(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

// ---------------------------------------------------------------------------
// SU2 element type name lookup
// ---------------------------------------------------------------------------
static const char* su2_element_type_name(int type_code) {
    switch (type_code) {
        case 1:  return "NODE";
        case 3:  return "LINE";
        case 5:  return "TRIANGLE";
        case 6:  return "QUADRILATERAL";
        case 9:  return "TETRAHEDRON";
        case 10: return "HEXAHEDRON";
        case 11: return "PRISM";
        case 12: return "PYRAMID";
        default: return "UNKNOWN";
    }
}

static int su2_element_node_count(int type_code) {
    switch (type_code) {
        case 1:  return 1;
        case 3:  return 2;
        case 5:  return 3;
        case 6:  return 4;
        case 9:  return 4;
        case 10: return 8;
        case 11: return 6;
        case 12: return 5;
        default: return 0;
    }
}

static std::string su2_face_type_name(int type_code) {
    switch (type_code) {
        case 3:  return "line";
        case 5:  return "triangle";
        case 6:  return "quad";
        default: return "unknown";
    }
}

// ---------------------------------------------------------------------------
// Source detection — checks for actual SU2 keywords (NPOIN=, NELEM=, NDIME=)
// ---------------------------------------------------------------------------
bool Su2Adapter::detect_source(const std::string& case_path,
                               SourceInfo& info) {
    std::string su2_file = case_path;
    if (su2_file.size() >= 5 && su2_file.substr(su2_file.size() - 5) != ".su2") {
        su2_file += ".su2";
    }

    std::ifstream test(su2_file);
    if (!test.is_open()) return false;

    // Check for SU2-specific section keywords in file content
    std::string content((std::istreambuf_iterator<char>(test)),
                        std::istreambuf_iterator<char>());
    if (content.find("NPOIN") != std::string::npos ||
        content.find("NELEM") != std::string::npos ||
        content.find("NDIME") != std::string::npos) {
        info.solver = "SU2";
        info.case_path = case_path;
        info.case_name = su2_file;
        info.format = "su2";
        info.version = "unknown";
        return true;
    }

    return false;
}

// ---------------------------------------------------------------------------
// Parse .su2 mesh file — actual SU2 format with NDIME=/NELEM=/NPOIN=/NMARK=
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_mesh(const std::string& su2_file) {
    std::ifstream in(su2_file);
    if (!in.is_open()) {
        std::cerr << "Su2Adapter: cannot open .su2 file: " << su2_file << "\n";
        return false;
    }

    su2_points_.clear();
    element_conn_.clear();
    element_types_.clear();
    element_offsets_.clear();
    boundaries_.clear();
    cfg_params_.clear();
    field_names_.clear();

    setup_.source.solver = "SU2";
    setup_.source.case_path = su2_file;
    setup_.source.case_name = su2_file;
    setup_.source.format = "su2";

    std::size_t n_points = 0;
    std::size_t n_elements = 0;
    std::size_t n_markers = 0;
    std::size_t n_dim = 3;

    std::string line;
    while (std::getline(in, line)) {
        std::string trimmed = trim_s(line);
        if (trimmed.empty() || trimmed[0] == '%') continue;

        std::string lower = to_lower(trimmed);

        // NDIME= <n>
        if (lower.rfind("ndime", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string key; iss >> key;
            iss >> n_dim;
            setup_.mesh_info.dimension = static_cast<int>(n_dim);
            continue;
        }

        // NELEM= <n>
        if (lower.rfind("nelem", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string key; iss >> key;
            iss >> n_elements;
            parse_su2_elements(in, n_elements);
            continue;
        }

        // NPOIN= <n>
        if (lower.rfind("npoin", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string key; iss >> key;
            iss >> n_points;
            parse_su2_points(in, n_points);
            continue;
        }

        // NMARK= <n>
        if (lower.rfind("nmark", 0) == 0) {
            std::istringstream iss(trimmed);
            std::string key; iss >> key;
            iss >> n_markers;
            parse_su2_boundary(in, n_markers);
            continue;
        }

        // Skip unrecognized lines (comments, etc.)
    }

    setup_.mesh_info.n_vertices = su2_points_.size();
    if (n_dim <= 2) setup_.mesh_info.dimension = 2;

    // Collect unique cell types
    std::set<std::string> types_found;
    for (int t : element_types_) {
        types_found.insert(su2_element_type_name(t));
    }
    setup_.mesh_info.cell_types.assign(types_found.begin(), types_found.end());
    setup_.mesh_info.has_polyhedral = types_found.count("POLYHEDRON") > 0;

    return !su2_points_.empty();
}

// ---------------------------------------------------------------------------
// SU2 .su2 POINTS section (NPOIN=)
// Each point line: x y z [point_id]
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_su2_points(std::istream& in, std::size_t n_points) {
    su2_points_.clear();
    su2_points_.reserve(n_points);

    for (std::size_t i = 0; i < n_points; ++i) {
        cfdx::core::Vec3 p;
        in >> p.x >> p.y >> p.z;
        // Read and discard the point index (if present)
        int dummy_id;
        if (in >> dummy_id) {
            // Successfully read the index — nothing more to do
        } else {
            // No index field — seek back
            in.clear();
        }
        su2_points_.push_back(p);
    }

    return true;
}

// ---------------------------------------------------------------------------
// SU2 .su2 ELEMENTS section (NELEM=)
// Each element line: <type> <n1> <n2> ... <nn> <cell_id>
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_su2_elements(std::istream& in, std::size_t n_elements) {
    element_conn_.clear();
    element_types_.clear();
    element_offsets_.clear();

    for (std::size_t i = 0; i < n_elements; ++i) {
        int elem_type;
        in >> elem_type;

        int n_nodes = su2_element_node_count(elem_type);
        if (n_nodes == 0) {
            std::string line;
            std::getline(in, line);
            continue;
        }

        element_offsets_.push_back(element_conn_.size());
        element_types_.push_back(elem_type);

        for (int n = 0; n < n_nodes; ++n) {
            std::uint32_t node_id;
            in >> node_id;
            // SU2 uses 0-based indexing
            element_conn_.push_back(node_id);
        }

        // Read and discard the cell_id (if present)
        int dummy_cell_id;
        if (in >> dummy_cell_id) {
            // Successfully read
        } else {
            in.clear();
        }
    }

    element_offsets_.push_back(element_conn_.size());
    return true;
}

// ---------------------------------------------------------------------------
// SU2 .su2 boundary section (NMARK=)
// Format:
//   MARKER_TAG <name>
//   MARKER_ELEMS <n_faces>
//   <type> <n1> <n2> ...
//   [repeat for each face]
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_su2_boundary(std::istream& in, std::size_t n_markers) {
    boundaries_.clear();
    boundaries_.reserve(n_markers);

    for (std::size_t b = 0; b < n_markers; ++b) {
        std::string token;
        in >> token;

        // Expect "MARKER_TAG"
        if (to_lower(token) != "marker_tag") {
            // Try reading as marker name directly
        }

        std::string marker_name;
        in >> marker_name;

        // Remove quotes if present
        if (!marker_name.empty() && marker_name.front() == '"')
            marker_name.erase(0, 1);
        if (!marker_name.empty() && marker_name.back() == '"')
            marker_name.pop_back();

        // Expect "MARKER_ELEMS"
        in >> token;  // should be "MARKER_ELEMS"
        std::size_t n_face_elements;
        in >> n_face_elements;

        Su2Boundary bound;
        bound.name = marker_name;
        bound.n_elements = n_face_elements;

        // Read face connectivity
        for (std::size_t i = 0; i < n_face_elements; ++i) {
            int face_type;
            in >> face_type;

            int n_nodes = su2_element_node_count(face_type);
            bound.element_type = su2_face_type_name(face_type);

            for (int n = 0; n < n_nodes; ++n) {
                std::uint32_t node_id;
                in >> node_id;
                bound.node_ids.push_back(node_id);
            }
        }

        boundaries_.push_back(bound);
    }

    return true;
}

// ---------------------------------------------------------------------------
// Parse .cfg key = value
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_cfg_value(const std::string& line,
                                 std::string& key,
                                 std::string& val) {
    size_t eq = line.find('=');
    if (eq == std::string::npos) return false;

    key = trim_s(line.substr(0, eq));
    val = trim_s(line.substr(eq + 1));

    // Remove inline comments
    size_t comma = val.find(';');
    if (comma != std::string::npos) val = val.substr(0, comma);
    val = trim_s(val);

    return true;
}

// ---------------------------------------------------------------------------
// Parse MARKER_BOUNDARY condition from config
// ---------------------------------------------------------------------------
void Su2Adapter::parse_cfg_marker_bc(const std::string& value) {
    // Format: ( marker_name, bc_type )
    // e.g. ( wall, euler_wall )
    std::string v = trim_s(value);
    if (v.empty()) return;

    // Remove parentheses
    if (!v.empty() && v.front() == '(') v.erase(0, 1);
    if (!v.empty() && v.back() == ')') v.pop_back();

    // Parse pairs
    std::istringstream iss(v);
    std::string token;
    std::vector<std::string> tokens;
    while (std::getline(iss, token, ',')) {
        token = trim_s(token);
        if (!token.empty()) tokens.push_back(token);
    }

    for (size_t i = 0; i + 1 < tokens.size(); i += 2) {
        const std::string& marker_name = tokens[i];
        const std::string& bc_type = tokens[i + 1];

        BoundarySpec bc;
        bc.patch_name = marker_name;
        bc.source_zone_name = marker_name;
        bc.source_type_name = bc_type;

        // Map SU2 BC types
        std::string lower_bc = to_lower(bc_type);
        std::string mapped_turb = rules_.map_turbulence(bc_type);
        bc.type = rules_.map_boundary_type(marker_name);
        bc.value_type = rules_.map_bc_value_type(bc_type);

        if (lower_bc.find("far") != std::string::npos) {
            bc.type = BCType::OUTLET;
            bc.value_type = BCValueType::FIXED;
        } else if (lower_bc.find("wall") != std::string::npos) {
            bc.type = BCType::WALL;
            bc.value_type = BCValueType::WALL_NO_SLIP;
        } else if (lower_bc.find("inlet") != std::string::npos) {
            bc.type = BCType::INLET;
            bc.value_type = BCValueType::FIXED;
        } else if (lower_bc.find("outlet") != std::string::npos) {
            bc.type = BCType::PRESSURE_OUTLET;
            bc.value_type = BCValueType::OUTLET_PRESSURE;
        } else if (lower_bc.find("symmetr") != std::string::npos) {
            bc.type = BCType::SYMMETRY;
            bc.value_type = BCValueType::ZERO_GRADIENT;
        }

        setup_.boundary_conditions.push_back(bc);

        // Also update matching Su2Boundary if present
        for (auto& bnd : boundaries_) {
            if (bnd.name == marker_name) {
                bnd.name = marker_name;  // already set
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Parse numerical scheme settings from cfg
// ---------------------------------------------------------------------------
void Su2Adapter::parse_cfg_numerics() {
    // MUSCL
    if (cfg_params_.count("MUSCL")) {
        std::string val = cfg_params_["MUSCL"];
        if (to_lower(val) == "yes" || to_lower(val) == "true") {
            setup_.numerics.momentum_scheme = "MUSCL";
        }
    }

    // Scheme type
    if (cfg_params_.count("NUM_METHOD_GRAD")) {
        std::string val = cfg_params_["NUM_METHOD_GRAD"];
        std::string mapped = rules_.map_scheme(val);
        if (!mapped.empty()) {
            setup_.numerics.gradient_operator = mapped;
        } else {
            setup_.numerics.gradient_operator = "least_squares";
        }
    }

    // Convergence
    if (cfg_params_.count("CONV_NUMB")) {
        try {
            setup_.numerics.max_iterations = std::stoi(cfg_params_["CONV_NUMB"]);
        } catch (...) {}
    }
    if (cfg_params_.count("RESIDUAL_REDUCTION")) {
        setup_.numerics.residual_target = cfg_params_["RESIDUAL_REDUCTION"];
    }
    if (cfg_params_.count("TIME_DISCRE_FEM")) {
        std::string val = cfg_params_["TIME_DISCRE_FEM"];
        if (!val.empty()) {
            setup_.numerics.transient_scheme = val;
            if (to_lower(val).find("runge") != std::string::npos) {
                setup_.transient = true;
            }
        }
    }

    setup_.numerics.momentum_interpolation = "linear";
    setup_.numerics.pressure_scheme = "standard";
    setup_.numerics.coupled_solver = "coupled";

    // Record unmapped settings
    for (const auto& [k, v] : cfg_params_) {
        setup_.numerics.raw_settings[k] = v;
    }
}

// ---------------------------------------------------------------------------
// Parse .cfg config file
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_config(const std::string& cfg_file) {
    std::ifstream in(cfg_file);
    if (!in.is_open()) {
        std::cerr << "Su2Adapter: cannot open .cfg file: " << cfg_file << "\n";
        return false;
    }

    setup_.source.case_path = cfg_file;

    std::string line;
    while (std::getline(in, line)) {
        // Skip comments and empty lines
        std::string trimmed = trim_s(line);
        if (trimmed.empty() || trimmed[0] == '%') continue;

        std::string key, val;
        if (parse_cfg_value(trimmed, key, val)) {
            cfg_params_[key] = val;

            // Map known keys
            std::string lower_key = to_lower(key);

            if (lower_key.find("mach") != std::string::npos) {
                try { setup_.initial_condition.velocity = std::stod(val) * 340.0; } catch (...) {}
            }
            if (lower_key.find("alpha") != std::string::npos &&
                lower_key.find("aoa") != std::string::npos) {
                try {
                    double alpha = std::stod(val) * M_PI / 180.0;
                    setup_.initial_condition.velocity_vector = {
                        std::cos(alpha), std::sin(alpha), 0.0
                    };
                } catch (...) {}
            }
            if (lower_key.find("marker_boundary") != std::string::npos &&
                lower_key.find("marker") != std::string::npos) {
                parse_cfg_marker_bc(val);
            }
            if (lower_key == "fluid_prandtl") {
                try {
                    MaterialSpec mat;
                    mat.name = "fluid";
                    mat.specific_heat = std::stod(val);
                    setup_.materials.push_back(mat);
                } catch (...) {}
            }
            if (lower_key == "freestream_pressure") {
                try { setup_.initial_condition.pressure = std::stod(val); } catch (...) {}
            }
            if (lower_key == "freestream_temperature") {
                try { setup_.initial_condition.temperature = std::stod(val); } catch (...) {}
            }
            if (lower_key == "visc") {
                // Viscosity
                try {
                    if (!setup_.materials.empty()) {
                        setup_.materials.back().dynamic_viscosity = std::stod(val);
                    }
                } catch (...) {}
            }
            if (lower_key.find("turbulence") != std::string::npos ||
                lower_key.find("model") != std::string::npos) {
                std::string mapped = rules_.map_turbulence(val);
                if (!mapped.empty()) {
                    setup_.turbulence_model = mapped;
                } else if (to_lower(val).find("sa") != std::string::npos) {
                    setup_.turbulence_model = "spalart_allmaras";
                } else if (to_lower(val).find("none") != std::string::npos) {
                    setup_.turbulence_model = "laminar";
                }
            }
            if (lower_key == "energy") {
                if (to_lower(val) == "yes" || to_lower(val) == "true") {
                    setup_.energy_model = "buoyant";
                } else {
                    setup_.energy_model = "isothermal";
                }
            }
            if (lower_key.find("time") != std::string::npos &&
                lower_key.find("stepping") != std::string::npos) {
                try { setup_.transient = true; } catch (...) {}
            }
        }
    }

    // Parse numerics from collected params
    parse_cfg_numerics();

    // Detect unknown/unmapped settings
    std::vector<std::string> known_keys = {
        "MUSCL", "NUM_METHOD_GRAD", "CONV_NUMB", "RESIDUAL_REDUCTION",
        "TIME_DISCRE_FEM", "MACH_NUMBER", "FREESTREAM_PRESSURE",
        "FREESTREAM_TEMPERATURE", "FLUID_PRANDTL", "VISCO",
        "MARKER_BOUNDARY", "MARKER_MONITOR", "MARKER_PLOT"
    };

    // Store all raw params — used for gap analysis
    for (const auto& [k, v] : cfg_params_) {
        bool is_known = false;
        for (const auto& known : known_keys) {
            if (to_lower(k).find(to_lower(known)) != std::string::npos) {
                is_known = true;
                break;
            }
        }
        if (!is_known) {
            setup_.numerics.raw_settings[k] = v;
        }
    }

    // Extract field names from MARKER_MONITOR or solution data
    if (cfg_params_.count("MARKER_MONITOR")) {
        std::string monitors = cfg_params_["MARKER_MONITOR"];
        // Parse field names: ( (x,y,z), marker, ... )
        std::istringstream iss(monitors);
        std::string tok;
        while (std::getline(iss, tok, ',')) {
            tok = trim_s(tok);
            tok.erase(std::remove(tok.begin(), tok.end(), '('), tok.end());
            tok.erase(std::remove(tok.begin(), tok.end(), ')'), tok.end());
            if (!tok.empty()) field_names_.push_back(tok);
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// Parse solution.csv for field data
// ---------------------------------------------------------------------------
bool Su2Adapter::parse_solution_csv(const std::string& csv_file,
                                    const std::string& field_name,
                                    cfdx::core::ScalarCellField& scalar_field,
                                    cfdx::core::Vec3CellField& vec_field) {
    std::ifstream in(csv_file);
    if (!in.is_open()) return false;

    std::string line;
    std::getline(in, line);  // header

    std::vector<double> values;
    double val;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        // First column is typically index/time
        int idx;
        std::string comma;
        if (!(iss >> idx >> comma)) continue;

        // Find the requested field column
        while (iss >> val) {
            values.push_back(val);
            char c;
            if (iss >> c && c != ',') {
                iss.unget();
            }
        }
    }

    if (values.empty()) return false;

    // If field_name suggests a vector field (e.g., "Velocity"), fill Vec3CellField
    std::string lower_fn = to_lower(field_name);
    if (lower_fn.find("vel") != std::string::npos || lower_fn.find("u") != std::string::npos) {
        std::size_t n_cells = values.size() / 3;
        if (n_cells * 3 != values.size()) return false;
        vec_field = cfdx::core::Vec3CellField(n_cells, field_name, "m/s", 3);
        for (std::size_t i = 0; i < n_cells; ++i) {
            vec_field.set(i, values[i * 3], values[i * 3 + 1], values[i * 3 + 2]);
        }
    } else {
        scalar_field = cfdx::core::ScalarCellField(values.size(), field_name, "Pa", 1);
        for (std::size_t i = 0; i < values.size(); ++i) {
            scalar_field(i) = values[i];
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// Results-only import
// ---------------------------------------------------------------------------
bool Su2Adapter::import_results(const std::string& results_path,
                                ConversionResult& result) {
    cfdx::core::ScalarCellField pressure;
    cfdx::core::Vec3CellField velocity;

    if (!parse_solution_csv(results_path, "Pressure", pressure, velocity)) {
        result.gap_report.unsupported_nonblocking(
            "results", "solution_csv",
            "Cannot parse SU2 solution.csv: " + results_path,
            "Ensure the file is a valid CSV with numerical data columns");
        return false;
    }

    result.gap_report.supported(
        "results", "solution_csv",
        "Imported SU2 solution.csv: " + results_path);

    if (pressure.size() > 0) {
        result.scalar_fields.push_back(pressure);
    }
    if (velocity.size() > 0) {
        result.vec_fields.push_back(velocity);
    }

    return true;
}

// ---------------------------------------------------------------------------
// Full conversion entry point
// ---------------------------------------------------------------------------
bool Su2Adapter::convert(const std::string& case_path,
                         ConversionResult& result) {
    // Determine file names
    std::string su2_file = case_path;
    std::string cfg_file = case_path;

    // Handle directory vs. file path
    if (su2_file.size() > 5 && su2_file.substr(su2_file.size() - 5) == ".su2") {
        // Direct file path
    } else if (cfg_file.size() > 4 && cfg_file.substr(cfg_file.size() - 4) == ".cfg") {
        // Direct cfg path
    } else {
        // Assume directory
        su2_file = case_path + "/mesh.su2";
        cfg_file = case_path + "/config.cfg";
    }

    bool mesh_ok = parse_mesh(su2_file);
    if (!mesh_ok) {
        result.gap_report.unsupported_blocking(
            "mesh", "su2_mesh",
            "Cannot parse SU2 mesh file: " + su2_file,
            "Ensure the .su2 file is a valid SU2 mesh format");
        return false;
    }

    bool cfg_ok = parse_config(cfg_file);
    if (!cfg_ok) {
        result.gap_report.unsupported_nonblocking(
            "config", "su2_cfg",
            "Cannot parse SU2 config file: " + cfg_file,
            "Config parsing is optional; mesh will still be imported");
    }

    result.source = setup_.source;

    // Populate mesh on result
    result.mesh.clear();
    result.mesh.points().resize(su2_points_.size());
    for (std::size_t i = 0; i < su2_points_.size(); ++i) {
        result.mesh.points().set(i, su2_points_[i].x, su2_points_[i].y, su2_points_[i].z);
    }

    // Build face topology from elements
    // For SU2, each element (tet, hex, prism, etc.) defines a cell;
    // we need to extract faces and build owner/neighbour
    // Simplified: triangulate all faces and build CSR connectivity

    // Build faces from elements
    std::vector<std::vector<cfdx::core::FaceIndex>> element_faces;
    for (std::size_t i = 0; i + 1 < element_offsets_.size(); ++i) {
        std::size_t off = element_offsets_[i];
        std::size_t n = element_offsets_[i + 1] - off;
        int elem_type = element_types_[i];

        // Get face connectivity for this element type
        // For tetra: 4 triangular faces
        // For hexa: 6 quadrilateral faces
        // etc.
        std::vector<std::vector<cfdx::core::FaceIndex>> faces;
        int n_nodes = su2_element_node_count(elem_type);

        if (elem_type == 9) {  // TETRAHEDRON
            // 4 triangular faces (consistent outward orientation)
            faces = {
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3])}
            };
        } else if (elem_type == 10) {  // HEXAHEDRON
            // 6 quadrilateral faces
            faces = {
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 4]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 5]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 6]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 7])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 5]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 4])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 6]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 5])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 7]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 6])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 4]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 7]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3])}
            };
        } else if (elem_type == 11) {  // PRISM
            faces = {
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 4]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 5])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 4])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 1]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 5]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 4])},
                {static_cast<cfdx::core::FaceIndex>(element_conn_[off + 0]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 2]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 5]),
                 static_cast<cfdx::core::FaceIndex>(element_conn_[off + 3])},
            };
        } else {
            // Skip unsupported element types
            continue;
        }

        for (const auto& face : faces) {
            // Canonical face (sorted vertex set for matching)
            std::vector<cfdx::core::FaceIndex> canonical = face;
            std::sort(canonical.begin(), canonical.end());

            bool found = false;
            for (std::size_t ef = 0; ef < element_faces.size(); ++ef) {
                std::vector<cfdx::core::FaceIndex> canonical_ef = element_faces[ef];
                std::sort(canonical_ef.begin(), canonical_ef.end());
                if (canonical_ef == canonical) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                element_faces.push_back(face);
            }
        }
    }

    // Build mesh faces and ownership
    mesh.faces().build_from_scratch(element_faces);

    result.mesh.ownership().resize(result.mesh.n_faces());
    // Build ownership from boundary info and element connectivity
    std::fill_n(result.mesh.ownership().neighbour_data(),
                result.mesh.n_faces(),
                cfdx::core::FaceOwnership::BOUNDARY);

    // Set owners from cell faces
    for (std::size_t c = 0; c < element_faces.size() && c < result.mesh.cells().n_cells(); ++c) {
        for (std::size_t k = 0; k < result.mesh.cells().cell_size(c); ++k) {
            auto f = result.mesh.cells().faces_data()[result.mesh.cells().cell_offset(c) + k];
            if (f < result.mesh.n_faces()) {
                result.mesh.ownership().set_owner(f, static_cast<cfdx::core::CellIndex>(c));
                result.mesh.ownership().set_neighbour(f, -1);  // boundary by default
            }
        }
    }

    // Boundary patches from SU2 boundary markers
    cfdx::core::BoundaryPatches bp;
    std::size_t face_offset = 0;
    for (const auto& bnd : boundaries_) {
        if (bnd.is_internal()) continue;

        cfdx::core::Patch p;
        p.name = bnd.name;
        std::string lc_name = to_lower(bnd.name);
        p.type = static_cast<cfdx::core::PatchType>(rules_.map_patch_type(lc_name));
        if (p.type == static_cast<cfdx::core::PatchType>(-1)) {
            p.type = cfdx::core::PatchType::UNKNOWN;
        }

        // Add face IDs (simplified — would need proper face matching in full impl)
        for (std::size_t k = 0; k < bnd.n_elements && face_offset < result.mesh.n_faces(); ++k) {
            p.face_ids.push_back(static_cast<std::uint32_t>(face_offset));
            face_offset++;
        }

        bp.add_patch(p);
    }
    result.mesh.set_boundary(bp);

    // --- Gap Analysis ---

    result.gap_report.supported(
        "mesh", "points",
        "Parsed " + std::to_string(su2_points_.size()) + " vertices from .su2 NPOIN section");

    result.gap_report.supported(
        "mesh", "elements",
        "Parsed " + std::to_string(element_types_.size()) + " elements from .su2 NELEM section");

    for (const auto& bnd : boundaries_) {
        result.gap_report.supported(
            "boundary", bnd.name,
            "Detected boundary marker '" + bnd.name + "' with " +
            std::to_string(bnd.n_elements) + " faces");
    }

    for (const auto& mat : setup_.materials) {
        result.gap_report.supported(
            "material", mat.name,
            "Parsed fluid material from config");
    }

    for (const auto& bc : setup_.boundary_conditions) {
        result.gap_report.supported(
            "boundary_condition", bc.patch_name,
            "Mapped SU2 marker '" + bc.source_zone_name + "' to CFDX BC type '" +
            to_string(bc.type) + "'");
    }

    result.gap_report.supported(
        "numerics", "scheme",
        "Mapped SU2 numerical scheme to CFDX: " + setup_.numerics.momentum_scheme);

    // Unsupported features
    result.gap_report.unsupported_nonblocking(
        "features", "turbomachinery",
        "SU2 turbomachinery extensions not supported",
        "Use base SU2 configuration without turbo machinery options");

    result.gap_report.unsupported_nonblocking(
        "features", "multiphysics",
        "SU2 multi-physics extensions not supported",
        "Use single-physics SU2 configuration");

    result.gap_report.unsupported_nonblocking(
        "features", "deforming_mesh",
        "SU2 dynamic mesh options not supported",
        "Use static mesh for CFDX conversion");

    // Check for unmapped raw settings
    if (!setup_.numerics.raw_settings.empty()) {
        for (const auto& [k, v] : setup_.numerics.raw_settings) {
            result.gap_report.unsupported_nonblocking(
                "config", k,
                "Unmapped SU2 config parameter: " + k + " = " + v,
                "Manually verify the parameter is applied in CFDX setup");
        }
    }

    return !result.gap_report.has_blocking();
}

// ---------------------------------------------------------------------------
// Standalone convenience functions
// ---------------------------------------------------------------------------
bool parse_su2_mesh(const std::string& su2_file,
                    Su2Adapter& adapter,
                    cfdx::core::Mesh& mesh,
                    GapAnalysis& gap,
                    MappingRules rules) {
    Su2Adapter local_adapter(rules);
    ConversionResult result;

    local_adapter.detect_source(su2_file, result.source);

    if (!local_adapter.parse_mesh(su2_file)) {
        return false;
    }

    local_adapter.convert(su2_file, result);
    gap = result.gap_report;

    // Transfer mesh
    mesh = std::move(result.mesh);
    adapter = std::move(local_adapter);
    return true;
}

bool parse_su2_config(const std::string& cfg_file,
                      Su2Adapter& adapter,
                      CaseSetup& setup,
                      GapAnalysis& gap) {
    if (!adapter.parse_config(cfg_file)) {
        gap.unsupported_blocking(
            "config", "su2_cfg",
            "Cannot parse SU2 config: " + cfg_file,
            "Ensure the file is a valid SU2 .cfg configuration");
        return false;
    }

    // Copy settings
    setup = adapter.setup();
    gap.supported("config", "su2_cfg", "Parsed SU2 config: " + cfg_file);
    return true;
}

}  // namespace su2
}  // namespace io
}  // namespace cfdx
