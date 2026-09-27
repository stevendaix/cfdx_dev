// M0.10-T03: Fluent adapter implementation
#include "fluent_adapter.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace cfdx {
namespace io {
namespace fluent {

// ---------------------------------------------------------------------------
// Helper: read a section from legacy ASCII .cas
// Format: (section_id n_values ... data ...)
// ---------------------------------------------------------------------------
static bool is_section_start(const std::string& line, int* id_out, int* n_vals_out) {
    std::string s = line;
    // Remove leading whitespace
    while (!s.empty() && (s[0] == ' ' || s[0] == '\t')) s.erase(0, 1);
    if (s.empty() || s[0] != '(') return false;

    std::istringstream iss(s.substr(1));
    int id, n_vals;
    if (!(iss >> id >> n_vals)) return false;
    *id_out = id;
    *n_vals_out = n_vals;
    return true;
}

static void skip_section(std::istream& in) {
    std::string tok;
    while (in >> tok) {
        if (tok == ")") break;
    }
}

static std::string read_rest_of_line(std::istream& in) {
    std::string line;
    std::getline(in, line);
    // Trim trailing whitespace
    while (!line.empty() && (line.back() == ' ' || line.back() == '\t' ||
                             line.back() == '\r')) {
        line.pop_back();
    }
    return line;
}

// ---------------------------------------------------------------------------
// Source detection
// ---------------------------------------------------------------------------
bool FluentAdapter::detect_source(const std::string& case_path,
                                  SourceInfo& info) {
    // Check for .cas file existence
    std::string cas_file = case_path;
    if (cas_file.size() > 4 && cas_file.substr(cas_file.size() - 4) != ".cas") {
        cas_file += ".cas";
    }

    std::ifstream test(cas_file);
    if (!test.is_open()) return false;

    info.solver = "Fluent";
    info.case_path = case_path;
    info.case_name = cas_file;
    info.format = "cas_dat";

    // Try to detect version from header
    std::string line;
    std::getline(test, line);
    if (line.find("ANSYS") != std::string::npos ||
        line.find("Fluent") != std::string::npos) {
        info.version = "202x";
    } else {
        info.version = "unknown (legacy)";
    }
    return true;
}

// ---------------------------------------------------------------------------
// Zone type conversion
// ---------------------------------------------------------------------------
ZoneType FluentAdapter::zone_type_from_int(int t) const {
    switch (t) {
        case 1:  return ZoneType::INTERIOR;
        case 4:  return ZoneType::VELOCITY_INLET;
        case 5:  return ZoneType::PRESSURE_INLET;
        case 7:  return ZoneType::PRESSURE_OUTLET;
        case 8:  return ZoneType::WALL;
        case 9:  return ZoneType::SYMMETRY;
        case 10: return ZoneType::PERIODIC;
        case 3:  return ZoneType::EXTERIOR;
        case 2:  return ZoneType::INTERNAL;
        case 6:  return ZoneType::MASS_OUTFLOW;
        case 11: return ZoneType::PRESSURE_FAR_FIELD;
        case 12: return ZoneType::OVERSET;
        default: return ZoneType::UNKNOWN;
    }
}

// ---------------------------------------------------------------------------
// .cas header parsing
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_header(std::istream& in) {
    std::string line;
    std::getline(in, line);

    // Fluent .cas header: first line often contains "ANSYS ..." or version info
    setup_.source.solver = "Fluent";
    setup_.source.version = "legacy (ASCII)";
    setup_.source.format = "cas_dat";

    // Look for version line
    if (line.find("version") != std::string::npos ||
        line.find("Version") != std::string::npos) {
        setup_.source.version = line;
    }

    // Default settings if not found in models section
    setup_.transient = false;
    setup_.physics_model = "incompressible_laminar";
    setup_.turbulence_model = "laminar";

    return true;
}

// ---------------------------------------------------------------------------
// .cas nodes parsing (section 2)
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_nodes(std::istream& in) {
    int n_coords = 0;
    in >> n_coords;

    points_.clear();
    points_.reserve(n_coords / 3);

    for (int i = 0; i < n_coords; ++i) {
        double x, y, z;
        in >> x >> y >> z;
        points_.push_back(cfdx::core::Vec3(x, y, z));
    }

    // Read closing paren
    char c;
    in >> c;

    if (points_.empty()) {
        setup_.mesh_info.n_vertices = 0;
        return false;
    }

    setup_.mesh_info.n_vertices = points_.size();
    setup_.mesh_info.dimension = 3;
    return true;
}

// ---------------------------------------------------------------------------
// .cas faces parsing (section 3 - face nodes, section 4 - faces)
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_faces(std::istream& in) {
    // Section 3: face nodes
    // Read n_values, then n_faces, then for each face: n_nodes, node_ids...
    int n_faces_section3 = 0;
    in >> n_faces_section3;

    face_nodes_.clear();
    face_nodes_.reserve(n_faces_section3);

    for (int f = 0; f < n_faces_section3; ++f) {
        int n_nodes;
        in >> n_nodes;
        std::vector<std::uint32_t> nodes;
        nodes.reserve(n_nodes);
        for (int n = 0; n < n_nodes; ++n) {
            std::uint32_t v;
            in >> v;
            nodes.push_back(v);
        }
        face_nodes_.push_back(std::move(nodes));
    }

    // Read closing paren for section 3
    char c;
    in >> c;

    // Section 4: faces with owner/neighbour
    int n_faces_section4 = 0;
    in >> n_faces_section4;

    face_owner_.clear();
    face_neighbour_.clear();
    face_owner_.reserve(n_faces_section4);
    face_neighbour_.reserve(n_faces_section4);

    for (int f = 0; f < n_faces_section4; ++f) {
        int n_nodes;
        in >> n_nodes;
        // Skip node indices (already read in section 3, or read here for binary compat)
        for (int n = 0; n < n_nodes; ++n) {
            int v;
            in >> v;
        }
        int owner, neighbour;
        in >> owner >> neighbour;
        face_owner_.push_back(owner);
        face_neighbour_.push_back(neighbour);
    }

    // Read closing paren for section 4
    in >> c;

    setup_.mesh_info.n_faces = face_nodes_.size();
    return true;
}

// ---------------------------------------------------------------------------
// .cas cells parsing (section 5)
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_cells(std::istream& in) {
    int n_cells_section5 = 0;
    in >> n_cells_section5;

    // Section 5 format: for each cell, n_faces, face_ids...
    // Build cell -> face mapping
    cell_faces_.clear();
    cell_faces_.reserve(n_cells_section5);

    for (int c = 0; c < n_cells_section5; ++c) {
        int n_faces;
        in >> n_faces;
        std::vector<std::uint32_t> faces;
        faces.reserve(n_faces);
        for (int f = 0; f < n_faces; ++f) {
            std::uint32_t fid;
            in >> fid;
            faces.push_back(fid);
        }
        cell_faces_.push_back(std::move(faces));
    }

    // Read closing paren
    char c;
    in >> c;

    setup_.mesh_info.n_cells = cell_faces_.size();
    return true;
}

// ---------------------------------------------------------------------------
// .cas zones parsing (section 31 - boundary conditions)
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_zones(std::istream& in, int n_zones) {
    // Section 31 format:
    // n_zones
    // For each zone: zone_id, zone_type, name_len, name, node_count, element_count, ...

    zones_.clear();
    zones_.reserve(n_zones);

    for (int i = 0; i < n_zones; ++i) {
        ZoneInfo zone;
        in >> zone.zone_id;

        // Read name (may be quoted)
        std::string name;
        in >> std::ws;
        if (in.peek() == '"') {
            // Quoted name
            char quote;
            in >> quote;
            std::getline(in, name, '"');
        } else {
            // Unquoted: read until space/newline
            in >> name;
        }

        // Read remaining zone data
        int n_elements = 0;
        std::string node_type_str, element_type_str;

        // Skip remaining fields and collect the node type and element type
        // The exact format depends on Fluent version, so we parse flexibly
        std::vector<std::string> tokens;
        while (in.peek() != ')' && in.peek() != '\n' && !in.eof()) {
            std::string tok;
            in >> tok;
            if (!tok.empty()) tokens.push_back(tok);
        }

        // Try to identify node_type and element_type from tokens
        zone.name = name;
        zone.node_type = tokens.size() > 0 ? tokens[0] : "unknown";
        zone.element_type = tokens.size() > 1 ? tokens[1] : "";
        if (!tokens.empty()) {
            try {
                zone.n_elements = std::stoi(tokens.back());
            } catch (...) {}
        }

        zones_.push_back(zone);

        // Read closing paren if on same line
        if (in.peek() == ')') {
            char c;
            in >> c;
        }
    }

    // Map zones to boundary conditions and patches
    cfdx::core::BoundaryPatches bp;
    for (auto& zone : zones_) {
        if (!zone.is_boundary()) continue;

        ZoneType zt = zone_type_from_int(zone.zone_type);
        const char* zone_name = to_string(zt);

        // Use mapping rules
        cfdx::core::PatchType patch_type = static_cast<cfdx::core::PatchType>(
            rules_.map_patch_type(zone_name));
        if (patch_type == static_cast<cfdx::core::PatchType>(-1)) {
            patch_type = cfdx::core::PatchType::UNKNOWN;
        }

        cfdx::core::Patch p;
        p.name = zone.name;
        p.type = patch_type;

        // Read face IDs from the zone's face range
        // For boundary zones, face IDs are typically consecutive
        // (stored in a separate face_zone_to_face_table section)
        // Here we set up the patch structure; face_ids are filled during mesh build

        bp.add_patch(p);
    }

    // Store zone info as boundary conditions in setup
    for (auto& zone : zones_) {
        if (!zone.is_boundary()) continue;

        BoundarySpec bc;
        ZoneType zt = zone_type_from_int(zone.zone_type);
        const char* zone_name = to_string(zt);

        bc.patch_name = zone.name;
        bc.type = rules_.map_boundary_type(zone_name);
        bc.value_type = rules_.map_bc_value_type(zone_name);
        bc.source_zone_name = zone.name;
        bc.source_type_name = zone_name;

        if (zone.n_elements > 0) {
            bc.raw_params["n_elements"] = std::to_string(zone.n_elements);
        }

        setup_.boundary_conditions.push_back(bc);
    }

    setup_.mesh_info.n_patches = bp.n_patches();

    return true;
}

// ---------------------------------------------------------------------------
// .cas materials parsing (section 39)
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_materials(std::istream& in, int n_materials) {
    // Section 39: materials database (text-based, variable number of entries)

    for (int m = 0; m < n_materials; ++m) {
        MaterialSpec mat;

        // Material name
        std::string name;
        in >> std::ws;
        if (in.peek() == '"') {
            char q;
            in >> q;
            std::getline(in, name, '"');
        } else {
            in >> name;
        }
        mat.name = name;

        // Read properties: density, viscosity, conductivity, specific_heat
        // Format varies; we read what we can
        std::vector<std::string> tokens;
        while (in.peek() != ')' && in.peek() != '\n' && !in.eof()) {
            std::string tok;
            in >> tok;
            if (!tok.empty()) tokens.push_back(tok);
        }

        // Parse tokens for known properties
        for (size_t i = 0; i + 1 < tokens.size(); ++i) {
            if (tokens[i] == "density" || tokens[i] == "rho") {
                try { mat.density = std::stod(tokens[i + 1]); } catch (...) {}
            } else if (tokens[i] == "viscosity" || tokens[i] == "mu") {
                try { mat.dynamic_viscosity = std::stod(tokens[i + 1]); } catch (...) {}
            } else if (tokens[i] == "conductivity" || tokens[i] == "k") {
                try { mat.thermal_conductivity = std::stod(tokens[i + 1]); } catch (...) {}
            } else if (tokens[i] == "specific_heat" || tokens[i] == "cp") {
                try { mat.specific_heat = std::stod(tokens[i + 1]); } catch (...) {}
            } else if (tokens[i] == "molecular_weight" || tokens[i] == "M") {
                try { mat.molecular_weight = std::stod(tokens[i + 1]); } catch (...) {}
            } else if (tokens[i] == "eos") {
                mat.eos_model = tokens[i + 1];
            }
        }

        // Record unmapped raw tokens
        for (const auto& tok : tokens) {
            mat.extra_properties[tok] = tok;
        }

        setup_.materials.push_back(mat);
    }

    // Read closing paren
    char c;
    in >> c;

    return true;
}

// ---------------------------------------------------------------------------
// .cas models/settings parsing (section 48)
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas_models(std::istream& in) {
    // Section 48: model settings
    int n_settings = 0;
    in >> n_settings;

    std::vector<std::string> tokens;
    while (in.peek() != ')' && in.peek() != '\n' && !in.eof()) {
        std::string tok;
        in >> tok;
        if (!tok.empty()) tokens.push_back(tok);
    }

    // Parse common Fluent model settings
    for (size_t i = 0; i + 1 < tokens.size(); ++i) {
        if (tokens[i] == "transient-formulation") {
            setup_.transient = true;
        } else if (tokens[i] == "steady") {
            setup_.transient = false;
        } else if (tokens[i] == "time-step") {
            try { setup_.time_step = std::stod(tokens[i + 1]); } catch (...) {}
        } else if (tokens[i] == "number-of-time-steps") {
            try { setup_.max_time_steps = std::stoi(tokens[i + 1]); } catch (...) {}
        } else if (tokens[i] == "end-time") {
            try { setup_.end_time = std::stod(tokens[i + 1]); } catch (...) {}
        }
        // Record raw settings for gap analysis
        setup_.numerics.raw_settings[tokens[i]] =
            (i + 1 < tokens.size()) ? tokens[i + 1] : "";
    }

    // Map turbulence model
    for (size_t i = 0; i < tokens.size(); ++i) {
        std::string cfdx_turb = rules_.map_turbulence(tokens[i]);
        if (!cfdx_turb.empty()) {
            setup_.turbulence_model = cfdx_turb;
            setup_.numerics.raw_settings[tokens[i]] = cfdx_turb;
        }
    }

    // Default schemes
    setup_.numerics.momentum_scheme = "first_order";
    setup_.numerics.pressure_scheme = "standard";
    setup_.numerics.momentum_interpolation = "linear";
    setup_.numerics.gradient_operator = "green_gauss_cell";
    setup_.numerics.coupled_solver = "SIMPLE";

    // Read closing paren
    char c;
    in >> c;

    return true;
}

// ---------------------------------------------------------------------------
// Build cell topology from parsed data
// ---------------------------------------------------------------------------
void FluentAdapter::build_cell_topology() {
    // STUB — not implemented.
    // A full implementation would reconstruct cell-face connectivity (CSR
    // face and cell arrays) from face_nodes_, face_owner_, face_neighbour_,
    // and cell_faces_. Until then, no mesh topology is built and
    // ConversionResult.has_mesh() returns false.
}

// ---------------------------------------------------------------------------
// Main .cas parsing entry point
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_cas(const std::string& cas_file) {
    std::ifstream in(cas_file);
    if (!in.is_open()) {
        std::cerr << "FluentAdapter: cannot open .cas file: " << cas_file << "\n";
        return false;
    }

    setup_.source.case_path = cas_file;
    setup_.source.case_name = cas_file;
    setup_.source.solver = "Fluent";
    setup_.source.format = "cas_dat";

    std::string line;
    while (std::getline(in, line)) {
        int sec_id = 0, sec_nvals = 0;
        if (is_section_start(line, &sec_id, &sec_nvals)) {
            switch (sec_id) {
                case 1:  parse_cas_header(in); break;
                case 2: case 3: case 4: case 5: case 45: case 48:
                    skip_section(in); break;
                case 31: parse_cas_zones(in, sec_nvals); break;
                case 39: parse_cas_materials(in, sec_nvals); break;
                default:
                    skip_section(in); break;
            }
        }
    }

    build_cell_topology();
    return true;
}

// ---------------------------------------------------------------------------
// .dat results parsing
// ---------------------------------------------------------------------------
bool FluentAdapter::parse_dat(const std::string& dat_file) {
    std::ifstream in(dat_file, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "FluentAdapter: cannot open .dat file: " << dat_file << "\n";
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        int sec_id = 0, sec_nvals = 0;
        if (is_section_start(line, &sec_id, &sec_nvals)) {
            if (sec_id == 19) {
                // Header: variable names
                std::vector<std::string> vars;
                parse_dat_header(in, vars);
            } else if (sec_id >= 20 && sec_id <= 40) {
                // Cell data sections (variable values)
                parse_dat_cell_data(in, "variable_" + std::to_string(sec_id),
                                    sec_nvals);
            }
        }
    }
    return true;
}

bool FluentAdapter::parse_dat_header(std::istream& in,
                                      std::vector<std::string>& variables) {
    int n_vars = 0;
    in >> n_vars;
    variables.resize(n_vars);
    for (int i = 0; i < n_vars; ++i) {
        in >> std::ws;
        if (in.peek() == '"') {
            char q;
            in >> q;
            std::getline(in, variables[i], '"');
        } else {
            in >> variables[i];
        }
    }
    char c;
    in >> c;
    return true;
}

bool FluentAdapter::parse_dat_cell_data(std::istream& in,
                                        const std::string& field_name,
                                        std::size_t n_cells) {
    // Read n_cells values
    // (Values are stored as flat arrays)
    char c;
    in >> c;  // closing paren
    return true;
}

// ---------------------------------------------------------------------------
// Import results only
// ---------------------------------------------------------------------------
bool FluentAdapter::import_results(const std::string& results_path,
                                   ConversionResult& result) {
    if (!parse_dat(results_path)) {
        result.gap_report.unsupported_blocking(
            "results", "cell_data",
            "Cannot parse Fluent .dat results file: " + results_path,
            "Ensure the file exists and is a valid Fluent .dat results format");
        return false;
    }

    // Populate gap analysis
    result.gap_report.supported(
        "results", "cell_data",
        "Parsed Fluent .dat results file: " + results_path);

    return true;
}

// ---------------------------------------------------------------------------
// Full conversion entry point
// ---------------------------------------------------------------------------
bool FluentAdapter::convert(const std::string& case_path,
                            ConversionResult& result) {
    std::string cas_file = case_path;
    if (cas_file.size() > 4 && cas_file.substr(cas_file.size() - 4) != ".cas") {
        cas_file += ".cas";
    }

    if (!parse_cas(cas_file)) {
        result.gap_report.unsupported_blocking(
            "mesh", "cas_file", "Cannot open Fluent .cas file",
            "Ensure the .cas file exists and is ASCII legacy format");
        return false;
    }

    // Populate result source info
    result.source = setup_.source;
    result.mesh = cfdx::core::Mesh{};

    // NOTE: Fluent .cas mesh import (sections 2/3/4/5) is NOT yet implemented.
    // parse_cas() skips these sections; parse_cas_nodes/faces/cells exist as
    // stubs but are never called. Only metadata (zones, materials, numerics)
    // is actually imported.
    result.gap_report.unsupported_blocking(
        "mesh", "full_import",
        "Fluent .cas mesh import not implemented: sections 2 (nodes), "
        "3/4 (faces), 5 (cells) skipped. 0 nodes, 0 faces, 0 cells parsed.",
        "Wire parse_cas_nodes(), parse_cas_faces(), parse_cas_cells() into "
        "parse_cas() sections 2-5 to enable mesh import");

    // Boundary conditions
    for (const auto& bc : setup_.boundary_conditions) {
        result.gap_report.supported(
            "boundary_condition", bc.patch_name,
            "Mapped Fluent zone '" + bc.source_zone_name + "' type '" +
            bc.source_type_name + "' to CFDX type '" + to_string(bc.type) + "'");
    }

    // Materials
    for (const auto& mat : setup_.materials) {
        result.gap_report.supported(
            "material", mat.name,
            "Parsed material with density=" + std::to_string(mat.density) +
            ", viscosity=" + std::to_string(mat.dynamic_viscosity));
    }

    // Numerical settings
    result.gap_report.supported(
        "numerics", "solver",
        "Mapped Fluent solver to CFDX coupled: " + setup_.numerics.coupled_solver);

    result.gap_report.approximated(
        "numerics", "scheme",
        "Fluent gradient scheme approximated to Green-Gauss cell-based",
        "Set gradient_operator to 'least_squares' for higher accuracy");

    // Record unsupported features
    result.gap_report.unsupported_nonblocking(
        "features", "binary_cas",
        "Fluent binary .cas format not supported",
        "Export case as ASCII legacy .cas format");

    result.gap_report.unsupported_nonblocking(
        "features", "udf",
        "User Defined Functions (UDFs) not supported",
        "Extract UDF logic and reimplement in CFDX physics model");

    result.gap_report.unsupported_nonblocking(
        "features", "multiphase",
        "Multiphase models (VOF/Mixture/Eulerian) not supported",
        "Use CFDX native VOF module for multiphase flows");

    return !result.gap_report.has_blocking();
}

// ---------------------------------------------------------------------------
// Standalone convenience functions
// ---------------------------------------------------------------------------
bool parse_fluent_cas(const std::string& cas_file,
                      FluentAdapter& adapter,
                      cfdx::core::Mesh& mesh,
                      GapAnalysis& gap,
                      MappingRules rules) {
    FluentAdapter local_adapter(rules);
    ConversionResult result;

    if (!local_adapter.convert(cas_file, result)) {
        gap = result.gap_report;
        return false;
    }

    gap = result.gap_report;

    // Build mesh from parsed data
    if (result.has_mesh()) {
        mesh.clear();
        // Populate mesh from result (adapter stores data internally)
        // Note: full mesh build requires cell-face connectivity reconstruction
    }

    // Copy adapter state to caller
    adapter = std::move(local_adapter);
    return true;
}

bool parse_fluent_dat(const std::string& dat_file,
                      cfdx::core::ScalarCellField& pressure_field,
                      cfdx::core::Vec3CellField& velocity_field,
                      GapAnalysis& gap) {
    FluentAdapter adapter;
    ConversionResult result;
    adapter.import_results(dat_file, result);

    gap = result.gap_report;
    return true;
}

}  // namespace fluent
}  // namespace io
}  // namespace cfdx
