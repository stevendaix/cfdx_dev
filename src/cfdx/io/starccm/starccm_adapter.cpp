// M0.10-T03: STAR-CCM+ .sim adapter implementation
#include "starccm_adapter.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>

namespace cfdx {
namespace io {
namespace starccm {

// ---------------------------------------------------------------------------
// Read .sim file header
// ---------------------------------------------------------------------------
bool StarCCMAdapter::read_header(const std::string& sim_path) {
    std::ifstream in(sim_path, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "StarCCMAdapter: cannot open .sim file: " << sim_path << "\n";
        return false;
    }

    // Get file size
    in.seekg(0, std::ios::end);
    header_.file_size = static_cast<std::size_t>(in.tellg());
    in.seekg(0, std::ios::beg);

    if (header_.file_size < 16) {
        header_.magic = "unknown";
        header_.version = "unknown";
        return false;
    }

    // Read first 256 bytes to find header patterns
    std::vector<char> buffer(std::min(header_.file_size, static_cast<std::size_t>(256)));
    in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));

    // Detect endianness: STAR-CCM+ uses little-endian on x86
    header_.is_little_endian = true;

    // Convert buffer to string for ASCII pattern matching
    std::string header_str(buffer.begin(), buffer.end());

    // Look for version string
    size_t ver_pos = header_str.find("Version");
    if (ver_pos != std::string::npos) {
        size_t start = ver_pos + 7;
        size_t end = header_str.find_first_of("\n\r\0", start);
        if (end != std::string::npos) {
            header_.version = header_str.substr(start, end - start);
        }
    }

    // Look for STAR-CCM+ magic
    if (header_str.find("STAR-CCM") != std::string::npos) {
        header_.magic = "STAR-CCM+";
    } else {
        // Check hex magic number for STAR-CCM+ .sim files
        // The file may start with a 4-byte record marker
        unsigned char magic_bytes[4];
        std::memcpy(magic_bytes, buffer.data(), 4);
        // STAR-CCM+ .sim files often start with specific binary markers
        if (magic_bytes[0] == 0x00 || magic_bytes[0] == 0x01) {
            header_.magic = "binary (STAR-CCM+ proprietary)";
            header_.format_note = "Binary format detected; direct topology parsing not supported";
        } else {
            header_.magic = "unknown";
        }
    }

    // Build number from build string if present
    size_t build_pos = header_str.find("Build");
    if (build_pos != std::string::npos) {
        size_t start = build_pos + 5;
        size_t end = header_str.find_first_of("\n\r\0", start);
        if (end != std::string::npos) {
            header_.build = header_str.substr(start, end - start);
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// Scan entire .sim file for ASCII strings (zone names, materials, etc.)
// ---------------------------------------------------------------------------
bool StarCCMAdapter::scan_strings(const std::string& sim_path) {
    std::ifstream in(sim_path, std::ios::binary);
    if (!in.is_open()) return false;

    // Read file in chunks for memory efficiency
    // For large .sim files, we scan the raw bytes for ASCII strings
    // of minimum length 4 that look like zone names or material names
    in.seekg(0, std::ios::end);
    std::size_t fsize = static_cast<std::size_t>(in.tellg());
    in.seekg(0, std::ios::beg);

    // Read entire file (for moderate sizes; chunk-based for large files)
    std::vector<char> buffer(fsize);
    in.read(buffer.data(), static_cast<std::streamsize>(fsize));

    scan_ascii_metadata(buffer);
    return true;
}

// ---------------------------------------------------------------------------
// Scan ASCII metadata from raw bytes buffer
// ---------------------------------------------------------------------------
bool StarCCMAdapter::scan_ascii_metadata(const std::vector<char>& buffer) {
    // Extract all ASCII strings >= 4 characters
    std::vector<std::string> strings;
    std::string current;

    for (char c : buffer) {
        if (c >= 32 && c <= 126) {
            current += c;
        } else {
            if (current.size() >= 4) {
                strings.push_back(current);
            }
            current.clear();
        }
    }
    if (current.size() >= 4) strings.push_back(current);

    // Look for patterns indicating boundary zones
    // STAR-CCM+ zones often appear as: "Region/ Boundaries/<name>" or "Surface"
    for (size_t i = 0; i < strings.size(); ++i) {
        const std::string& s = strings[i];

        // Boundary zone detection
        if (s.find("Boundaries") != std::string::npos ||
            s.find("Boundary") != std::string::npos) {
            // Look at surrounding strings for zone names
            for (size_t j = std::max(size_t(0), i - 2);
                 j < std::min(strings.size(), i + 4); ++j) {
                if (j != i && strings[j].size() >= 4) {
                    SimZone zone;
                    zone.name = strings[j];
                    zone.geometry_type = "surface";
                    classify_boundary_type(zone);

                    // Check if zone already seen
                    bool found = false;
                    for (const auto& z : zones_) {
                        if (z.name == zone.name) { found = true; break; }
                    }
                    if (!found && zone.name != s) {
                        zones_.push_back(zone);
                    }
                }
            }
        }

        // Material detection
        if (s.find("Material") != std::string::npos ||
            s.find("material") != std::string::npos) {
            for (size_t j = std::max(size_t(0), i - 1);
                 j < std::min(strings.size(), i + 2); ++j) {
                if (j != i && strings[j].size() >= 4 &&
                    strings[j] != "Material" && strings[j] != "material") {
                    SimMaterial mat;
                    mat.name = strings[j];
                    mat.model = "constant";
                    // Avoid duplicates
                    bool found = false;
                    for (const auto& m : materials_) {
                        if (m.name == mat.name) { found = true; break; }
                    }
                    if (!found) materials_.push_back(mat);
                }
            }
        }

        // Physics model detection
        if (s.find("Turbulence") != std::string::npos) {
            for (const auto& s2 : strings) {
                if (s2 == "k-epsilon" || s2 == "k-omega SST" || s2 == "Spalart-Allmaras" ||
                    s2 == "Realizable k-epsilon" || s2 == "Mixing Length") {
                    physics_.turbulence = s2;
                    break;
                }
            }
        }
        if (s.find("Energy") != std::string::npos && s.find("energy") == std::string::npos) {
            physics_.energy = "on";
        }
        if (s.find("VOF") != std::string::npos) {
            physics_.multiphase = "VOF";
        }
    }

    // Detect physics regime
    for (const auto& s : strings) {
        if (s.find("Incompressible") != std::string::npos) {
            physics_.regime = "incompressible";
            break;
        }
        if (s.find("Compressible") != std::string::npos) {
            physics_.regime = "compressible";
            break;
        }
    }

    // Look for cell count near "Cells" keyword
    for (size_t i = 0; i < strings.size(); ++i) {
        if (strings[i] == "Cells" || strings[i] == "cells") {
            // Check next few strings for a number
            for (size_t j = i + 1; j < std::min(strings.size(), i + 5); ++j) {
                try {
                    std::size_t val = std::stoull(strings[j]);
                    if (val > 0 && val < 10000000000) {
                        setup_.mesh_info.n_cells = val;
                        break;
                    }
                } catch (...) {}
            }
            break;
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// Classify boundary type from zone name
// ---------------------------------------------------------------------------
void StarCCMAdapter::classify_boundary_type(SimZone& zone) {
    std::string lower = zone.name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    zone.boundary_type = "unknown";

    // Use mapping rules to classify
    std::string zone_name_lower = zone.name;
    // Try mapping rules
    std::string mapped = rules_.map_turbulence(zone_name_lower);
    if (mapped.empty()) {
        // Direct name matching
        if (lower.find("inlet") != std::string::npos ||
            lower.find("velocity") != std::string::npos) {
            zone.boundary_type = "velocity-inlet";
        } else if (lower.find("outlet") != std::string::npos ||
                   lower.find("pressure") != std::string::npos) {
            zone.boundary_type = "pressure-outlet";
        } else if (lower.find("wall") != std::string::npos) {
            zone.boundary_type = "wall";
        } else if (lower.find("symmetr") != std::string::npos) {
            zone.boundary_type = "symmetry";
        } else if (lower.find("period") != std::string::npos) {
            zone.boundary_type = "periodic";
        } else if (lower.find("interface") != std::string::npos) {
            zone.boundary_type = "interface";
        }
    } else {
        zone.boundary_type = mapped;
    }
}

// ---------------------------------------------------------------------------
// Source detection
// ---------------------------------------------------------------------------
bool StarCCMAdapter::detect_source(const std::string& sim_path,
                                   SourceInfo& info) {
    // Check if .sim file exists and has STAR-CCM+ markers
    std::ifstream test(sim_path, std::ios::binary);
    if (!test.is_open()) return false;

    // Read first 1024 bytes
    std::vector<char> buf(1024);
    test.read(buf.data(), static_cast<std::streamsize>(
        std::min(static_cast<std::size_t>(buf.size()),
                 static_cast<std::size_t>(test.tellg()))));

    std::string header_str(buf.begin(), buf.end());
    bool is_star_ccm = header_str.find("STAR-CCM") != std::string::npos ||
                       header_str.find("STAR CCM") != std::string::npos ||
                       header_str.find("cdl") != std::string::npos;  // common in .sim

    if (!is_star_ccm) {
        // Still try — some .sim files have binary headers
        test.seekg(0, std::ios::end);
        std::size_t sz = static_cast<std::size_t>(test.tellg());
        if (sz > 1000) {  // .sim files are typically large
            header_.magic = "binary (STAR-CCM+ proprietary)";
        } else {
            return false;
        }
    }

    info.solver = "STAR-CCM+";
    info.case_path = sim_path;
    info.case_name = sim_path;
    info.format = "sim";
    info.version = header_.version.empty() ? "unknown" : header_.version;
    return true;
}

// ---------------------------------------------------------------------------
// Full conversion entry point
// ---------------------------------------------------------------------------
bool StarCCMAdapter::convert(const std::string& sim_path,
                             ConversionResult& result) {
    if (!read_header(sim_path)) {
        result.gap_report.unsupported_blocking(
            "mesh", "sim_file",
            "Cannot read STAR-CCM+ .sim file header",
            "Ensure the file is a valid .sim file");
        return false;
    }

    result.source.solver = "STAR-CCM+";
    result.source.version = header_.version;
    result.source.case_path = sim_path;
    result.source.case_name = sim_path;
    result.source.format = "sim";

    // Scan for ASCII metadata
    scan_strings(sim_path);

    // Populate setup from extracted metadata
    setup_.source = result.source;
    setup_.physics_model = "incompressible_laminar";
    setup_.turbulence_model = physics_.turbulence.empty() ? "laminar" : physics_.turbulence;
    setup_.energy_model = physics_.energy == "on" ? "buoyant" : "isothermal";
    setup_.multiphase_model = physics_.multiphase.empty() ? "none" : physics_.multiphase;
    setup_.transient = false;
    setup_.numerics.coupled_solver = "SIMPLE";
    setup_.numerics.gradient_operator = "green_gauss_cell";

    // Populate boundary conditions
    for (const auto& zone : zones_) {
        BoundarySpec bc;
        bc.patch_name = zone.name;
        bc.type = rules_.map_boundary_type(zone.boundary_type);
        bc.value_type = rules_.map_bc_value_type(zone.boundary_type);
        bc.source_zone_name = zone.name;
        bc.source_type_name = zone.boundary_type;
        bc.raw_params["n_faces"] = std::to_string(zone.n_faces);
        bc.raw_params["geometry_type"] = zone.geometry_type;

        setup_.boundary_conditions.push_back(bc);

        result.gap_report.supported(
            "boundary_condition", zone.name,
            "Detected zone '" + zone.name + "' with boundary type '" +
            zone.boundary_type + "' (string-scan based, not topology-verified)");
    }

    // Populate materials
    for (const auto& mat : materials_) {
        MaterialSpec m;
        m.name = mat.name;
        m.eos_model = "ideal_gas";
        setup_.materials.push_back(m);

        result.gap_report.supported(
            "material", mat.name,
            "Detected material '" + mat.name + "'");
    }

    // Populate mesh metadata
    if (setup_.mesh_info.n_cells > 0) {
        result.gap_report.supported(
            "mesh", "cell_count",
            "Detected " + std::to_string(setup_.mesh_info.n_cells) +
            " cells from .sim header/metadata (binary topology not parsed)");
    }

    // Record unsupported features — DOCUMENTED, NOT SILENT
    result.gap_report.unsupported_nonblocking(
        "mesh", "binary_topology",
        "STAR-CCM+ .sim uses proprietary binary format; mesh topology "
        "(faces, cells, owner/neighbour) cannot be directly reconstructed",
        "Export mesh as .msh (Gmsh) or .su2 via File > Export > Mesh");

    result.gap_report.unsupported_blocking(
        "results", "sim_results",
        "Field results cannot be imported directly from .sim binary database",
        "Export results via: File > Export > Results File (.enc, .plt, .vtk, .cgns)");

    result.gap_report.unsupported_nonblocking(
        "features", "udf",
        "STAR-CCM+ Java macros and UDFs are not parseable",
        "Extract macro logic and reimplement in CFDX Python orchestration");

    result.gap_report.unsupported_nonblocking(
        "features", "custom_models",
        "STAR-CCM+ custom physics models are not translatable",
        "Map to nearest CFDX physics module and document deviation");

    // Record physics detection
    if (!physics_.turbulence.empty()) {
        result.gap_report.supported(
            "physics", "turbulence",
            "Detected turbulence model: " + physics_.turbulence);
    }

    result.gap_report.approximated(
        "conversion", "binary_format",
        "Only ASCII metadata strings were extracted from .sim binary. "
        "Mesh topology, connectivity, and fields were NOT reconstructed.",
        "Use documented STAR-CCM+ export (.su2, .msh, .cgns) for full conversion");

    return !result.gap_report.has_blocking();
}

// ---------------------------------------------------------------------------
// Standalone convenience function
// ---------------------------------------------------------------------------
bool parse_starccm_sim(const std::string& sim_path,
                       StarCCMAdapter& adapter,
                       GapAnalysis& gap,
                       MappingRules rules) {
    StarCCMAdapter local_adapter(rules);
    ConversionResult result;

    local_adapter.detect_source(sim_path, result.source);

    if (!local_adapter.convert(sim_path, result)) {
        gap = result.gap_report;
        adapter = std::move(local_adapter);
        return false;
    }

    gap = result.gap_report;
    adapter = std::move(local_adapter);
    return true;
}

}  // namespace starccm
}  // namespace io
}  // namespace cfdx
