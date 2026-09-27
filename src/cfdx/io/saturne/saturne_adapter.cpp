// M0.10-T03: Code_Saturne adapter implementation
#include "saturne_adapter.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <regex>
#include <filesystem>

namespace cfdx {
namespace io {
namespace saturne {

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Helper: trim string
// ---------------------------------------------------------------------------
static std::string trim_sat(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static std::string to_lower_sat(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

// ---------------------------------------------------------------------------
// Source detection
// ---------------------------------------------------------------------------
bool SaturneAdapter::detect_source(const std::string& case_path,
                                   SourceInfo& info) {
    // Code_Saturne cases have a specific directory structure:
    //   <case>.xml        — setup file
    //   <case>.py         — Python case definition
    //   SRC/              — mesh sources
    //   POST/             — results
    //   RESU/             — results directory

    // Look for XML setup file
    std::string xml_file = case_path;
    if (xml_file.size() > 4 && xml_file.substr(xml_file.size() - 4) == ".xml") {
        // Direct XML file
    } else {
        // Search for .xml in directory
        fs::path case_dir(case_path);
        if (fs::is_directory(case_dir)) {
            for (const auto& entry : fs::directory_iterator(case_dir)) {
                if (entry.path().extension() == ".xml") {
                    xml_file = entry.path().string();
                    break;
                }
            }
        }
    }

    std::ifstream test(xml_file);
    if (!test.is_open()) {
        // Try Python file
        std::string py_file = case_path;
        if (py_file.size() > 3 && py_file.substr(py_file.size() - 3) != ".py") {
            py_file += ".py";
        }
        std::ifstream test_py(py_file);
        if (!test_py.is_open()) return false;

        info.solver = "Code_Saturne";
        info.case_path = case_path;
        info.case_name = py_file;
        info.format = "saturne_py";
        info.version = "unknown";
        return true;
    }

    // Verify XML content
    std::string first_check;
    std::getline(test, first_check);
    if (first_check.find("<?xml") != std::string::npos ||
        first_check.find("mesh") != std::string::npos ||
        first_check.find("case") != std::string::npos) {
        info.solver = "Code_Saturne";
        info.case_path = case_path;
        info.case_name = xml_file;
        info.format = "saturne_xml";
        info.version = "unknown";
        return true;
    }

    return false;
}

// ---------------------------------------------------------------------------
// XML tag parsing (minimal, no external XML library)
// ---------------------------------------------------------------------------
bool SaturneAdapter::parse_xml_tag(const std::string& line,
                                   std::string& tag,
                                   std::map<std::string, std::string>& attrs,
                                   bool& is_opening,
                                   bool& is_closing) {
    std::string s = trim_sat(line);
    if (s.empty() || s[0] != '<') return false;

    is_opening = true;
    is_closing = false;

    // Check for closing tag
    if (s.size() > 1 && s[1] == '/') {
        is_closing = true;
        is_opening = false;
        s = s.substr(2);
    } else if (s.size() > 1 && s[1] != '?') {
        s = s.substr(1);
    } else {
        return false;  // skip comments/xml declarations
    }

    // Extract tag name
    std::istringstream iss(s);
    iss >> tag;

    // Parse attributes: key="value" pairs
    std::string token;
    while (iss >> token) {
        // Check for self-closing tag
        if (token.find("/>") != std::string::npos) {
            is_closing = true;
            token = token.substr(0, token.size() - 2);
        }

        size_t eq = token.find('=');
        if (eq != std::string::npos) {
            std::string key = token.substr(0, eq);
            std::string val = token.substr(eq + 1);

            // Strip quotes
            if (!val.empty() && val.front() == '"') val.erase(0, 1);
            if (!val.empty() && val.back() == '"') val.pop_back();

            attrs[key] = val;
        }
    }

    // Remove trailing '>' from tag if present
    size_t gt = tag.find('>');
    if (gt != std::string::npos) tag = tag.substr(0, gt);

    return !tag.empty();
}

// ---------------------------------------------------------------------------
// Process XML content based on tag
// ---------------------------------------------------------------------------
void SaturneAdapter::parse_xml_content(const std::string& tag,
                                       const std::map<std::string, std::string>& attrs,
                                       const std::string& content) {
    SaturneSetupEntry entry;
    entry.section = tag;

    // Map XML tags to setup entries
    if (tag == "model") {
        entry.key = "physics";
        entry.value = attrs.count("type") ? attrs.at("type") : "unknown";
    } else if (tag == "boundary") {
        entry.key = "boundary_type";
        entry.value = attrs.count("type") ? attrs.at("type") : "unknown";
        if (attrs.count("name")) entry.attributes["name"] = attrs.at("name");
        if (attrs.count("marker")) entry.attributes["marker"] = attrs.at("marker");
    } else if (tag == "material") {
        entry.key = "material";
        entry.value = attrs.count("name") ? attrs.at("name") : "unknown";
    } else if (tag == "turbulence") {
        entry.key = "turbulence_model";
        entry.value = attrs.count("model") ? attrs.at("model") : "unknown";
    } else if (tag == "mesh") {
        entry.key = "mesh";
        if (attrs.count("file")) entry.value = attrs.at("file");
    } else if (tag == "results") {
        entry.key = "results";
        if (attrs.count("format")) entry.attributes["format"] = attrs.at("format");
        if (attrs.count("file")) entry.attributes["file"] = attrs.at("file");
    }

    // Store all attributes as raw parameters
    for (const auto& [k, v] : attrs) {
        if (k != "type" && k != "name" && k != "model") {
            entry.attributes[k] = v;
        }
    }

    if (!entry.key.empty()) {
        setup_entries_.push_back(entry);
    }

    // Map to CFDX schema
    if (tag == "model" && attrs.count("type")) {
        std::string model = to_lower_sat(attrs.at("type"));
        if (model.find("incompressible") != std::string::npos) {
            setup_.physics_model = "incompressible_laminar";
        } else if (model.find("compressible") != std::string::npos) {
            setup_.physics_model = "compressible";
        }
    }

    if (tag == "turbulence" && attrs.count("model")) {
        std::string turb = attrs.at("model");
        std::string mapped = rules_.map_turbulence(turb);
        setup_.turbulence_model = mapped.empty() ? to_lower_sat(turb) : mapped;
    }

    if (tag == "material" && attrs.count("name")) {
        MaterialSpec mat;
        mat.name = attrs.at("name");
        if (attrs.count("density")) {
            try { mat.density = std::stod(attrs.at("density")); } catch (...) {}
        }
        if (attrs.count("viscosity")) {
            try { mat.dynamic_viscosity = std::stod(attrs.at("viscosity")); } catch (...) {}
        }
        if (attrs.count("conductivity")) {
            try { mat.thermal_conductivity = std::stod(attrs.at("conductivity")); } catch (...) {}
        }
        if (attrs.count("specific_heat")) {
            try { mat.specific_heat = std::stod(attrs.at("specific_heat")); } catch (...) {}
        }
        setup_.materials.push_back(mat);
    }

    if (tag == "boundary" && attrs.count("type")) {
        BoundarySpec bc;
        bc.source_zone_name = attrs.count("name") ? attrs.at("name") : "";
        bc.source_type_name = attrs.at("type");
        std::string mapped_type = to_lower_sat(attrs.at("type"));
        bc.type = rules_.map_boundary_type(mapped_type);
        bc.value_type = rules_.map_bc_value_type(mapped_type);

        if (attrs.count("name")) bc.patch_name = attrs.at("name");
        if (attrs.count("marker")) bc.raw_params["marker"] = attrs.at("marker");

        setup_.boundary_conditions.push_back(bc);
    }

    if (tag == "results" && attrs.count("format")) {
        neutral_files_.push_back(attrs.at("format"));
    }
}

// ---------------------------------------------------------------------------
// Parse XML setup file
// ---------------------------------------------------------------------------
bool SaturneAdapter::parse_xml_setup(const std::string& xml_file) {
    std::ifstream in(xml_file);
    if (!in.is_open()) {
        std::cerr << "SaturneAdapter: cannot open XML file: " << xml_file << "\n";
        return false;
    }

    setup_.source.solver = "Code_Saturne";
    setup_.source.case_path = xml_file;
    setup_.source.case_name = xml_file;
    setup_.source.format = "saturne_xml";

    std::string line;
    std::string current_tag;
    std::map<std::string, std::string> current_attrs;
    bool tag_open = false;

    while (std::getline(in, line)) {
        // Skip XML declaration
        if (line.find("<?xml") != std::string::npos) continue;

        std::string tag;
        std::map<std::string, std::string> attrs;
        bool is_opening = false;
        bool is_closing = false;

        if (parse_xml_tag(line, tag, attrs, is_opening, is_closing)) {
            if (is_opening && !is_closing) {
                current_tag = tag;
                current_attrs = attrs;
                tag_open = true;
            } else if (is_closing) {
                // Self-closing tag — process immediately
                if (tag_open) {
                    parse_xml_content(tag, attrs, "");
                    tag_open = false;
                }
                // Also handle standalone self-closing tags
                if (line.find("/>") != std::string::npos && line.find("</") == std::string::npos) {
                    parse_xml_content(tag, attrs, "");
                }
            }
        }
    }

    return !setup_entries_.empty();
}

// ---------------------------------------------------------------------------
// Python string extraction
// ---------------------------------------------------------------------------
std::string SaturneAdapter::extract_python_string(const std::string& expr) {
    std::string v = trim_sat(expr);
    if (v.empty()) return "";
    if (v.front() == '\'' && v.back() == '\'') return v.substr(1, v.size() - 2);
    if (v.front() == '"' && v.back() == '"') return v.substr(1, v.size() - 2);
    return v;
}

double SaturneAdapter::extract_python_number(const std::string& expr) {
    try {
        return std::stod(trim_sat(expr));
    } catch (...) {
        return 0.0;
    }
}

// ---------------------------------------------------------------------------
// Parse Python case.set() call
// ---------------------------------------------------------------------------
void SaturneAdapter::parse_python_set_call(const std::string& line) {
    // Match patterns like: case.set('key', 'value') or case.set('key', value)
    std::regex pattern(R"(set\s*\(\s*['"]([^'"]+)['"]\s*,\s*['"]?([^'")]+)['"]?\s*\))");
    std::smatch matches;

    if (std::regex_search(line, matches, pattern)) {
        if (matches.size() >= 3) {
            SaturneSetupEntry entry;
            entry.key = trim_sat(matches[1].str());
            entry.value = trim_sat(matches[2].str());

            // Classify by key
            std::string lower_key = to_lower_sat(entry.key);

            if (lower_key.find("physic") != std::string::npos) {
                entry.section = "physics";
                setup_.physics_model = entry.value;
            } else if (lower_key.find("turbul") != std::string::npos) {
                entry.section = "turbulence";
                std::string mapped = rules_.map_turbulence(entry.value);
                setup_.turbulence_model = mapped.empty() ? entry.value : mapped;
            } else if (lower_key.find("mater") != std::string::npos ||
                       lower_key.find("fluid") != std::string::npos) {
                entry.section = "material";
                if (!setup_.materials.empty()) {
                    setup_.materials.back().name = entry.value;
                } else {
                    MaterialSpec mat;
                    mat.name = entry.value;
                    setup_.materials.push_back(mat);
                }
            } else if (lower_key.find("bound") != std::string::npos) {
                entry.section = "boundary";
                // Parse boundary type
                BoundarySpec bc;
                std::string mapped = to_lower_sat(entry.value);
                bc.type = rules_.map_boundary_type(mapped);
                bc.value_type = rules_.map_bc_value_type(mapped);
                bc.source_type_name = entry.value;
                bc.patch_name = entry.value;
                setup_.boundary_conditions.push_back(bc);
            } else if (lower_key.find("time") != std::string::npos) {
                entry.section = "time";
                try { setup_.time_step = std::stod(entry.value); } catch (...) {}
            } else if (lower_key.find("energy") != std::string::npos) {
                entry.section = "energy";
                setup_.energy_model = (to_lower_sat(entry.value).find("yes") != std::string::npos)
                    ? "buoyant" : "isothermal";
            }

            setup_entries_.push_back(entry);
        }
    }
}

// ---------------------------------------------------------------------------
// Parse Python case definition
// ---------------------------------------------------------------------------
bool SaturneAdapter::parse_python_setup(const std::string& py_file) {
    std::ifstream in(py_file);
    if (!in.is_open()) {
        std::cerr << "SaturneAdapter: cannot open Python file: " << py_file << "\n";
        return false;
    }

    setup_.source.solver = "Code_Saturne";
    setup_.source.case_path = py_file;
    setup_.source.case_name = py_file;
    setup_.source.format = "saturne_py";

    std::string line;
    while (std::getline(in, line)) {
        std::string trimmed = trim_sat(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;

        if (trimmed.find("case.set(") != std::string::npos ||
            trimmed.find("cs_.set(") != std::string::npos ||
            trimmed.find(".set(") != std::string::npos) {
            parse_python_set_call(trimmed);
        }
    }

    // Default to at least laminar if no turbulence found
    if (setup_.turbulence_model.empty()) {
        setup_.turbulence_model = "laminar";
    }

    return !setup_entries_.empty();
}

// ---------------------------------------------------------------------------
// Detect neutral result exports in a case directory
// ---------------------------------------------------------------------------
static void detect_neutral_files(const std::string& case_dir,
                                 std::vector<std::string>& neutral_files) {
    fs::path dir(case_dir);
    if (!fs::is_directory(dir)) return;

    for (const auto& entry : fs::recursive_directory_iterator(dir)) {
        std::string ext = entry.path().extension();
        std::string lower_ext = to_lower_sat(ext);

        if (lower_ext == ".med" || lower_ext == ".cgns" ||
            lower_ext == ".vtk" || lower_ext == ".vtu" ||
            lower_ext == ".case" || lower_ext == ".geo") {
            neutral_files.push_back(entry.path().string());
        }
    }
}

// ---------------------------------------------------------------------------
// Import results from neutral export
// ---------------------------------------------------------------------------
bool SaturneAdapter::import_results(const std::string& results_path,
                                    ConversionResult& result) {
    std::string lower_path = to_lower_sat(results_path);

    if (lower_path.find(".med") != std::string::npos) {
        result.gap_report.supported(
            "results", "med_import",
            "Code_Saturne MED result export detected: " + results_path);
        result.gap_report.unsupported_nonblocking(
            "results", "med_parser",
            "Direct MED parser not implemented in this adapter",
            "Import MED via meshio adapter or HDF5 reader");
    } else if (lower_path.find(".cgns") != std::string::npos) {
        result.gap_report.supported(
            "results", "cgns_import",
            "Code_Saturne CGNS result export detected: " + results_path);
        result.gap_report.unsupported_nonblocking(
            "results", "cgns_parser",
            "Direct CGNS parser not implemented in this adapter",
            "Import CGNS via HDF5 reader or meshio adapter");
    } else if (lower_path.find(".vtk") != std::string::npos ||
               lower_path.find(".vtu") != std::string::npos) {
        result.gap_report.supported(
            "results", "vtk_import",
            "VTK/VTU result export detected: " + results_path);
        result.gap_report.unsupported_nonblocking(
            "results", "vtk_parser",
            "VTK parser delegates to meshio adapter",
            "Use meshio adapter for field extraction from VTK");
    } else {
        result.gap_report.unsupported_nonblocking(
            "results", "unknown_format",
            "Unrecognised results format: " + results_path,
            "Supported: .med, .cgns, .vtk, .vtu, .case");
    }

    return true;
}

// ---------------------------------------------------------------------------
// Full conversion entry point
// ---------------------------------------------------------------------------
bool SaturneAdapter::convert(const std::string& case_path,
                             ConversionResult& result) {
    // Detect source
    SourceInfo info;
    if (!detect_source(case_path, info)) {
        result.gap_report.unsupported_blocking(
            "source", "saturne_case",
            "Cannot identify Code_Saturne case at: " + case_path,
            "Ensure the path contains an XML or Python case definition");
        return false;
    }

    result.source = info;

    bool parsed = false;

    // Try XML first
    if (info.format == "saturne_xml" || info.format == "saturne_py") {
        parsed = parse_xml_setup(info.case_path);
        if (!parsed) {
            // Try Python
            std::string py_file = case_path;
            if (py_file.size() > 3 && py_file.substr(py_file.size() - 3) != ".py") {
                py_file += ".py";
            }
            parsed = parse_python_setup(py_file);
        }
    } else {
        // Directory path: look for XML or Python
        std::string xml_candidate = case_path;
        if (xml_candidate.back() != '/') xml_candidate += "/";
        // Try to find .xml
        fs::path dir(case_path);
        if (fs::is_directory(dir)) {
            for (const auto& entry : fs::directory_iterator(dir)) {
                if (entry.path().extension() == ".xml") {
                    parsed = parse_xml_setup(entry.path().string());
                    break;
                }
                if (entry.path().extension() == ".py") {
                    parsed = parse_python_setup(entry.path().string());
                    break;
                }
            }
        }
    }

    if (!parsed) {
        result.gap_report.unsupported_nonblocking(
            "setup", "parse",
            "Could not parse Code_Saturne setup (XML or Python)",
            "Ensure the case has a valid .xml or .py setup file");
    }

    // Copy setup to result context
    // Result source already set

    // --- Gap Analysis ---

    result.gap_report.supported(
        "setup", "saturne_xml_py",
        "Parsed Code_Saturne setup from " + info.format + " format");

    // Boundary conditions
    for (const auto& bc : setup_.boundary_conditions) {
        result.gap_report.supported(
            "boundary_condition", bc.patch_name,
            "Mapped Code_Saturne BC: " + bc.source_type_name +
            " → CFDX type " + to_string(bc.type));
    }

    // Materials
    for (const auto& mat : setup_.materials) {
        result.gap_report.supported(
            "material", mat.name,
            "Parsed material: " + mat.name +
            " (density=" + std::to_string(mat.density) + ")");
    }

    // Physics
    if (!setup_.turbulence_model.empty() && setup_.turbulence_model != "laminar") {
        result.gap_report.supported(
            "physics", "turbulence",
            "Mapped turbulence model: " + setup_.turbulence_model);
    } else {
        result.gap_report.supported(
            "physics", "turbulence", "Laminar model (no turbulence)");
    }

    if (!setup_.energy_model.empty() && setup_.energy_model != "isothermal") {
        result.gap_report.supported(
            "physics", "energy",
            "Energy model: " + setup_.energy_model);
    }

    // Check for neutral mesh/result exports
    detect_neutral_files(case_path, neutral_files_);
    for (const auto& f : neutral_files_) {
        result.gap_report.supported(
            "mesh", "neutral_export",
            "Detected neutral export: " + f);
    }

    // Document limitations
    result.gap_report.unsupported_nonblocking(
        "mesh", "direct_topology",
        "Code_Saturne does not store mesh topology in a single parseable file; "
        "mesh import requires neutral export",
        "Export mesh as .msh, .su2, or .cgns from the Code_Saturne GUI");

    result.gap_report.unsupported_nonblocking(
        "features", "user_subroutines",
        "Code_Saturne user subroutines (.c files in SRC/) are not translatable",
        "Extract subroutine logic and reimplement in CFDX physics model");

    result.gap_report.unsupported_nonblocking(
        "features", "coupling",
        "Code_Saturne coupling options (MPI, multi-code coupling) not supported",
        "Use single-code simulation for CFDX conversion");

    // Check for unmapped raw parameters
    for (const auto& entry : setup_entries_) {
        if (entry.section == "unknown" || entry.section.empty()) {
            result.gap_report.unsupported_nonblocking(
                "config", entry.key,
                "Unmapped Code_Saturne parameter: " + entry.key + " = " + entry.value,
                "Verify the parameter is applied in CFDX setup");
        }
    }

    return !result.gap_report.has_blocking();
}

// ---------------------------------------------------------------------------
// Standalone convenience functions
// ---------------------------------------------------------------------------
bool parse_saturne_xml(const std::string& xml_file,
                       SaturneAdapter& adapter,
                       CaseSetup& setup,
                       GapAnalysis& gap,
                       MappingRules rules) {
    SaturneAdapter local_adapter(rules);
    local_adapter.detect_source(xml_file, setup.source);

    if (!local_adapter.parse_xml_setup(xml_file)) {
        gap.unsupported_blocking(
            "setup", "saturne_xml",
            "Cannot parse Code_Saturne XML: " + xml_file,
            "Ensure the file is a valid Code_Saturne XML case definition");
        adapter = std::move(local_adapter);
        return false;
    }

    setup = local_adapter.setup();
    gap = GapAnalysis{};
    gap.supported("setup", "xml_parse", "Parsed Code_Saturne XML: " + xml_file);
    adapter = std::move(local_adapter);
    return true;
}

bool parse_saturne_python(const std::string& py_file,
                          SaturneAdapter& adapter,
                          CaseSetup& setup,
                          GapAnalysis& gap) {
    SaturneAdapter local_adapter;

    if (!local_adapter.parse_python_setup(py_file)) {
        gap.unsupported_nonblocking(
            "setup", "saturne_py",
            "Could not extract any case.set() calls from Python file: " + py_file,
            "Ensure the Python file contains case.set() calls");
        adapter = std::move(local_adapter);
        return false;
    }

    setup = local_adapter.setup();
    gap = GapAnalysis{};
    gap.supported("setup", "python_parse", "Parsed Code_Saturne Python: " + py_file);
    adapter = std::move(local_adapter);
    return true;
}

}  // namespace saturne
}  // namespace io
}  // namespace cfdx
