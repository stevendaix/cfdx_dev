// M0.10-T03: Code_Saturne adapter
//
// Spécification CFDX v0.7 issue #425 §Code_Saturne
// ---------------------------------------------------------------------------
// Parses Code_Saturne case files: XML setup and Python case definitions.
//
// Supported Code_Saturne formats:
//   - XML case setup (.xml) — decomposition of physics, BCs, materials
//   - Python case definition (.py) — case.set() calls
//   - MED/CGNS/EnSight/VTK result exports
//
// ARCHITECTURAL NOTE (per issue #425):
//   "Evaluate XML/Python setup extraction"
//   "Support documented neutral mesh/result exports such as MED/CGNS/EnSight/VTK"
//
// This adapter:
//   - Parses the Code_Saturne XML case file for physics model, BCs, materials
//   - Parses Python case scripts for case.set() calls (key-value extraction)
//   - Recognises when a neutral mesh export (MED/CGNS/VTK) is present
//   - Documents all limitations in the GapAnalysis
// ---------------------------------------------------------------------------
#pragma once

#include "cfdx/io/cfdx_io/io_interface.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/io/cfdx_io/mapping_rules.h"
#include <string>
#include <vector>
#include <map>

namespace cfdx {
namespace io {
namespace saturne {

// ---------------------------------------------------------------------------
// Code_Saturne setup element (XML or Python)
// ---------------------------------------------------------------------------
struct SaturneSetupEntry {
    std::string key;        // e.g. "physics", "boundary_type", "material"
    std::string value;      // e.g. "incompressible_laminar"
    std::string section;    // e.g. "model", "boundary_conditions"
    std::map<std::string, std::string> attributes;
};

// ---------------------------------------------------------------------------
// SaturneAdapter — implements SolverAdapter for Code_Saturne
// ---------------------------------------------------------------------------
class SaturneAdapter : public SolverAdapter {
public:
    SaturneAdapter() = default;
    explicit SaturneAdapter(const MappingRules& rules) : rules_(rules) {}

    const char* solver_name() const override { return "Code_Saturne"; }
    const char* format_name() const override { return "saturne_xml_py"; }

    // Full conversion: parse XML/Python setup + detect neutral exports
    bool convert(const std::string& case_path,
                 ConversionResult& result) override;

    // Results import from neutral exports (MED/CGNS/VTK/EnSight)
    bool import_results(const std::string& results_path,
                        ConversionResult& result) override;

    // Source detection
    bool detect_source(const std::string& case_path,
                       SourceInfo& info) override;

    // --- Individual parsers ---
    bool parse_xml_setup(const std::string& xml_file);
    bool parse_python_setup(const std::string& py_file);

    // --- Accessors ---
    const std::vector<SaturneSetupEntry>& entries() const { return setup_entries_; }
    const CaseSetup& setup() const { return setup_; }

private:
    // --- XML parsing helpers (minimal, no external dependency) ---
    bool parse_xml_tag(const std::string& line,
                       std::string& tag,
                       std::map<std::string, std::string>& attrs,
                       bool& is_opening,
                       bool& is_closing);
    void parse_xml_content(const std::string& tag,
                           const std::map<std::string, std::string>& attrs,
                           const std::string& content);

    // --- Python parsing helpers ---
    void parse_python_set_call(const std::string& line);
    std::string extract_python_string(const std::string& expr);
    double extract_python_number(const std::string& expr);

    // --- State ---
    MappingRules rules_ = MappingRules::load_default();
    CaseSetup setup_;
    std::vector<SaturneSetupEntry> setup_entries_;
    std::vector<std::string> neutral_files_;  // detected MED/CGNS/VTK exports
};

// Convenience API
bool parse_saturne_xml(const std::string& xml_file,
                       SaturneAdapter& adapter,
                       CaseSetup& setup,
                       GapAnalysis& gap,
                       MappingRules rules = MappingRules::load_default());

bool parse_saturne_python(const std::string& py_file,
                          SaturneAdapter& adapter,
                          CaseSetup& setup,
                          GapAnalysis& gap);

}  // namespace saturne
}  // namespace io
}  // namespace cfdx
