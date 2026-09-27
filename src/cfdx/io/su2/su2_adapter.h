// M0.10-T03: SU2 adapter
//
// Spécification CFDX v0.7 issue #425 §SU2
// ---------------------------------------------------------------------------
// Parses SU2 case files: .su2 mesh and .cfg configuration.
//
// Supported SU2 mesh format (.su2):
//   - POINTS section (vertex coordinates)
//   - ELEMENTS section (element types: tetra, hexa, etc.)
//   - BOUNDARY sections (patch names, face lists)
//
// Supported SU2 config format (.cfg):
//   - Physical parameters (Mach, AoA, etc.)
//   - Boundary condition markers
//   - Numerical methods (MUSCL, scheme type)
//   - Solver settings
//
// Supported SU2 solution format:
//   - solution.csv (scalar/vector fields per vertex)
//   - history.csv (convergence history)
//   - VTK output (through meshio adapter)
// ---------------------------------------------------------------------------
#pragma once

#include "cfdx/io/cfdx_io/io_interface.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/io/cfdx_io/mapping_rules.h"
#include <string>
#include <vector>
#include <cstdint>

namespace cfdx {
namespace io {
namespace su2 {

// ---------------------------------------------------------------------------
// SU2 element type codes (from SU2 format specification)
// ---------------------------------------------------------------------------
enum class Su2ElementType : int {
    LINE          = 3,
    TRIANGLE      = 5,
    QUADRILATERAL = 6,
    TETRAHEDRON   = 9,
    HEXAHEDRON    = 10,
    PRISM         = 11,
    PYRAMID       = 12,
    POINT         = 1,
    BAR           = 3,
};

// ---------------------------------------------------------------------------
// SU2 boundary marker
// ---------------------------------------------------------------------------
struct Su2Boundary {
    std::string name;
    std::string element_type;  // "triangle", "quad"
    std::size_t n_elements = 0;
    std::vector<std::uint32_t> node_ids;
    bool is_internal() const {
        return name.find("INTERNAL") != std::string::npos ||
               name == "inner_faces";
    }
};

// ---------------------------------------------------------------------------
// Su2Adapter — implements SolverAdapter for SU2
// ---------------------------------------------------------------------------
class Su2Adapter : public SolverAdapter {
public:
    Su2Adapter() = default;
    explicit Su2Adapter(const MappingRules& rules) : rules_(rules) {}

    const char* solver_name() const override { return "SU2"; }
    const char* format_name() const override { return "su2"; }

    // Full conversion: parse .cfg + .su2 mesh + results
    bool convert(const std::string& case_path,
                 ConversionResult& result) override;

    // Results-only import (solution.csv)
    bool import_results(const std::string& results_path,
                        ConversionResult& result) override;

    // Source detection
    bool detect_source(const std::string& case_path,
                       SourceInfo& info) override;

    // --- Individual file parsers ---
    bool parse_mesh(const std::string& su2_file);
    bool parse_config(const std::string& cfg_file);
    bool parse_solution_csv(const std::string& csv_file,
                            const std::string& field_name,
                            cfdx::core::ScalarCellField& scalar_field,
                            cfdx::core::Vec3CellField& vec_field);

    // --- Accessors ---
    const std::vector<cfdx::core::Vec3>& su2_points() const { return su2_points_; }
    const std::vector<Su2Boundary>& su2_boundaries() const { return boundaries_; }
    const CaseSetup& setup() const { return setup_; }

private:
    // --- .su2 mesh parsing — actual SU2 format uses NDIME=/NELEM=/NPOIN=/NMARK=
    bool parse_su2_points(std::istream& in, std::size_t n_points);
    bool parse_su2_elements(std::istream& in, std::size_t n_elements);
    bool parse_su2_boundary(std::istream& in, std::size_t n_markers);

    // --- .cfg parsing ---
    bool parse_cfg_value(const std::string& line, std::string& key, std::string& val);
    void parse_cfg_marker_bc(const std::string& value);
    void parse_cfg_numerics();

    // --- State ---
    MappingRules rules_ = MappingRules::load_default();
    CaseSetup setup_;

    // Parsed mesh data
    std::vector<cfdx::core::Vec3> su2_points_;
    std::vector<std::uint32_t> element_conn_;
    std::vector<int> element_types_;
    std::vector<std::size_t> element_offsets_;
    std::vector<Su2Boundary> boundaries_;

    // Parsed config data
    std::map<std::string, std::string> cfg_params_;
    std::vector<std::string> field_names_;
};

// Convenience API
bool parse_su2_mesh(const std::string& su2_file,
                    Su2Adapter& adapter,
                    cfdx::core::Mesh& mesh,
                    GapAnalysis& gap,
                    MappingRules rules = MappingRules::load_default());

bool parse_su2_config(const std::string& cfg_file,
                      Su2Adapter& adapter,
                      CaseSetup& setup,
                      GapAnalysis& gap);

}  // namespace su2
}  // namespace io
}  // namespace cfdx
