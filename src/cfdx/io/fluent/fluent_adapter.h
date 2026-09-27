// M0.10-T03: Fluent .cas/.dat adapter
//
// Spécification CFDX v0.7 issue #425 §Fluent adapter
// ---------------------------------------------------------------------------
// Parses Fluent legacy ASCII .cas (mesh + setup) and .dat (results) files.
//
// Supported .cas sections (legacy ASCII format):
//   - Section 1:  Header
//   - Section 2:  Nodes (coordinates)
//   - Section 3:  Face nodes (with cell tree)
//   - Section 4:  Faces (with owner/neighbour)
//   - Section 5:  Cells (cell-node connectivity)
//   - Section 31: Boundary conditions (zone + BC type info)
//   - Section 39: Materials
//   - Section 45: Cell zone conditions
//   - Section 48: Model settings (turbulence, transient, etc.)
//
// Supported .dat sections (results):
//   - Section 19:  Header (variable count)
//   - Cell data:   Pressure, velocity components, temperature
//
// Unsupported: .cas binary format, DPM, UDFs, custom field functions,
//   multiphase (VOF/Mixture/Eulerian), sliding mesh, dynamic mesh,
//   Fluent 202x new-native binary .cas.h5 (HDF5-based).
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
namespace fluent {

// ---------------------------------------------------------------------------
// Fluent .cas section identifier (legacy format)
// ---------------------------------------------------------------------------
enum class CasSection : int {
    HEADER           = 1,
    FACE_NODES       = 3,
    FACES            = 4,
    CELLS            = 5,
    BOUNDARY_ZONES   = 31,
    MATERIALS        = 39,
    CELL_ZONE_COND   = 45,
    MODEL_SETTINGS   = 48,
    CUSTOM_FUNCTIONS = 64,
};

// ---------------------------------------------------------------------------
// Fluent zone type (maps to boundary conditions)
// ---------------------------------------------------------------------------
enum class ZoneType : int {
    INTERIOR       = 0,
    VELOCITY_INLET = 1,
    PRESSURE_INLET = 2,
    PRESSURE_OUTLET= 3,
    WALL           = 4,
    SYMMETRY       = 5,
    PERIODIC       = 6,
    EXTERIOR       = 7,
    INTERNAL       = 8,
    MASS_OUTFLOW   = 9,
    PRESSURE_FAR_FIELD = 10,
    OUTLET_VENT    = 11,
    OVERSET        = 12,
    DEGENERATE     = 13,
    SLIP_WALL      = 14,
    UNKNOWN        = 99,
};

inline const char* to_string(ZoneType t) {
    switch (t) {
        case ZoneType::INTERIOR:          return "interior";
        case ZoneType::VELOCITY_INLET:    return "velocity-inlet";
        case ZoneType::PRESSURE_INLET:    return "pressure-inlet";
        case ZoneType::PRESSURE_OUTLET:   return "pressure-outlet";
        case ZoneType::WALL:              return "wall";
        case ZoneType::SYMMETRY:          return "symmetry";
        case ZoneType::PERIODIC:          return "periodic";
        case ZoneType::EXTERIOR:          return "exterior";
        case ZoneType::INTERNAL:          return "internal";
        case ZoneType::MASS_OUTFLOW:      return "mass-outflow";
        case ZoneType::PRESSURE_FAR_FIELD:return "pressure-far-field";
        case ZoneType::OUTLET_VENT:       return "outlet-vent";
        case ZoneType::OVERSET:           return "overset";
        case ZoneType::DEGENERATE:        return "degenerate";
        case ZoneType::SLIP_WALL:         return "slip-wall";
        default:                          return "unknown";
    }
}

// ---------------------------------------------------------------------------
// Fluent zone descriptor (read from section 31)
// ---------------------------------------------------------------------------
struct ZoneInfo {
    int zone_id = 0;
    int zone_type = 0;       // numeric fluent zone type
    std::string name;
    std::string node_type;   // "node", "face", "cell", "particle"
    std::string element_type; // "tri", "quad", "tet", "hex", etc.
    int n_elements = 0;
    std::vector<std::uint32_t> face_ids;   // for boundary face zones
    bool is_boundary() const {
        return node_type == "face";
    }
    bool is_cell_zone() const {
        return node_type == "cell";
    }
};

// ---------------------------------------------------------------------------
// FluentAdapter — implements SolverAdapter for Fluent .cas/.dat
// ---------------------------------------------------------------------------
class FluentAdapter : public SolverAdapter {
public:
    FluentAdapter() = default;
    explicit FluentAdapter(const MappingRules& rules) : rules_(rules) {}

    const char* solver_name() const override { return "Fluent"; }
    const char* format_name() const override { return "cas_dat"; }

    // Full case + results import
    bool convert(const std::string& case_path,
                 ConversionResult& result) override;

    // Results-only import (after mesh loaded)
    bool import_results(const std::string& results_path,
                        ConversionResult& result) override;

    // Source detection
    bool detect_source(const std::string& case_path,
                       SourceInfo& info) override;

    // --- Accessors for parsed zones (useful for tests) ---
    const std::vector<ZoneInfo>& zones() const { return zones_; }
    const CaseSetup& setup() const { return setup_; }

private:
    // --- .cas parsing ---
    bool parse_cas(const std::string& cas_file);
    bool parse_cas_header(std::istream& in);
    bool parse_cas_nodes(std::istream& in);
    bool parse_cas_faces(std::istream& in);
    bool parse_cas_cells(std::istream& in);
    bool parse_cas_zones(std::istream& in);
    bool parse_cas_materials(std::istream& in);
    bool parse_cas_models(std::istream& in);

    // --- .dat parsing (results) ---
    bool parse_dat(const std::string& dat_file);
    bool parse_dat_header(std::istream& in, std::vector<std::string>& variables);
    bool parse_dat_cell_data(std::istream& in,
                             const std::string& field_name,
                             std::size_t n_cells);

    // --- Helpers ---
    ZoneType zone_type_from_int(int t) const;
    void build_cell_topology();

    // --- State ---
    MappingRules rules_ = MappingRules::load_default();
    CaseSetup setup_;
    std::vector<ZoneInfo> zones_;
    std::vector<cfdx::core::Vec3> points_;
    std::vector<std::vector<std::uint32_t>> face_nodes_;
    std::vector<int> face_owner_;
    std::vector<int> face_neighbour_;
    std::vector<std::vector<std::uint32_t>> cell_faces_;
};

// Convenience API for standalone .cas parsing
bool parse_fluent_cas(const std::string& cas_file,
                      FluentAdapter& adapter,
                      cfdx::core::Mesh& mesh,
                      GapAnalysis& gap,
                      MappingRules rules = MappingRules::load_default());

// Convenience API for standalone .dat parsing
bool parse_fluent_dat(const std::string& dat_file,
                      cfdx::core::ScalarCellField& pressure_field,
                      cfdx::core::Vec3CellField& velocity_field,
                      GapAnalysis& gap);

}  // namespace fluent
}  // namespace io
}  // namespace cfdx
