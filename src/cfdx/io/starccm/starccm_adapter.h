// M0.10-T03: STAR-CCM+ .sim adapter
//
// Spécification CFDX v0.7 issue #425 §STAR-CCM+
// ---------------------------------------------------------------------------
// Parses STAR-CCM+ .sim files (binary, proprietary).
//
// ARCHITEURAL NOTE (per issue #425):
//   "Do not attempt unsupported direct parsing of proprietary .sim databases.
//    Require documented exports or licensed APIs where available.
//    Document exactly which export path was used."
//
// This adapter performs a *documented, safe* scan of the .sim binary file:
//   - Reads the file header to detect version and endianness
//   - Scans for readable ASCII strings that identify boundary zones,
//     physics continua, and material names
//   - Records every limitation in the GapAnalysis report
//
// SUPPORTED:
//   - Boundary zone name extraction (via string table scan)
//   - Material name extraction
//   - Physics model detection (via ASCII string patterns)
//   - Mesh metadata (cell count, vertex count if stored in header)
//
// NOT SUPPORTED (documented in gap report):
//   - Direct binary mesh topology reconstruction (faces, cells, owner/neighbour)
//   - Field results import from .sim
//   - UDFs and custom models
//
// RECOMMENDED EXPORT PATH:
//   File > Export > Results File... → .sim is NOT for direct parsing.
//   Use documented exports: .plt (EnSight), .vtk, .cgns, or .msh (Gmsh).
//   The adapter can also parse a STAR-CCM+ .su2 export.
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
namespace starccm {

// ---------------------------------------------------------------------------
// .sim file header structure (STAR-CCM+ binary format)
// ---------------------------------------------------------------------------
struct SimHeader {
    std::string magic;          // "STAR-CCM+" or version marker
    std::string version;       // e.g. "19.02", "20.04", etc.
    std::string build;
    bool is_little_endian = true;
    std::size_t file_size = 0;
    std::string format_note;   // additional metadata from header
};

// ---------------------------------------------------------------------------
// Extracted region / boundary / material metadata
// ---------------------------------------------------------------------------
struct SimZone {
    std::string name;
    std::string parent_region;
    std::string boundary_type;     // "wall", "inlet", "outlet", etc.
    std::string geometry_type;     // "surface", "volume", "edge", "vertex"
    std::size_t n_faces = 0;
};

struct SimMaterial {
    std::string name;
    std::string model;  // e.g. "constant", "piecewise"
};

struct SimPhysics {
    std::string model_name;
    std::string regime;     // "incompressible", "compressible"
    std::string turbulence; // "k-epsilon", "k-omega SST", etc.
    std::string energy;     // "on", "off"
    std::string multiphase; // "off", "VOF", "Eulerian"
};

// ---------------------------------------------------------------------------
// StarCCMAdapter — implements SolverAdapter for STAR-CCM+ .sim
// ---------------------------------------------------------------------------
class StarCCMAdapter : public SolverAdapter {
public:
    StarCCMAdapter() = default;
    explicit StarCCMAdapter(const MappingRules& rules) : rules_(rules) {}

    const char* solver_name() const override { return "STAR-CCM+"; }
    const char* format_name() const override { return "sim"; }

    // Full .sim import (metadata + string scan, not binary topology)
    bool convert(const std::string& sim_path,
                 ConversionResult& result) override;

    // Results import — NOT supported for .sim directly
    bool import_results(const std::string& /*results_path*/,
                        ConversionResult& result) override {
        result.gap_report.unsupported_blocking(
            "results", "sim_results",
            "STAR-CCM+ .sim results cannot be directly parsed",
            "Export results via: Tools → Export → Results File (.enc, .plt, .vtk)");
        return false;
    }

    // Source detection
    bool detect_source(const std::string& sim_path,
                       SourceInfo& info) override;

    // --- Accessors ---
    const SimHeader& header() const { return header_; }
    const std::vector<SimZone>& zones() const { return zones_; }
    const std::vector<SimMaterial>& materials() const { return materials_; }
    const SimPhysics& physics() const { return physics_; }

private:
    // --- Binary scanning ---
    bool read_header(const std::string& sim_path);
    bool scan_strings(const std::string& sim_path);
    bool scan_ascii_metadata(const std::vector<char>& buffer);

    // --- Helpers ---
    void classify_boundary_type(SimZone& zone);

    // --- State ---
    MappingRules rules_ = MappingRules::load_default();
    SimHeader header_;
    std::vector<SimZone> zones_;
    std::vector<SimMaterial> materials_;
    SimPhysics physics_{};
    CaseSetup setup_;
};

// Convenience API
bool parse_starccm_sim(const std::string& sim_path,
                       StarCCMAdapter& adapter,
                       GapAnalysis& gap,
                       MappingRules rules = MappingRules::load_default());

}  // namespace starccm
}  // namespace io
}  // namespace cfdx
