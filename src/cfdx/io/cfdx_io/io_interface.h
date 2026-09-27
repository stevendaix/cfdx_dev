// M0.10-T03: Unified solver I/O converter interfaces
//
// Spécification CFDX v0.7 §21/22/23 (extension issue #425)
// ---------------------------------------------------------------------------
// Common abstraction layer for translating external solver file formats
// into the stable CFDX intermediate representation.
//
// Each external adapter implements SolverAdapter and produces a ConversionResult.
// Every unsupported or approximated feature is recorded in the GapAnalysis.
// No silent substitutions/fallbacks.
#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/boundary/boundary_field.h"
#include <string>
#include <vector>
#include <memory>
#include <cstddef>

namespace cfdx {
namespace io {

// ---------------------------------------------------------------------------
// Source identification (solver, version, case path)
// ---------------------------------------------------------------------------
struct SourceInfo {
    std::string solver;        // "Fluent", "OpenFOAM", "SU2", "STAR-CCM+", ...
    std::string version;       // version when known, empty otherwise
    std::string case_path;     // full path to the case directory or file
    std::string case_name;     // case identifier
    std::string format;        // "cas_dat", "sim", "su2_cfg", ...
};

// ---------------------------------------------------------------------------
// Severity levels for GapAnalysis entries
// ---------------------------------------------------------------------------
enum class Severity {
    SUPPORTED,          // mapped exactly
    APPROXIMATED,       // mapped with documented approximation
    UNSUPPORTED_NONBLOCK,  // not mapped, does not block execution
    UNSUPPORTED_BLOCK,     // not mapped, blocks execution
    UNAVAILABLE          // source information not present in file
};

inline const char* to_string(Severity s) {
    switch (s) {
        case Severity::SUPPORTED:             return "supported";
        case Severity::APPROXIMATED:          return "approximated";
        case Severity::UNSUPPORTED_NONBLOCK:  return "unsupported_nonblocking";
        case Severity::UNSUPPORTED_BLOCK:     return "unsupported_blocking";
        case Severity::UNAVAILABLE:           return "unavailable";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// A single gap-analysis finding
// ---------------------------------------------------------------------------
struct Finding {
    Severity severity;
    std::string category;   // e.g. "boundary_condition", "material", "numerical_scheme"
    std::string feature;    // e.g. "pressure outlet", "turbulence kepsilon"
    std::string detail;     // human-readable description
    std::string suggestion; // remediation
};

// ---------------------------------------------------------------------------
// GapAnalysis report — accumulates all findings from a conversion
// ---------------------------------------------------------------------------
class GapAnalysis {
public:
    GapAnalysis() = default;

    GapAnalysis& add(Severity sev, const std::string& cat,
                     const std::string& feat, const std::string& detail,
                     const std::string& suggestion = "") {
        findings_.push_back({sev, cat, feat, detail, suggestion});
        return *this;
    }

    // Convenience helpers
    GapAnalysis& supported(const std::string& cat, const std::string& feat,
                          const std::string& detail = "") {
        return add(Severity::SUPPORTED, cat, feat, detail);
    }
    GapAnalysis& approximated(const std::string& cat, const std::string& feat,
                              const std::string& detail, const std::string& sug = "") {
        return add(Severity::APPROXIMATED, cat, feat, detail, sug);
    }
    GapAnalysis& unsupported_nonblocking(const std::string& cat, const std::string& feat,
                                         const std::string& detail, const std::string& sug = "") {
        return add(Severity::UNSUPPORTED_NONBLOCK, cat, feat, detail, sug);
    }
    GapAnalysis& unsupported_blocking(const std::string& cat, const std::string& feat,
                                      const std::string& detail, const std::string& sug = "") {
        return add(Severity::UNSUPPORTED_BLOCK, cat, feat, detail, sug);
    }
    GapAnalysis& unavailable(const std::string& cat, const std::string& feat,
                             const std::string& detail) {
        return add(Severity::UNAVAILABLE, cat, feat, detail);
    }

    // Aggregators
    bool has_blocking() const {
        for (const auto& f : findings_)
            if (f.severity == Severity::UNSUPPORTED_BLOCK) return true;
        return false;
    }

    std::size_t n_supported() const {
        return count(Severity::SUPPORTED) + count(Severity::APPROXIMATED);
    }
    std::size_t n_approximated() const { return count(Severity::APPROXIMATED); }
    std::size_t n_unsupported_nonblocking() const { return count(Severity::UNSUPPORTED_NONBLOCK); }
    std::size_t n_unsupported_blocking() const { return count(Severity::UNSUPPORTED_BLOCK); }
    std::size_t n_unavailable() const { return count(Severity::UNAVAILABLE); }

    const std::vector<Finding>& findings() const { return findings_; }

    // Serialisation
    std::string to_json() const;
    std::string to_markdown() const;

private:
    std::size_t count(Severity s) const {
        std::size_t c = 0;
        for (const auto& f : findings_)
            if (f.severity == s) ++c;
        return c;
    }

    std::vector<Finding> findings_;
};

// ---------------------------------------------------------------------------
// ConversionResult — everything an adapter produces
// ---------------------------------------------------------------------------
struct ConversionResult {
    SourceInfo source;
    cfdx::core::Mesh mesh;
    std::vector<cfdx::core::ScalarCellField> scalar_fields;   // e.g. p, T
    std::vector<cfdx::core::Vec3CellField>   vec_fields;      // e.g. U
    GapAnalysis gap_report;

    bool success() const { return !gap_report.has_blocking(); }
    bool has_mesh() const { return mesh.n_cells() > 0; }
};

// ---------------------------------------------------------------------------
// SolverAdapter — abstract interface (pure virtual)
// ---------------------------------------------------------------------------
// Each external solver provides a concrete adapter.
// Adapters never touch CFDX numerical-core internals; they only populate
// the common CFDX Mesh, Fields, and GapAnalysis.
class SolverAdapter {
public:
    virtual ~SolverAdapter() = default;

    // --- Identification ---
    virtual const char* solver_name() const = 0;
    virtual const char* format_name() const = 0;

    // --- Core conversion ---
    // Parses the source files and populates result (mesh + fields + gap report).
    // Returns true if the conversion completed without blocking incompatibilities.
    virtual bool convert(const std::string& case_path,
                         ConversionResult& result) = 0;

    // --- Optional: results-only import (after mesh already loaded) ---
    virtual bool import_results(const std::string& results_path,
                                ConversionResult& result) {
        (void)results_path; (void)result;
        return false;  // not all adapters support results-only import
    }

    // --- Source/version detection ---
    virtual bool detect_source(const std::string& case_path, SourceInfo& info) {
        (void)case_path; (void)info;
        return false;
    }
};

}  // namespace io
}  // namespace cfdx
