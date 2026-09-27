// M0.10-T03: Mapping rules loader and lookup
//
// Spécification CFDX v0.7 issue #425 §Mapping rules
// ---------------------------------------------------------------------------
// Mapping rules translate vendor-specific strings to CFDX internal enums.
// Rules are externalised so they can be updated without recompilation.
#pragma once

#include "case_schema.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <optional>

namespace cfdx {
namespace io {

// ---------------------------------------------------------------------------
// RuleSet — maps vendor-specific strings to CFDX internal types
// ---------------------------------------------------------------------------
class MappingRules {
public:
    // Load from a YAML file (simple key: value parser — no external YAML dependency)
    static MappingRules load_from_file(const std::string& path);
    static MappingRules load_default();

    // --- Lookups ---

    // Boundary condition type mapping
    BCType map_boundary_type(const std::string& vendor_type) const {
        if (auto it = bc_map_.find(vendor_type); it != bc_map_.end())
            return it->second;
        return BCType::UNKNOWN;
    }

    // BC value type mapping
    BCValueType map_bc_value_type(const std::string& vendor_type) const {
        if (auto it = bcvt_map_.find(vendor_type); it != bcvt_map_.end())
            return it->second;
        return BCValueType::UNKNOWN;
    }

    // Turbulence model mapping
    std::string map_turbulence(const std::string& vendor_turbulence) const {
        if (auto it = turbulence_map_.find(vendor_turbulence); it != turbulence_map_.end())
            return it->second;
        return "";
    }

    // Numerical scheme mapping
    std::string map_scheme(const std::string& vendor_scheme) const {
        if (auto it = scheme_map_.find(vendor_scheme); it != scheme_map_.end())
            return it->second;
        return "";
    }

    // Patch type mapping (returns CFDX PatchType from core/mesh)
    int map_patch_type(const std::string& vendor_patch_type) const {
        if (auto it = patch_map_.find(vendor_patch_type); it != patch_map_.end())
            return it->second;
        return -1;  // UNKNOWN (use PatchType::UNKNOWN)
    }

    // --- Modifiers ---
    void add_bc_mapping(const std::string& vendor, BCType cfdx_type) {
        bc_map_[vendor] = cfdx_type;
    }
    void add_bcvt_mapping(const std::string& vendor, BCValueType vtype) {
        bcvt_map_[vendor] = vtype;
    }
    void add_turbulence_mapping(const std::string& vendor, const std::string& cfdx_name) {
        turbulence_map_[vendor] = cfdx_name;
    }
    void add_scheme_mapping(const std::string& vendor, const std::string& cfdx_name) {
        scheme_map_[vendor] = cfdx_name;
    }
    void add_patch_mapping(const std::string& vendor, int cfdx_patch_type) {
        patch_map_[vendor] = cfdx_patch_type;
    }

    bool has_bc_mapping(const std::string& vendor) const {
        return bc_map_.find(vendor) != bc_map_.end();
    }

    std::size_t n_bc_mappings() const { return bc_map_.size(); }
    std::size_t n_turbulence_mappings() const { return turbulence_map_.size(); }
    std::size_t n_scheme_mappings() const { return scheme_map_.size(); }

private:
    std::unordered_map<std::string, BCType> bc_map_;
    std::unordered_map<std::string, BCValueType> bcvt_map_;
    std::unordered_map<std::string, std::string> turbulence_map_;
    std::unordered_map<std::string, std::string> scheme_map_;
    std::unordered_map<std::string, int> patch_map_;
};

}  // namespace io
}  // namespace cfdx
