// M0.10-T03: Mapping rules implementation
#include "mapping_rules.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

namespace cfdx {
namespace io {

// ---------------------------------------------------------------------------
// Simple YAML-like parser (no external dependency)
// Recognised syntax:
//   bc:
//     "velocity-inlet": INLET
//     "pressure-outlet": PRESSURE_OUTLET
//   bcvt:
//     "velocity-inlet": FIXED
//   turbulence:
//     "kepsilon": k_epsilon
//   schemes:
//     "second-order-upwind": second_order_upwind
//   patches:
//     "wall": 0
// ---------------------------------------------------------------------------

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static std::string strip_quotes(const std::string& s) {
    std::string t = trim(s);
    if (t.size() >= 2) {
        if ((t.front() == '"' && t.back() == '"') ||
            (t.front() == '\'' && t.back() == '\'')) {
            return t.substr(1, t.size() - 2);
        }
    }
    return t;
}

MappingRules MappingRules::load_from_file(const std::string& path) {
    MappingRules rules = load_default();

    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[MappingRules] Warning: cannot open " << path
                  << ", using defaults only\n";
        return rules;
    }

    std::string line;
    std::string section = "";
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;

        // Check for section headers (indented or not)
        if (line.find("bc:") != std::string::npos) { section = "bc"; continue; }
        if (line.find("bcvt:") != std::string::npos) { section = "bcvt"; continue; }
        if (line.find("turbulence:") != std::string::npos) { section = "turbulence"; continue; }
        if (line.find("schemes:") != std::string::npos) { section = "schemes"; continue; }
        if (line.find("patches:") != std::string::npos) { section = "patches"; continue; }

        // Parse key: value
        size_t colon = line.find(':');
        if (colon == std::string::npos) continue;

        std::string key = strip_quotes(line.substr(0, colon));
        std::string val = strip_quotes(line.substr(colon + 1));
        if (key.empty() || val.empty()) continue;

        if (section == "bc") {
            if (val == "WALL") rules.add_bc_mapping(key, BCType::WALL);
            else if (val == "INLET") rules.add_bc_mapping(key, BCType::INLET);
            else if (val == "OUTLET") rules.add_bc_mapping(key, BCType::OUTLET);
            else if (val == "PRESSURE_OUTLET") rules.add_bc_mapping(key, BCType::PRESSURE_OUTLET);
            else if (val == "SYMMETRY") rules.add_bc_mapping(key, BCType::SYMMETRY);
            else if (val == "PERIODIC") rules.add_bc_mapping(key, BCType::PERIODIC);
            else if (val == "INTERFACE") rules.add_bc_mapping(key, BCType::INTERFACE);
            else if (val == "EMPTY") rules.add_bc_mapping(key, BCType::EMPTY);
            else if (val == "INTERNAL") rules.add_bc_mapping(key, BCType::INTERNAL);
        } else if (section == "bcvt") {
            if (val == "FIXED") rules.add_bcvt_mapping(key, BCValueType::FIXED);
            else if (val == "ZERO_GRADIENT") rules.add_bcvt_mapping(key, BCValueType::ZERO_GRADIENT);
            else if (val == "MIXED") rules.add_bcvt_mapping(key, BCValueType::MIXED);
            else if (val == "OUTLET_PRESSURE") rules.add_bcvt_mapping(key, BCValueType::OUTLET_PRESSURE);
            else if (val == "WALL_NO_SLIP") rules.add_bcvt_mapping(key, BCValueType::WALL_NO_SLIP);
            else if (val == "WALL_SLIP") rules.add_bcvt_mapping(key, BCValueType::WALL_SLIP);
            else if (val == "WALL_THERMAL") rules.add_bcvt_mapping(key, BCValueType::WALL_THERMAL);
        } else if (section == "turbulence") {
            rules.add_turbulence_mapping(key, val);
        } else if (section == "schemes") {
            rules.add_scheme_mapping(key, val);
        } else if (section == "patches") {
            try {
                rules.add_patch_mapping(key, std::stoi(val));
            } catch (...) {}
        }
    }

    return rules;
}

// ---------------------------------------------------------------------------
// Default mappings — Fluent, OpenFOAM, SU2 vendor strings
// ---------------------------------------------------------------------------
MappingRules MappingRules::load_default() {
    MappingRules rules;

    // --- Fluent boundary types ---
    rules.add_bc_mapping("velocity-inlet", BCType::INLET);
    rules.add_bc_mapping("mass-flow-inlet", BCType::INLET);
    rules.add_bc_mapping("pressure-inlet", BCType::INLET);
    rules.add_bc_mapping("pressure-outlet", BCType::PRESSURE_OUTLET);
    rules.add_bc_mapping("outflow", BCType::OUTLET);
    rules.add_bc_mapping("wall", BCType::WALL);
    rules.add_bc_mapping("symmetry", BCType::SYMMETRY);
    rules.add_bc_mapping("periodic", BCType::PERIODIC);
    rules.add_bc_mapping("interface", BCType::INTERFACE);
    rules.add_bc_mapping("degenerate", BCType::EMPTY);
    rules.add_bc_mapping("internal", BCType::INTERNAL);
    rules.add_bc_mapping(" porous-jump", BCType::WALL);
    rules.add_bc_mapping("radial-inlet", BCType::INLET);
    rules.add_bc_mapping("vacuum-pressure", BCType::PRESSURE_OUTLET);

    // Fluent BC value types
    rules.add_bcvt_mapping("velocity-inlet", BCValueType::FIXED);
    rules.add_bcvt_mapping("mass-flow-inlet", BCValueType::FIXED);
    rules.add_bcvt_mapping("pressure-inlet", BCValueType::FIXED);
    rules.add_bcvt_mapping("pressure-outlet", BCValueType::OUTLET_PRESSURE);
    rules.add_bcvt_mapping("outflow", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("wall", BCValueType::WALL_NO_SLIP);
    rules.add_bcvt_mapping("symmetry", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("periodic", BCValueType::MIXED);
    rules.add_bcvt_mapping("interface", BCValueType::MIXED);

    // --- OpenFOAM boundary types ---
    rules.add_bc_mapping("fixedValue", BCType::INLET);
    rules.add_bc_mapping("fixedValueOutlet", BCType::OUTLET);
    rules.add_bc_mapping("zeroGradient", BCType::INLET);
    rules.add_bc_mapping("inlet", BCType::INLET);
    rules.add_bc_mapping("outlet", BCType::OUTLET);
    rules.add_bc_mapping("wall", BCType::WALL);
    rules.add_bc_mapping("symmetryPlane", BCType::SYMMETRY);
    rules.add_bc_mapping("symmetry", BCType::SYMMETRY);
    rules.add_bc_mapping("empty", BCType::EMPTY);
    rules.add_bc_mapping("processor", BCType::INTERFACE);
    rules.add_bc_mapping("overset", BCType::INTERFACE);

    // OpenFOAM BC value types
    rules.add_bcvt_mapping("fixedValue", BCValueType::FIXED);
    rules.add_bcvt_mapping("zeroGradient", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("mixed", BCValueType::MIXED);
    rules.add_bcvt_mapping("outletInlet", BCValueType::OUTLET_PRESSURE);
    rules.add_bcvt_mapping("inletOutlet", BCValueType::FIXED);

    // --- SU2 boundary markers ---
    rules.add_bc_mapping(" inlet", BCType::INLET);
    rules.add_bc_mapping("outlet", BCType::PRESSURE_OUTLET);
    rules.add_bc_mapping("wall", BCType::WALL);
    rules.add_bc_mapping("symmetry", BCType::SYMMETRY);
    rules.add_bc_mapping("periodic", BCType::PERIODIC);
    rules.add_bc_mapping("interface", BCType::INTERFACE);
    rules.add_bc_mapping("far", BCType::OUTLET);
    rules.add_bc_mapping("neumann", BCType::OUTLET);

    // SU2 BC value types
    rules.add_bcvt_mapping("inlet", BCValueType::FIXED);
    rules.add_bcvt_mapping("outlet", BCValueType::OUTLET_PRESSURE);
    rules.add_bcvt_mapping("wall", BCValueType::WALL_NO_SLIP);
    rules.add_bcvt_mapping("symmetry", BCValueType::ZERO_GRADIENT);

    // --- STAR-CCM+ region/boundary types ---
    rules.add_bc_mapping("Inlet", BCType::INLET);
    rules.add_bc_mapping("Pressure Outlet", BCType::PRESSURE_OUTLET);
    rules.add_bc_mapping("Outlet", BCType::OUTLET);
    rules.add_bc_mapping("Wall", BCType::WALL);
    rules.add_bc_mapping("Symmetry", BCType::SYMMETRY);
    rules.add_bc_mapping("Periodic", BCType::PERIODIC);
    rules.add_bc_mapping("Interface", BCType::INTERFACE);
    rules.add_bc_mapping("Free Stream", BCType::INLET);
    rules.add_bc_mapping("Opening", BCType::PRESSURE_OUTLET);

    // STAR-CCM+ BC value types
    rules.add_bcvt_mapping("Inlet", BCValueType::FIXED);
    rules.add_bcvt_mapping("Pressure Outlet", BCValueType::OUTLET_PRESSURE);
    rules.add_bcvt_mapping("Outlet", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("Wall", BCValueType::WALL_NO_SLIP);
    rules.add_bcvt_mapping("Symmetry", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("Periodic", BCValueType::MIXED);
    rules.add_bcvt_mapping("Interface", BCValueType::MIXED);

    // --- Code_Saturne boundary types ---
    rules.add_bc_mapping("inlet", BCType::INLET);
    rules.add_bc_mapping("outlet", BCType::OUTLET);
    rules.add_bc_mapping("wall", BCType::WALL);
    rules.add_bc_mapping("symmetry", BCType::SYMMETRY);
    rules.add_bc_mapping("periodic", BCType::PERIODIC);
    rules.add_bc_mapping("interface", BCType::INTERFACE);
    rules.add_bc_mapping("far_field", BCType::OUTLET);
    rules.add_bc_mapping("opening", BCType::PRESSURE_OUTLET);

    // Code_Saturne BC value types
    rules.add_bcvt_mapping("inlet", BCValueType::FIXED);
    rules.add_bcvt_mapping("outlet", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("wall", BCValueType::WALL_NO_SLIP);
    rules.add_bcvt_mapping("symmetry", BCValueType::ZERO_GRADIENT);
    rules.add_bcvt_mapping("far_field", BCValueType::FIXED);
    rules.add_bcvt_mapping("opening", BCValueType::OUTLET_PRESSURE);

    // Code_Saturne patch types
    rules.add_patch_mapping("wall", 0);
    rules.add_patch_mapping("inlet", 1);
    rules.add_patch_mapping("outlet", 2);
    rules.add_patch_mapping("symmetry", 3);
    rules.add_patch_mapping("periodic", 4);
    rules.add_patch_mapping("interface", 5);
    rules.add_patch_mapping("far_field", 7);  // UNKNOWN
    rules.add_patch_mapping("opening", 2);

    // --- Turbulence model mappings ---
    // Fluent
    rules.add_turbulence_mapping("ke-realizable-viscous", "realizable_ke");
    rules.add_turbulence_mapping("ke-standard-viscous", "k_epsilon");
    rules.add_turbulence_mapping("ke-pseudo-viscous", "k_epsilon");
    rules.add_turbulence_mapping("kw-viscous", "k_omega_sst");
    rules.add_turbulence_mapping("mixing-length", "mixing_length");
    rules.add_turbulence_mapping("transition-to-turbulence", "transition");
    rules.add_turbulence_mapping("sST-ko-Transport", "k_omega_sst");

    // OpenFOAM
    rules.add_turbulence_mapping("kEpsilon", "k_epsilon");
    rules.add_turbulence_mapping("kOmegaSST", "k_omega_sst");
    rules.add_turbulence_mapping("kOmega", "k_omega");
    rules.add_turbulence_mapping("realizableKE", "realizable_ke");
    rules.add_turbulence_mapping("LRR", "reynolds_stress");
    rules.add_turbulence_mapping("Boussinesq", "buoyant");

    // SU2
    rules.add_turbulence_mapping("NONE", "laminar");
    rules.add_turbulence_mapping("SA", "spalart_allmaras");
    rules.add_turbulence_mapping("SST_SUTHERLAND", "k_omega_sst");
    rules.add_turbulence_mapping("SA_NTS", "spalart_allmaras_nts");

    // STAR-CCM+
    rules.add_turbulence_mapping("k-Epsilon", "k_epsilon");
    rules.add_turbulence_mapping("k-Omega SST", "k_omega_sst");
    rules.add_turbulence_mapping("Spalart-Allmaras", "spalart_allmaras");
    rules.add_turbulence_mapping("Laminar", "laminar");
    rules.add_turbulence_mapping("Realizable k-Epsilon", "realizable_ke");

    // --- Code_Saturne turbulence models ---
    rules.add_turbulence_mapping("KEPSILON", "k_epsilon");
    rules.add_turbulence_mapping("KEPSIL", "k_epsilon");
    rules.add_turbulence_mapping("KOMEGA", "k_omega");
    rules.add_turbulence_mapping("K-OMEGA-SST", "k_omega_sst");
    rules.add_turbulence_mapping("SPALART-ALLMARAS", "spalart_allmaras");
    rules.add_turbulence_mapping("Laminar", "laminar");
    rules.add_turbulence_mapping("laminar", "laminar");
    rules.add_turbulence_mapping("LES", "LES");
    rules.add_turbulence_mapping("RANS", "RANS");

    // --- Scheme mappings ---
    // Fluent
    rules.add_scheme_mapping("first-order", "first_order");
    rules.add_scheme_mapping("second-order-upwind", "second_order_upwind");
    rules.add_scheme_mapping("third-order-muscle", "third_order_muscle");
    rules.add_scheme_mapping("MUSCL", "MUSCL");
    rules.add_scheme_mapping("PRESTO!", "PRESTO");
    rules.add_scheme_mapping("STANDARD", "standard");
    rules.add_scheme_mapping("body-force-weighted", "body_force_weighted");
    rules.add_scheme_mapping("phase-coupled SIMPLE", "coupled");

    // OpenFOAM
    rules.add_scheme_mapping("Gauss linear", "linear");
    rules.add_scheme_mapping("Gauss linearUpwind", "linear_upwind");
    rules.add_scheme_mapping("Gauss upwind", "upwind");
    rules.add_scheme_mapping("Gauss limited", "limited");
    rules.add_scheme_mapping("Gauss LUST", "LUST");

    // SU2
    rules.add_scheme_mapping("GREEN_GAUSS_CELL_VOLUME", "green_gauss_cell");
    rules.add_scheme_mapping("LEAST_CONSERVATIVE", "least_squares");

    // --- Patch type mappings ---
    // Fluent zone types (use PatchType enum values from boundary.h)
    // WALL=0, INLET=1, OUTLET=2, SYMMETRY=3, PERIODIC=4, INTERFACE=5, EMPTY=6, UNKNOWN=7
    rules.add_patch_mapping("wall", 0);
    rules.add_patch_mapping("velocity-inlet", 1);
    rules.add_patch_mapping("pressure-inlet", 1);
    rules.add_patch_mapping("mass-flow-inlet", 1);
    rules.add_patch_mapping("pressure-outlet", 2);
    rules.add_patch_mapping("outflow", 2);
    rules.add_patch_mapping("symmetry", 3);
    rules.add_patch_mapping("periodic", 4);
    rules.add_patch_mapping("interface", 5);
    rules.add_patch_mapping("degenerate", 6);
    rules.add_patch_mapping("internal", 5);

    // OpenFOAM
    rules.add_patch_mapping("wall", 0);
    rules.add_patch_mapping("mappedWall", 0);
    rules.add_patch_mapping("inlet", 1);
    rules.add_patch_mapping("outlet", 2);
    rules.add_patch_mapping("symmetry", 3);
    rules.add_patch_mapping("symmetryPlane", 3);
    rules.add_patch_mapping("empty", 6);
    rules.add_patch_mapping("processor", 5);
    rules.add_patch_mapping("overset", 5);

    // STAR-CCM+
    rules.add_patch_mapping("Wall", 0);
    rules.add_patch_mapping("Inlet", 1);
    rules.add_patch_mapping("Pressure Outlet", 2);
    rules.add_patch_mapping("Outlet", 2);
    rules.add_patch_mapping("Symmetry", 3);
    rules.add_patch_mapping("Periodic", 4);
    rules.add_patch_mapping("Interface", 5);

    return rules;
}

}  // namespace io
}  // namespace cfdx
