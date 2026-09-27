// M0.10-T03: Typed case setup/configuration schema
//
// Spécification CFDX v0.7 §6, §8 (case configuration) + issue #425
// ---------------------------------------------------------------------------
// Versioned, strongly-typed internal schema for case setup.
// Adapters populate this schema from external formats.
// The schema is versioned for forward/backward compatibility.
#pragma once

#include <string>
#include <vector>
#include <map>
#include <cstddef>
#include <cstdint>
#include <optional>

namespace cfdx {
namespace io {

// Schema version — incremented when the schema changes in a non-backward-compatible way
static constexpr int CFDX_SCHEMA_VERSION = 1;

// ---------------------------------------------------------------------------
// Boundary condition types (CFDX internal taxonomy)
// ---------------------------------------------------------------------------
enum class BCType : std::uint8_t {
    WALL = 0,
    INLET,
    OUTLET,
    PRESSURE_OUTLET,
    SYMMETRY,
    PERIODIC,
    INTERFACE,
    EMPTY,
    INTERNAL,       // for internal faces/zones
    UNKNOWN
};

inline const char* to_string(BCType t) {
    switch (t) {
        case BCType::WALL:            return "wall";
        case BCType::INLET:           return "inlet";
        case BCType::OUTLET:          return "outlet";
        case BCType::PRESSURE_OUTLET: return "pressure_outlet";
        case BCType::SYMMETRY:        return "symmetry";
        case BCType::PERIODIC:        return "periodic";
        case BCType::INTERFACE:       return "interface";
        case BCType::EMPTY:           return "empty";
        case BCType::INTERNAL:        return "internal";
        default:                      return "unknown";
    }
}

// ---------------------------------------------------------------------------
// Boundary condition value specification
// ---------------------------------------------------------------------------
enum class BCValueType : std::uint8_t {
    FIXED,            // fixed value
    ZERO_GRADIENT,    // zero gradient / Neumann
    MIXED,            // mixed / combined
    OUTLET_PRESSURE,  // pressure outlet
    WALL_NO_SLIP,     // no-slip wall
    WALL_SLIP,        // slip wall
    WALL_THERMAL,     // thermally coupled wall
    UNKNOWN
};

// ---------------------------------------------------------------------------
// Material property specification
// ---------------------------------------------------------------------------
struct MaterialSpec {
    std::string name;                // e.g. "air", "water"
    double density = 0.0;            // kg/m³ (or reference value)
    double dynamic_viscosity = 0.0;  // Pa·s
    double thermal_conductivity = 0.0; // W/(m·K)
    double specific_heat = 0.0;      // J/(kg·K)
    double molecular_weight = 0.0;   // kg/mol
    std::string eos_model;           // "ideal_gas", "BWR", ...
    std::string turbulence_model;    // "laminar", "kepsilon", "komega_sst", ...
    std::map<std::string, std::string> extra_properties; // vendor-specific extras
};

// ---------------------------------------------------------------------------
// Boundary condition specification
// ---------------------------------------------------------------------------
struct BoundarySpec {
    std::string patch_name;
    BCType type = BCType::UNKNOWN;
    BCValueType value_type = BCValueType::UNKNOWN;

    // For inlet
    double velocity_magnitude = 0.0;
    std::vector<double> velocity_vector;  // [vx, vy, vz]
    double temperature = 0.0;
    double pressure = 0.0;
    double turbulence_intensity = 0.0;
    double turbulence_length_scale = 0.0;

    // For wall
    double roughness_height = 0.0;
    double wall_temperature = 0.0;

    // For outlet / pressure outlet
    double back_pressure = 0.0;

    // Metadata from source
    std::string source_zone_name;    // original zone/patch name in source
    std::string source_type_name;    // original type string in source (e.g. "velocity-inlet")
    std::map<std::string, std::string> raw_params;  // unmapped raw parameters

    bool has_raw_unmapped() const { return !raw_params.empty(); }
};

// ---------------------------------------------------------------------------
// Initial condition specification
// ---------------------------------------------------------------------------
struct InitialCondition {
    double velocity = 0.0;
    double pressure = 0.0;
    double temperature = 0.0;
    std::vector<double> velocity_vector;
    std::map<std::string, double> scalar_fields;  // e.g. "turbulentKE": 1e-4
};

// ---------------------------------------------------------------------------
// Numerical scheme settings
// ---------------------------------------------------------------------------
struct NumericalScheme {
    std::string momentum_scheme;     // "first_order", "second_order_upwind", "MUSCL"
    std::string pressure_scheme;     // "PRESTO!", "STANDARD", "body_force_weighted"
    std::string momentum_interpolation;  // "linear", "upwind", "SecondOrderUpwind"
    std::string transient_scheme;    // "steady", "unsteady_first_order", "unsteady_second_order_RungeKutta"
    std::string gradient_operator;   // "Green-Gauss cell based", "Green-Gauss node based", "least squares"
    double under_relaxation_momentum = 0.7;
    double under_relaxation_pressure = 0.3;
    std::string coupled_solver;      // "SIMPLE", "SIMPLEC", "PISO", "PIMPLE", "coupled"
    std::string preconditioner;      // for linear solvers
    std::string residual_target = "1e-5";
    int max_iterations = 500;
    std::map<std::string, std::string> raw_settings;
};

// ---------------------------------------------------------------------------
// Mesh metadata (quality findings, topology summary)
// ---------------------------------------------------------------------------
struct MeshMetadata {
    std::size_t n_vertices = 0;
    std::size_t n_faces = 0;
    std::size_t n_cells = 0;
    std::size_t n_boundary_faces = 0;
    std::size_t n_internal_faces = 0;
    std::size_t n_patches = 0;
    int dimension = 3;
    bool has_polyhedral = false;
    bool has_nonplanar_faces = false;
    std::string mesh_type;  // "structured", "unstructured", "polyhedral"
    std::vector<std::string> cell_types;  // e.g. "tetra", "hex", "prism", "pyramid"
};

// ---------------------------------------------------------------------------
// Complete typed case setup (the CFDX intermediate representation)
// ---------------------------------------------------------------------------
struct CaseSetup {
    // Schema version
    int schema_version = CFDX_SCHEMA_VERSION;

    // Case identification
    SourceInfo source;

    // Mesh metadata (populated by adapter, validated against CFDX Mesh)
    MeshMetadata mesh_info;

    // Physics
    std::string physics_model;            // "incompressible_laminar", "compressible", "buoyant"
    std::string turbulence_model;         // "laminar", "k_epsilon", "k_omega_sst", "realizable_ke"
    std::string energy_model;             // "none", "isothermal", "buoyant", "conjugate"
    std::string multiphase_model;         // "none", "VOF", "Eulerian", "mixture"
    std::string radiation_model;          // "none", "DO", "solar", "S2S"
    bool transient = false;
    double time_step = 0.0;
    double end_time = 0.0;
    int max_time_steps = 0;

    // Materials
    std::vector<MaterialSpec> materials;

    // Boundary conditions
    std::vector<BoundarySpec> boundary_conditions;

    // Initial conditions
    InitialCondition initial_condition;

    // Numerical schemes
    NumericalScheme numerics;

    // Reference values
    double ref_length = 1.0;
    double ref_density = 1.0;
    double ref_velocity = 1.0;
    double ref_temperature = 293.15;
    double ref_pressure = 101325.0;

    std::string gravity_vector;  // e.g. "0 0 -9.81"
    std::string units;           // "SI"
    std::string solver_mode;     // "steady", "unsteady"

    // Extra metadata preserved from source
    std::map<std::string, std::string> source_metadata;

    // --- Accessors ---

    const BoundarySpec* find_boundary(const std::string& name) const {
        for (const auto& bc : boundary_conditions)
            if (bc.patch_name == name) return &bc;
        return nullptr;
    }

    const MaterialSpec* find_material(const std::string& name) const {
        for (const auto& m : materials)
            if (m.name == name) return &m;
        return nullptr;
    }
};

}  // namespace io
}  // namespace cfdx
