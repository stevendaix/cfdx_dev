// Fail-closed translation guard for the HDF5 -> incompressible production path.
// This does not implement new mathematical boundary operators. It prevents
// physical BC labels from silently downgrading to a different BC.
#pragma once

#include "cfdx/core/mesh/boundary.h"
#include "cfdx/io/cfdx_io/case_schema.h"

#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <string>

namespace cfdx::physics {

// Reject duplicate names: CaseSetup::find_boundary() otherwise returns the
// first specification and silently ignores later edits in the same HDF5 case.
inline void validate_production_boundary_unique_names(
    const cfdx::io::CaseSetup& setup)
{
    std::unordered_set<std::string> names;
    for (const auto& bc : setup.boundary_conditions) {
        if (bc.patch_name.empty())
            throw std::invalid_argument("production boundary patch name cannot be empty");
        if (!names.insert(bc.patch_name).second)
            throw std::invalid_argument(
                "production boundary '" + bc.patch_name + "': duplicate case boundary definition");
    }
}

// The Workbench can persist independent U and p specifications. Before
// mapping them to FV operators, reject conflicting or incomplete requests;
// silently falling back to the legacy BC would solve a different problem.
inline void validate_native_flow_field_overrides(const cfdx::io::BoundarySpec& bc)
{
    using cfdx::io::BCType;
    auto check = [&](const char* field, std::size_t components) {
        const auto* native = bc.find_field(field);
        if (!native) return;
        const std::string where = "production boundary '" + bc.patch_name +
            "' native " + field + ": ";
        if (native->type != "FIXED_VALUE" && native->type != "ZERO_GRADIENT" &&
            !(components == 1 && native->type == "FIXED_GRADIENT"))
            throw std::invalid_argument(where + "unsupported operator");
        if (native->type == "FIXED_VALUE" && native->value.size() != components)
            throw std::invalid_argument(where + "wrong number of Dirichlet components");
        for (const double value : native->value)
            if (!std::isfinite(value))
                throw std::invalid_argument(where + "non-finite component");
        if (!std::isfinite(native->gradient))
            throw std::invalid_argument(where + "non-finite gradient");
        if (native->type == "FIXED_VALUE" && std::abs(native->gradient) > 0.0)
            throw std::invalid_argument(where + "Dirichlet and gradient cannot both be imposed");
    };
    check("velocity", 3);
    check("pressure", 1);
    if (const auto* u = bc.find_field("velocity")) {
        const bool prescribed = u->type == "FIXED_VALUE";
        if (prescribed && bc.type != BCType::INLET && bc.type != BCType::WALL)
            throw std::invalid_argument("fixed velocity requires an inlet or wall patch");
        if (!prescribed && (bc.type == BCType::INLET || bc.type == BCType::WALL))
            throw std::invalid_argument("inlet/wall velocity cannot be zeroGradient in the current FV adapter");
        if (bc.type == BCType::WALL && bc.value_type != cfdx::io::BCValueType::FIXED) {
            for (const auto value : u->value)
                if (value != 0.0)
                    throw std::invalid_argument("no-slip wall cannot have moving explicit velocity");
        }
    }
    if (const auto* pressure = bc.find_field("pressure")) {
        if (pressure->type == "FIXED_VALUE" && bc.type != BCType::PRESSURE_OUTLET)
            throw std::invalid_argument("fixed pressure requires pressure_outlet in the current FV adapter");
        if (pressure->type != "FIXED_VALUE" && bc.type == BCType::PRESSURE_OUTLET)
            throw std::invalid_argument("pressure_outlet cannot discard its fixed pressure");
        if (pressure->type == "FIXED_GRADIENT" &&
            bc.type != BCType::WALL && bc.type != BCType::OUTLET)
            throw std::invalid_argument("fixed pressure gradient unsupported on this physical patch");
    }
}

// The Workbench can persist independent U and p specifications. Before
// mapping them to FV operators, reject conflicting or incomplete requests;
// silently falling back to the legacy BC would solve a different problem.
inline void validate_native_flow_field_overrides(const cfdx::io::BoundarySpec& bc)
{
    using cfdx::io::BCType;
    auto check = [&](const char* field, std::size_t components) {
        const auto* native = bc.find_field(field);
        if (!native) return;
        const std::string where = "production boundary '" + bc.patch_name +
            "' native " + field + ": ";
        if (native->type != "FIXED_VALUE" && native->type != "ZERO_GRADIENT" &&
            !(components == 1 && native->type == "FIXED_GRADIENT"))
            throw std::invalid_argument(where + "unsupported operator");
        if (native->type == "FIXED_VALUE" && native->value.size() != components)
            throw std::invalid_argument(where + "wrong number of Dirichlet components");
        for (const double value : native->value)
            if (!std::isfinite(value))
                throw std::invalid_argument(where + "non-finite component");
        if (!std::isfinite(native->gradient))
            throw std::invalid_argument(where + "non-finite gradient");
        if (native->type == "FIXED_VALUE" && std::abs(native->gradient) > 0.0)
            throw std::invalid_argument(where + "Dirichlet and gradient cannot both be imposed");
    };
    check("velocity", 3);
    check("pressure", 1);
    if (const auto* u = bc.find_field("velocity")) {
        const bool prescribed = u->type == "FIXED_VALUE";
        if (prescribed && bc.type != BCType::INLET && bc.type != BCType::WALL)
            throw std::invalid_argument("fixed velocity requires an inlet or wall patch");
        if (!prescribed && (bc.type == BCType::INLET || bc.type == BCType::WALL))
            throw std::invalid_argument("inlet/wall velocity cannot be zeroGradient in the current FV adapter");
        if (bc.type == BCType::WALL && bc.value_type != cfdx::io::BCValueType::FIXED) {
            for (const auto value : u->value)
                if (value != 0.0)
                    throw std::invalid_argument("no-slip wall cannot have moving explicit velocity");
        }
    }
    if (const auto* pressure = bc.find_field("pressure")) {
        if (pressure->type == "FIXED_VALUE" && bc.type != BCType::PRESSURE_OUTLET)
            throw std::invalid_argument("fixed pressure requires pressure_outlet in the current FV adapter");
        if (pressure->type != "FIXED_VALUE" && bc.type == BCType::PRESSURE_OUTLET)
            throw std::invalid_argument("pressure_outlet cannot discard its fixed pressure");
        if (pressure->type == "FIXED_GRADIENT" &&
            bc.type != BCType::WALL && bc.type != BCType::OUTLET)
            throw std::invalid_argument("fixed pressure gradient unsupported on this physical patch");
    }
}

inline void validate_production_incompressible_boundary(
    const cfdx::io::BoundarySpec& bc, cfdx::core::PatchType mesh_patch_type)
{
    using cfdx::core::PatchType;
    using cfdx::io::BCType;
    using cfdx::io::BCValueType;

    const std::string where = "production boundary '" + bc.patch_name + "': ";
    auto reject = [&](const char* message) {
        throw std::invalid_argument(where + message);
    };
    if (bc.patch_name.empty()) reject("patch name cannot be empty");
    if (!std::isfinite(bc.velocity_magnitude) ||
        !std::isfinite(bc.pressure) || !std::isfinite(bc.back_pressure) ||
        !std::isfinite(bc.wall_temperature))
        reject("non-finite boundary value");
    if (!bc.velocity_vector.empty() && bc.velocity_vector.size() != 3)
        reject("velocity_vector must contain exactly 3 components");
    for (const double v : bc.velocity_vector)
        if (!std::isfinite(v)) reject("non-finite velocity component");

    const bool zero_velocity = bc.velocity_magnitude == 0.0 &&
        (bc.velocity_vector.empty() ||
         (bc.velocity_vector[0] == 0.0 &&
          bc.velocity_vector[1] == 0.0 &&
          bc.velocity_vector[2] == 0.0));

    switch (bc.type) {
    case BCType::WALL:
        if (mesh_patch_type != PatchType::WALL)
            reject("WALL case type requires a WALL mesh patch");
        if (bc.value_type == BCValueType::WALL_SLIP)
            reject("slip wall requires a tangential-stress/zero-normal-velocity operator; not supported by this production adapter");
        if (bc.value_type != BCValueType::WALL_NO_SLIP &&
            bc.value_type != BCValueType::WALL_THERMAL &&
            bc.value_type != BCValueType::FIXED)
            reject("wall requires wall_no_slip, wall_thermal or fixed moving-wall velocity");
        if ((bc.value_type == BCValueType::WALL_NO_SLIP ||
             bc.value_type == BCValueType::WALL_THERMAL) && !zero_velocity)
            reject("no-slip wall cannot prescribe non-zero motion; use FIXED moving-wall BC explicitly");
        return;
    case BCType::INLET:
        if (mesh_patch_type != PatchType::INLET)
            reject("INLET case type requires an INLET mesh patch");
        if (bc.source_type_name == "pressure-inlet" ||
            bc.source_type_name == "mass-flow-inlet" ||
            bc.source_type_name == "pressure-far-field")
            reject("imported inlet uses a physical operator not implemented by the incompressible production adapter");
        if (bc.value_type != BCValueType::FIXED)
            reject("velocity inlet requires FIXED velocity; pressure/debit inlet needs a dedicated operator");
        return;
    case BCType::PRESSURE_OUTLET:
        if (mesh_patch_type != PatchType::OUTLET)
            reject("PRESSURE_OUTLET requires an OUTLET mesh patch");
        if (bc.value_type != BCValueType::OUTLET_PRESSURE)
            reject("pressure outlet requires OUTLET_PRESSURE: otherwise p becomes zeroGradient");
        // `pressure` is the actual fixed pressure consumed by the solver.
        // Distinct `back_pressure` in the case would be silently ignored.
        if (bc.back_pressure != 0.0 && bc.back_pressure != bc.pressure)
            reject("back_pressure differs from pressure; solver only consumes pressure");
        return;
    case BCType::OUTLET:
        if (mesh_patch_type != PatchType::OUTLET)
            reject("OUTLET case type requires an OUTLET mesh patch");
        if (bc.source_type_name == "pressure-far-field")
            reject("Fluent compressible pressure-far-field is not a generic incompressible outlet");
        if (bc.value_type != BCValueType::ZERO_GRADIENT)
            reject("generic outlet requires ZERO_GRADIENT; use PRESSURE_OUTLET to fix pressure");
        return;
    case BCType::SYMMETRY:
        // The FV solver enforces phi=0 from PatchType::SYMMETRY, *not*
        // from its generic U zero-gradient patch-field data.
        if (mesh_patch_type != PatchType::SYMMETRY)
            reject("SYMMETRY case type requires a SYMMETRY mesh patch to enforce zero normal flux");
        if (bc.value_type != BCValueType::ZERO_GRADIENT)
            reject("symmetry requires ZERO_GRADIENT tangential fields");
        return;
    case BCType::EMPTY:
        if (mesh_patch_type != PatchType::EMPTY)
            reject("EMPTY case type requires an EMPTY mesh patch");
        if (bc.value_type != BCValueType::ZERO_GRADIENT &&
            bc.value_type != BCValueType::UNKNOWN)
            reject("EMPTY patch cannot impose an external flow condition");
        return;
    case BCType::PERIODIC:
    case BCType::INTERFACE:
        reject("periodic/interface face coupling is not implemented in this production path");
        return;
    case BCType::INTERNAL:
    case BCType::UNKNOWN:
    default:
        reject("unsupported physical boundary type");
    }
}

} // namespace cfdx::physics
