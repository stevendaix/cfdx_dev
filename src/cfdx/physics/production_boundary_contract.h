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
        // INLET identifies the mesh patch. Its field operator may be either
        // an explicitly prescribed velocity or a zero-gradient open boundary
        // (used by body-force-driven channel cases). Imported pressure/mass-flow
        // inlets remain rejected above because they require distinct operators.
        if (bc.value_type != BCValueType::FIXED &&
            bc.value_type != BCValueType::ZERO_GRADIENT)
            reject("inlet requires FIXED velocity or ZERO_GRADIENT; pressure/debit inlet needs a dedicated operator");
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