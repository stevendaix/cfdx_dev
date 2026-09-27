#pragma once

#include "cfdx/core/boundary/boundary.h"
#include "cfdx/core/boundary/boundary_constraint.h"
#include "cfdx/core/boundary/value_provider.h"
#include "cfdx/physics/steady_incompressible_solver.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

/// Transitional bridge from the new backend-independent BC contract to the
/// existing incompressible solver boundary maps.
///
/// This adapter is deliberately conservative:
/// - only constant providers are lowered to the legacy patch-wide maps;
/// - time/spatial/table providers are rejected instead of being silently
///   collapsed to a single value;
/// - mathematical semantics are preserved (Dirichlet/Neumann);
/// - patch names never infer a physical boundary condition.
///
/// The bridge is temporary. Face-aware provider evaluation and direct FVM
/// constraint assembly belong to the next migration stage.
struct BoundaryConstraintBridge {
    VelocityBoundaryConditions velocity;
    ScalarBoundaryConditions pressure;
};

inline double constant_provider_value(
    const std::shared_ptr<const cfdx::core::ValueProvider>& provider,
    const char* field)
{
    if (!provider)
        throw std::invalid_argument("boundary bridge: null provider");

    const cfdx::core::ValueContext context{};
    const double value = provider->evaluate(context);

    // Do not accept a non-constant provider merely because it can be evaluated
    // at the default context. This would silently destroy time/spatial semantics.
    if (dynamic_cast<const cfdx::core::ConstantValueProvider*>(provider.get()) == nullptr)
        throw std::invalid_argument(
            std::string("boundary bridge: non-constant provider for ") + field +
            " requires face/time-aware assembly");

    if (!std::isfinite(value))
        throw std::invalid_argument(
            std::string("boundary bridge: non-finite value for ") + field);
    return value;
}

inline void add_constraint(
    BoundaryConstraintBridge& out,
    const cfdx::core::BoundaryConstraint& constraint)
{
    using MC = cfdx::core::MathematicalCondition;

    if (constraint.field == "U.x" || constraint.field == "U.y" ||
        constraint.field == "U.z") {
        const double value = std::visit(
            [&](const auto& condition) -> double {
                using T = std::decay_t<decltype(condition)>;
                if constexpr (std::is_same_v<T, cfdx::core::Dirichlet>) {
                    return constant_provider_value(condition.value, constraint.field.c_str());
                } else if constexpr (std::is_same_v<T, cfdx::core::Neumann>) {
                    throw std::invalid_argument(
                        "boundary bridge: velocity Neumann requires direct FVM assembly");
                } else {
                    throw std::invalid_argument(
                        "boundary bridge: unsupported velocity mathematical condition");
                }
            }, constraint.condition);

        auto& bc = out.velocity[constraint.field == "U.x" ? "x" :
                                 constraint.field == "U.y" ? "y" : "z"];
        bc.type = VelocityBoundaryCondition::Type::FIXED_VALUE;
        if (constraint.field == "U.x") bc.value.x = value;
        if (constraint.field == "U.y") bc.value.y = value;
        if (constraint.field == "U.z") bc.value.z = value;
        return;
    }

    if (constraint.field == "p") {
        std::visit(
            [&](const auto& condition) {
                using T = std::decay_t<decltype(condition)>;
                if constexpr (std::is_same_v<T, cfdx::core::Dirichlet>) {
                    out.pressure["default"] = {
                        ScalarBoundaryType::FIXED_VALUE,
                        constant_provider_value(condition.value, "p")};
                } else if constexpr (std::is_same_v<T, cfdx::core::Neumann>) {
                    out.pressure["default"] = {
                        ScalarBoundaryType::ZERO_GRADIENT, 0.0};
                } else {
                    throw std::invalid_argument(
                        "boundary bridge: unsupported pressure mathematical condition");
                }
            }, constraint.condition);
        return;
    }

    throw std::invalid_argument(
        "boundary bridge: unsupported field '" + constraint.field + "'");
}

/// Lower a physical BC's constraints for the transitional legacy solver path.
///
/// The returned maps use the reserved key "default" for scalar pressure and
/// are intended only as an intermediate contract. Patch-aware application is
/// deliberately left to the solver migration PR so no patch-name inference is
/// introduced here.
inline BoundaryConstraintBridge lower_boundary_constraints(
    const cfdx::core::Boundary& boundary,
    const cfdx::core::BoundaryCondition& condition)
{
    BoundaryConstraintBridge out;
    const auto constraints = condition.constraints(boundary);
    for (const auto& c : constraints)
        add_constraint(out, c);
    return out;
}

} // namespace cfdx::physics
