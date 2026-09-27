#pragma once

#include "cfdx/core/boundary/boundary_constraint.h"
#include "cfdx/core/boundary/mathematical_condition.h"
#include "cfdx/core/boundary/boundary.h"
#include "cfdx/core/mesh/mesh.h"
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <memory>
#include <type_traits>
#include <variant>
#include <vector>

namespace cfdx::physics {

// Face-resolved mathematical constraints used by FVM assembly. The solver
// receives one resolved condition per global boundary face; no patch-name
// inference is performed here.
struct ScalarBoundaryFaceCondition {
    enum class Type { FIXED_VALUE, FIXED_GRADIENT };
    Type type = Type::FIXED_GRADIENT;
    double value = 0.0;
};

struct ScalarBoundaryFaceConditions {
    std::vector<ScalarBoundaryFaceCondition> conditions;
    std::vector<unsigned char> specified;

    explicit ScalarBoundaryFaceConditions(std::size_t n_faces = 0)
        : conditions(n_faces), specified(n_faces, 0) {}
    bool has(std::size_t face) const noexcept {
        return face < specified.size() && specified[face] != 0;
    }
};

using BoundaryConstraintMap =
    std::map<std::string, std::vector<cfdx::core::BoundaryConstraint>>;

inline double evaluate_boundary_provider(
    const std::shared_ptr<const cfdx::core::ValueProvider>& provider,
    const cfdx::core::Vec3& centre, std::size_t face_id, double time)
{
    if (!provider)
        throw std::invalid_argument("boundary constraint requires a value provider");
    cfdx::core::ValueContext context;
    context.time = time;
    context.x = centre.x;
    context.y = centre.y;
    context.z = centre.z;
    context.face_id = face_id;
    const double value = provider->evaluate(context);
    if (!std::isfinite(value))
        throw std::runtime_error("boundary constraint provider returned a non-finite value");
    return value;
}

// Resolve mathematical constraints onto the actual boundary faces of the mesh.
// This is deliberately strict: every requested field/patch must have exactly
// one supported scalar condition, and unsupported mathematical conditions are
// rejected rather than silently downgraded to a historical PatchField.
inline ScalarBoundaryFaceConditions resolve_scalar_boundary_constraints(
    const cfdx::core::Mesh& mesh,
    const std::vector<cfdx::core::Vec3>& face_centres,
    const BoundaryConstraintMap& constraints,
    const std::string& field,
    double time = 0.0)
{
    if (face_centres.size() != mesh.n_faces())
        throw std::invalid_argument("resolve_scalar_boundary_constraints: geometry/mesh mismatch");

    ScalarBoundaryFaceConditions resolved(mesh.n_faces());
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        const auto& patch = mesh.boundary().patch(p);
        if (patch.type == cfdx::core::PatchType::EMPTY)
            continue;
        const auto it = constraints.find(patch.name);
        if (it == constraints.end())
            throw std::invalid_argument(
                "missing mathematical boundary constraint for patch '" + patch.name + "'");

        const cfdx::core::BoundaryConstraint* selected = nullptr;
        for (const auto& constraint : it->second) {
            if (constraint.field != field) continue;
            if (selected != nullptr)
                throw std::invalid_argument(
                    "multiple mathematical boundary constraints for field '" + field +
                    "' on patch '" + patch.name + "'");
            selected = &constraint;
        }
        if (selected == nullptr)
            throw std::invalid_argument(
                "missing mathematical boundary constraint for field '" + field +
                "' on patch '" + patch.name + "'");

        std::visit([&](const auto& condition) {
            using T = std::decay_t<decltype(condition)>;
            if constexpr (std::is_same_v<T, cfdx::core::Dirichlet>) {
                for (const std::size_t face : patch.face_ids) {
                    if (face >= mesh.n_faces())
                        throw std::out_of_range("boundary patch contains an invalid face id");
                    resolved.conditions[face] = {
                        ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                        evaluate_boundary_provider(condition.value, face_centres[face], face, time)};
                    resolved.specified[face] = 1;
                }
            } else if constexpr (std::is_same_v<T, cfdx::core::Neumann>) {
                for (const std::size_t face : patch.face_ids) {
                    if (face >= mesh.n_faces())
                        throw std::out_of_range("boundary patch contains an invalid face id");
                    resolved.conditions[face] = {
                        ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,
                        evaluate_boundary_provider(condition.gradient, geometry.face_centres[face], face, time)};
                    resolved.specified[face] = 1;
                }
            } else if constexpr (std::is_same_v<T, cfdx::core::Robin>) {
                throw std::invalid_argument(
                    "Robin mathematical condition is not yet supported by direct scalar FVM assembly");
            } else if constexpr (std::is_same_v<T, cfdx::Flux>) {
                throw std::invalid_argument(
                    "Flux mathematical condition is deferred until flux-dependent assembly is implemented");
            } else if constexpr (std::is_same_v<T, cfdx::Mixed>) {
                throw std::invalid_argument(
                    "Mixed mathematical condition is deferred until mixed assembly is implemented");
            } else if constexpr (std::is_same_v<T, cfdx::Coupled>) {
                throw std::invalid_argument(
                    "Coupled mathematical condition is not valid for scalar local FVM assembly");
            } else if constexpr (std::is_same_v<T, cfdx::Periodic>) {
                throw std::invalid_argument(
                    "Periodic mathematical condition requires paired-face assembly");
            }
        }, selected->condition);
    }
    return resolved;
}

} // namespace cfdx::physics