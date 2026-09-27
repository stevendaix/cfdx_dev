#pragma once
#include "cfdx/core/boundary/boundary_constraint.h"
#include "cfdx/core/boundary/boundary_field.h"
#include <cmath>
#include <memory>
#include <stdexcept>
namespace cfdx::core {

class LegacyPatchFieldAdapter {
public:
    static BoundaryConstraint to_constraint(const PatchField& patch, const std::string& field) {
        if (field.empty()) throw std::invalid_argument("legacy boundary adapter requires a field");
        if (patch.size() == 0) throw std::invalid_argument("legacy boundary adapter requires a non-empty patch");

        auto uniform_value = [&]() {
            const double value = patch(0);
            for (std::size_t i = 1; i < patch.size(); ++i) {
                if (patch(i) != value)
                    throw std::invalid_argument(
                        "legacy fixedValue patch is spatially varying; use a mapped/profile provider before migration");
            }
            return value;
        };

        if (patch.type() == "fixedValue" || patch.type() == "FIXED_VALUE") {
            return BoundaryConstraint(
                field, Dirichlet{std::make_shared<ConstantValueProvider>(uniform_value())});
        }

        if (patch.type() == "zero" || patch.type() == "zeroGradient" ||
            patch.type() == "ZERO_GRADIENT") {
            return BoundaryConstraint(
                field, Neumann{std::make_shared<ConstantValueProvider>(0.0)});
        }

        if (patch.type() == "robin" || patch.type() == "ROBIN") {
            const auto& r = patch.robin_coefficients();
            return BoundaryConstraint(
                field, Robin{r.alpha, r.beta,
                    std::make_shared<ConstantValueProvider>(r.gamma)});
        }

        throw std::invalid_argument(
            "unsupported legacy PatchField type '" + patch.type() + "'");
    }
};

} // namespace cfdx::core
