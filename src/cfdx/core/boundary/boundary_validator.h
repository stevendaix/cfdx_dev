#pragma once
#include "cfdx/core/boundary/boundary_constraint.h"
#include <algorithm>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cfdx::core {

struct BoundaryRequirement {
    std::string field;
    bool require_constraint = true;
};

class BoundaryValidator {
public:
    void require_field(std::string field) { requirements_.push_back({std::move(field), true}); }

    std::vector<std::string> validate(const std::vector<BoundaryConstraint>& constraints) const {
        std::vector<std::string> errors;
        std::unordered_set<std::string> seen;

        for (const auto& constraint : constraints) {
            if (constraint.field.empty()) {
                errors.push_back("boundary constraint has an empty field name");
                continue;
            }
            if (!seen.insert(constraint.field).second)
                errors.push_back("duplicate boundary constraint for field '" + constraint.field + "'");
            try {
                validate_condition(constraint.condition);
            } catch (const std::exception& e) {
                errors.push_back("invalid boundary constraint for field '" +
                                 constraint.field + "': " + e.what());
            }
        }

        for (const auto& req : requirements_) {
            const bool found = std::any_of(
                constraints.begin(), constraints.end(),
                [&](const auto& c) { return c.field == req.field; });
            if (req.require_constraint && !found)
                errors.push_back("missing boundary constraint for field '" + req.field + "'");
        }
        return errors;
    }

    bool valid(const std::vector<BoundaryConstraint>& constraints) const {
        return validate(constraints).empty();
    }

private:
    std::vector<BoundaryRequirement> requirements_;
};

} // namespace cfdx::core
