#pragma once
#include "cfdx/core/boundary/boundary_constraint.h"
#include <string>
#include <vector>
namespace cfdx::core {
struct BoundaryRequirement { std::string field; bool require_constraint = true; };
class BoundaryValidator {
public:
    void require_field(std::string field) { requirements_.push_back({std::move(field), true}); }
    std::vector<std::string> validate(const std::vector<BoundaryConstraint>& constraints) const {
        std::vector<std::string> errors;
        for (const auto& req : requirements_) {
            bool found = false;
            for (const auto& c : constraints) if (c.field == req.field) { found = true; break; }
            if (req.require_constraint && !found) errors.push_back("missing boundary constraint for field '" + req.field + "'");
        }
        return errors;
    }
private:
    std::vector<BoundaryRequirement> requirements_;
};
}