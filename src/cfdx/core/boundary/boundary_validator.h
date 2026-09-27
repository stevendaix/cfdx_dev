#pragma once
#include "cfdx/core/boundary/boundary_constraint.h"
#include <algorithm>
#include <map>
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

    // Validate cross-patch references without assuming any physical meaning
    // from patch names. Periodic links must be reciprocal; coupled links must
    // resolve to a declared patch in the same boundary scope.
    std::vector<std::string> validate_references(
        const std::map<std::string, std::vector<BoundaryConstraint>>& by_patch) const {
        std::vector<std::string> errors;
        for (const auto& [patch, constraints] : by_patch) {
            for (const auto& constraint : constraints) {
                std::visit([&](const auto& condition) {
                    using T = std::decay_t<decltype(condition)>;
                    if constexpr (std::is_same_v<T, Periodic>) {
                        const auto partner = by_patch.find(condition.partner_boundary);
                        if (partner == by_patch.end()) {
                            errors.push_back("periodic boundary '" + patch +
                                             "' references missing partner '" +
                                             condition.partner_boundary + "'");
                            return;
                        }
                        bool reciprocal = false;
                        for (const auto& pc : partner->second) {
                            if (pc.field != constraint.field) continue;
                            if (const auto* p = std::get_if<Periodic>(&pc.condition))
                                reciprocal = p->partner_boundary == patch;
                        }
                        if (!reciprocal)
                            errors.push_back("periodic boundary '" + patch +
                                             "' is not reciprocally paired with '" +
                                             condition.partner_boundary + "'");
                    } else if constexpr (std::is_same_v<T, Coupled>) {
                        if (by_patch.find(condition.target_boundary) == by_patch.end())
                            errors.push_back("coupled boundary '" + patch +
                                             "' references missing target '" +
                                             condition.target_boundary + "'");
                    }
                }, constraint.condition);
            }
        }
        return errors;
    }

private:
    std::vector<BoundaryRequirement> requirements_;
};

} // namespace cfdx::core
