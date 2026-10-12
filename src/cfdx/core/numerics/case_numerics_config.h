#pragma once

#include "cfdx/core/numerics/numerical_method_selection.h"

#include <algorithm>
#include <string>
#include <vector>

namespace cfdx::core {

// Explicit case-level numerical configuration. The case declares the
// selections that are applicable to it and, separately, the families that
// are mandatory for that case. This avoids both hidden defaults and a
// hard-coded global list of families (e.g. a scalar diffusion case need not
// select a convection scheme).
struct CaseNumericsConfig {
    std::vector<CaseNumericsEntry> entries;
    std::vector<NumericalMethodFamily> required_families;
};

inline CaseNumericsReport resolve_case_numerics(
    const std::string& case_name,
    const CaseNumericsConfig& config)
{
    CaseNumericsReport report = resolve_case_numerics(case_name, config.entries);

    // A numerical configuration with no explicit selections is never valid.
    if (config.entries.empty()) {
        report.errors.push_back(
            "case numerical configuration must contain at least one explicit selection");
        return report;
    }

    // A family can have at most one selected method. Reject ambiguity before
    // solver setup instead of allowing a later consumer to choose a winner.
    for (const auto family : config.required_families) {
        std::size_t count = 0;
        for (const auto& entry : config.entries)
            if (entry.family == family) ++count;

        if (count == 0) {
            report.errors.push_back(
                "required numerical method family is not selected: " +
                std::string(to_string(family)));
        } else if (count > 1) {
            report.errors.push_back(
                "numerical method family has multiple selections: " +
                std::string(to_string(family)));
        }
    }

    for (const auto& entry : config.entries) {
        std::size_t count = 0;
        for (const auto& other : config.entries)
            if (other.family == entry.family) ++count;
        if (count > 1) {
            const std::string error =
                "numerical method family has multiple selections: " +
                std::string(to_string(entry.family));
            if (std::find(report.errors.begin(), report.errors.end(), error) ==
                report.errors.end()) {
                report.errors.push_back(error);
            }
        }
    }

    return report;
}

} // namespace cfdx::core
