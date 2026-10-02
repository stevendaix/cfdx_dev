// Numerical Methods Maturity — Issue #461 (N1)
// Explicit, auditable numerical-scheme selection from case/numerics config.
//
// The registry in numerical_method_registry.h declares, for every method, the
// configuration key under which a case selects it. That key was previously
// metadata only: nothing resolved it, so a case could not be proven to have
// selected a scheme and a missing key could silently fall back to an implicit
// default. This header closes that N1 gap.
//
// Selection is explicit and total:
//   - a family resolves to exactly one registered method via its key;
//   - an empty, unknown, ambiguous or wrong-family key throws;
//   - the resolved selection carries the method id and verification status so
//     it can be recorded verbatim in a validation report.
#pragma once

#include "cfdx/core/numerics/numerical_method_contract.h"
#include "cfdx/core/numerics/numerical_method_registry.h"

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::core {

// A resolved scheme selection suitable for case configuration and reports.
struct NumericalSchemeSelection {
    NumericalMethodFamily family = NumericalMethodFamily::Convection;
    std::string configuration_key;
    std::string method_id;
    std::string method_name;
    VerificationStatus status = VerificationStatus::Planned;
    bool resolved = false;
};

inline std::optional<NumericalMethodContract> find_method_by_id(const std::string& id) {
    for (const auto& method : numerical_method_registry()) {
        if (method.id == id) return method;
    }
    return std::nullopt;
}

// Resolve a configuration key to its unique registered method. A key that maps
// to more than one method is a registry defect and throws rather than picking
// an arbitrary winner.
inline std::optional<NumericalMethodContract> find_method_by_configuration_key(
    const std::string& configuration_key) {
    std::optional<NumericalMethodContract> found;
    for (const auto& method : numerical_method_registry()) {
        if (method.configuration_key == configuration_key) {
            if (found) {
                throw std::invalid_argument(
                    "ambiguous numerical method configuration key: " + configuration_key);
            }
            found = method;
        }
    }
    return found;
}

// Resolve an explicit scheme selection. Throws when the key is empty, unknown
// or does not belong to the requested family, so a case can never silently fall
// back to an implicit scheme.
inline NumericalMethodContract select_numerical_method(
    NumericalMethodFamily family, const std::string& configuration_key) {
    if (configuration_key.empty()) {
        throw std::invalid_argument(
            "numerical method selection requires an explicit configuration key");
    }
    const auto method = find_method_by_configuration_key(configuration_key);
    if (!method) {
        throw std::invalid_argument(
            "unknown numerical method configuration key: " + configuration_key);
    }
    if (method->family != family) {
        throw std::invalid_argument(
            "numerical method configuration key '" + configuration_key +
            "' does not belong to family '" + std::string(to_string(family)) + "'");
    }
    return *method;
}

inline NumericalSchemeSelection make_scheme_selection(
    NumericalMethodFamily family, const std::string& configuration_key) {
    const auto method = select_numerical_method(family, configuration_key);
    NumericalSchemeSelection selection;
    selection.family = family;
    selection.configuration_key = configuration_key;
    selection.method_id = method.id;
    selection.method_name = method.name;
    selection.status = method.status;
    selection.resolved = true;
    return selection;
}

// Report-ready representation. Recording this string in a validation report
// makes the selected scheme and its verification status auditable without
// re-deriving them from the case file.
inline std::string format_scheme_selection(const NumericalSchemeSelection& selection) {
    return std::string("scheme[") + to_string(selection.family) + "]=" +
           selection.configuration_key + "->" + selection.method_id +
           "(" + to_string(selection.status) + ")";
}

// Every configuration key a case may use to select a method in this family.
inline std::vector<std::string> configuration_keys_for_family(
    NumericalMethodFamily family) {
    std::vector<std::string> keys;
    for (const auto& method : numerical_method_registry()) {
        if (method.family == family) keys.push_back(method.configuration_key);
    }
    return keys;
}

// One explicit selection request of a case numerics configuration.
struct CaseNumericsEntry {
    NumericalMethodFamily family = NumericalMethodFamily::Convection;
    std::string configuration_key;
};

// Result of resolving a whole case numerics configuration. Every request is
// either resolved (with auditable fields) or recorded as an error — there is
// no silent fallback and no throw on a malformed case: the report is what the
// case-loading and validation pipeline consumes.
struct CaseNumericsReport {
    std::string case_name;
    std::vector<NumericalSchemeSelection> resolved;
    std::vector<std::string> errors;

    bool valid() const { return errors.empty(); }
};

inline CaseNumericsReport resolve_case_numerics(
    const std::string& case_name,
    const std::vector<CaseNumericsEntry>& config)
{
    CaseNumericsReport report;
    report.case_name = case_name;
    for (const auto& entry : config) {
        try {
            report.resolved.push_back(
                make_scheme_selection(entry.family, entry.configuration_key));
        } catch (const std::invalid_argument& e) {
            report.errors.push_back(e.what());
        }
    }
    return report;
}

// Deterministic plain-text validation report: one error line per rejected
// entry and one resolved-scheme line per accepted entry, in input order.
inline std::string format_numerics_report(const CaseNumericsReport& report)
{
    std::string out = "case=" + report.case_name + "\n";
    for (const auto& e : report.errors)
        out += std::string("error: ") + e + "\n";
    for (const auto& s : report.resolved)
        out += format_scheme_selection(s) + "\n";
    return out;
}

} // namespace cfdx::core
