#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"

#include <cctype>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

// The production incompressible driver currently advances velocity and pressure
// only. Until a coupled RANS outer loop is implemented, accepting a declared
// RANS model would silently execute laminar equations for a turbulent case.
inline std::string normalize_turbulence_model_key(const std::string& declared)
{
    std::string normalized;
    normalized.reserve(declared.size());
    for (const char raw_character : declared) {
        const auto character = static_cast<unsigned char>(raw_character);
        if (std::isalnum(character)) {
            normalized.push_back(static_cast<char>(std::tolower(character)));
        }
    }
    return normalized;
}

inline bool is_explicit_sst_model_key(const std::string& declared)
{
    const auto normalized = normalize_turbulence_model_key(declared);
    return normalized == "sst" || normalized == "komegasst";
}

inline bool is_explicit_kepsilon_model_key(const std::string& declared)
{
    const auto normalized = normalize_turbulence_model_key(declared);
    return normalized == "kepsilon" || normalized == "ke"
        || normalized == "realizablekepsilon" || normalized == "realizableke"
        || normalized == "realizable";
}

// General turbulence selection remains fail-closed. The only opt-in exceptions
// are the explicitly experimental SST production path and the experimental
// k-epsilon path; neither is a qualification claim and neither must silently
// turn an arbitrary RANS declaration on.
inline void validate_production_turbulence_selection(
    const cfdx::io::CaseSetup& setup,
    bool allow_experimental_sst = false,
    bool allow_experimental_kepsilon = false)
{
    const bool case_sst = is_explicit_sst_model_key(setup.turbulence_model);
    const bool case_ke = is_explicit_kepsilon_model_key(setup.turbulence_model);
    if (allow_experimental_sst && !case_sst && !case_ke)
        throw std::invalid_argument(
            "production solver: --experimental-sst requires case turbulence_model "
            "to explicitly select k-omega SST");
    if (allow_experimental_kepsilon && !case_ke && !case_sst)
        throw std::invalid_argument(
            "production solver: --experimental-kepsilon requires case turbulence_model "
            "to explicitly select k-epsilon");

    const auto validate = [allow_experimental_sst, allow_experimental_kepsilon](
        const std::string& declared, const std::string& source) {
        const std::string normalized = normalize_turbulence_model_key(declared);
        const bool unspecified =
            declared.find_first_not_of(" \t\r\n") == std::string::npos;
        if (unspecified || normalized == "laminar" ||
            normalized == "none" || normalized == "off" ||
            normalized == "disabled" || normalized == "incompressiblelaminar")
            return;
        if (allow_experimental_sst && is_explicit_sst_model_key(declared))
            return;
        if (allow_experimental_kepsilon && is_explicit_kepsilon_model_key(declared))
            return;
        throw std::invalid_argument(
            std::string("production solver: declared turbulence model '") +
            declared + "' from " + source +
            " is not coupled to the production momentum solver; only laminar flow "
            "is supported unless the case explicitly opts into the experimental SST "
            "or k-epsilon path");
    };

    const std::string physics =
        normalize_turbulence_model_key(setup.physics_model);
    if (physics.find("turbulent") != std::string::npos ||
        physics == "rans" || physics.find("rans") != std::string::npos) {
        if (!(allow_experimental_sst && case_sst) &&
            !(allow_experimental_kepsilon && case_ke))
            throw std::invalid_argument(
                "production solver: physics_model '" + setup.physics_model +
                "' declares a turbulent flow model that is not coupled to the "
                "production momentum solver; only laminar flow is supported without "
                "--experimental-sst or --experimental-kepsilon");
    }

    validate(setup.turbulence_model, "case setup");
    for (const auto& material : setup.materials)
        validate(material.turbulence_model, "material '" + material.name + "'");
}
}  // namespace cfdx::physics
