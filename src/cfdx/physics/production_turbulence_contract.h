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
    for (const unsigned char character : declared) {
        if (std::isalnum(character)) {
            normalized.push_back(static_cast<char>(std::tolower(character)));
        }
    }
    return normalized;
}

inline void validate_production_turbulence_selection(
    const cfdx::io::CaseSetup& setup)
{
    const auto validate = [](const std::string& declared, const std::string& source) {
        const std::string normalized = normalize_turbulence_model_key(declared);
        if (normalized.empty() || normalized == "laminar" ||
            normalized == "none" || normalized == "off" ||
            normalized == "disabled" || normalized == "incompressiblelaminar") {
            return;
        }
        throw std::invalid_argument(
            std::string("production solver: declared turbulence model '") +
            declared + "' from " + source +
            " is not coupled to the production momentum solver yet; "
            "only laminar flow is currently supported by this executable");
    };

    const std::string physics =
        normalize_turbulence_model_key(setup.physics_model);
    if ((physics.find("turbulent") != std::string::npos ||
         physics == "rans" || physics.find("rans") != std::string::npos) &&
        normalize_turbulence_model_key(setup.turbulence_model).empty()) {
        throw std::invalid_argument(
            "production solver: physics_model '" + setup.physics_model +
            "' declares turbulence but no turbulence_model is specified; "
            "the production solver currently supports laminar flow only");
    }

    validate(setup.turbulence_model, "case setup");
    for (const auto& material : setup.materials) {
        validate(material.turbulence_model,
                 "material '" + material.name + "'");
    }
}
}  // namespace cfdx::physics
