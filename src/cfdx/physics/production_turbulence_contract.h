#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"

#include <cctype>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

// The production incompressible driver currently advances velocity and pressure
// only. Until a coupled RANS outer loop is implemented, accepting a declared
// RANS model would silently execute laminar equations for a turbulent case.
inline void validate_production_turbulence_selection(
    const cfdx::io::CaseSetup& setup)
{
    const auto validate = [](const std::string& declared, const char* source) {
        std::string normalized;
        normalized.reserve(declared.size());
        for (const unsigned char character : declared) {
            if (std::isalnum(character)) {
                normalized.push_back(
                    static_cast<char>(std::tolower(character)));
            }
        }
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

    validate(setup.turbulence_model, "case setup");
    for (const auto& material : setup.materials) {
        validate(material.turbulence_model,
                 ("material '" + material.name + "'").c_str());
    }
}

}  // namespace cfdx::physics
