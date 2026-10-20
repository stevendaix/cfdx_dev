#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"

#include <cctype>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

// Fail closed until a radiation model is actually orchestrated by the
// production driver. The case schema already accepts radiation declarations
// (fiche methode docs/audits/11-rayonnement.md: Rosseland, P1, S2S, DOM);
// silently running a declared case without its radiative source/flux
// contribution is physically wrong, so any declared model is rejected
// instead of being treated as radiation-free.
inline std::string normalize_radiation_model_key(const std::string& declared)
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

inline void validate_production_radiation_selection(
    const cfdx::io::CaseSetup& setup)
{
    const std::string& model = setup.radiation_model;
    // Whitespace-only declaration = not declared (legacy case files may
    // omit the block); same contract as the turbulence selection guard.
    const bool unspecified =
        model.find_first_not_of(" \t\r\n") == std::string::npos;
    if (unspecified) {
        return;
    }
    const std::string normalized = normalize_radiation_model_key(model);
    if (normalized == "none" || normalized == "off" ||
        normalized == "disabled") {
        return;
    }

    throw std::invalid_argument(
        "production solver: declared radiation_model '" + model +
        "' is not coupled to the production energy equation yet; "
        "supported declarations are none/off only until the production "
        "radiation stage is integrated (S2S, DO/DOM, P1, Rosseland and "
        "solar are not silently treated as radiation-free)");
}

}  // namespace cfdx::physics
