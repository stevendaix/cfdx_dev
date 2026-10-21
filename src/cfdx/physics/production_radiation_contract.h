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
// instead of being treated as radiation-free. S2S is the first coupled
// exception (issue #885 stage 2): the surface-to-surface radiosity stage
// feeds the production energy equation, so it is accepted only when the
// conjugate energy model carries it.
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

// Routing predicate of the production driver (issue #885 stage 2): true
// when the case declares S2S radiation, whatever the energy model is.
// The driver selects the coupled S2S->energy stage with it; the stage
// itself re-checks the conjugate energy model fail-closed.
inline bool production_s2s_radiation_requested(const cfdx::io::CaseSetup& setup)
{
    return normalize_radiation_model_key(setup.radiation_model) == "s2s";
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
    if (normalized == "s2s") {
        if (setup.energy_model == "conjugate" ||
            setup.energy_model == "buoyant") {
            return;
        }
        throw std::invalid_argument(
            "production solver: declared radiation_model '" + model +
            "' requires energy_model 'conjugate' or 'buoyant' (declared: '" +
            setup.energy_model + "'); the S2S radiosity stage couples into "
            "the production energy equation, so a radiation declaration "
            "without a coupled energy model is rejected instead of "
            "silently running radiation-free");
    }

    throw std::invalid_argument(
        "production solver: declared radiation_model '" + model +
        "' is not coupled to the production energy equation yet; "
        "supported declarations are none/off and S2S with energy_model "
        "'conjugate' or 'buoyant' (DO/DOM, P1, Rosseland and solar are "
        "not silently treated as radiation-free)");
}

}  // namespace cfdx::physics
