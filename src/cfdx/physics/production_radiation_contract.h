#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"

#include <stdexcept>
#include <string>

namespace cfdx::physics {

// Fail closed until a radiation model is actually orchestrated by the
// production driver. The case schema already accepts radiation declarations;
// silently running without their source/flux contribution is physically wrong.
inline void validate_production_radiation_selection(
    const cfdx::io::CaseSetup& setup)
{
    const std::string& model = setup.radiation_model;
    if (model.empty() || model == "none" || model == "off" ||
        model == "disabled") {
        return;
    }

    throw std::invalid_argument(
        "production solver: declared radiation_model '" + model +
        "' is not coupled to the production energy equation yet; "
        "supported declarations are none/off only until the production "
        "radiation stage is integrated (S2S, DO/DOM and solar are not "
        "silently treated as radiation-free)");
}

}  // namespace cfdx::physics
