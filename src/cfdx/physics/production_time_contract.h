#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/physics/production_turbulence_contract.h"

#include <stdexcept>
#include <string>

namespace cfdx::physics {

inline void validate_production_steady_time_selection(
    const cfdx::io::CaseSetup& setup, const std::string& executable)
{
    const auto mode = normalize_turbulence_model_key(setup.solver_mode);
    const auto temporal = normalize_turbulence_model_key(setup.numerics.transient_scheme);
    const auto declared = [](const std::string& value) {
        return value.find_first_not_of(" \t\r\n") != std::string::npos;
    };
    bool explicit_temporal = false;
    for (const auto& entry : setup.numerical_config.entries)
        explicit_temporal = explicit_temporal ||
            entry.family == cfdx::core::NumericalMethodFamily::Temporal;
    for (const auto& entry : setup.numerical_report.resolved)
        explicit_temporal = explicit_temporal ||
            entry.family == cfdx::core::NumericalMethodFamily::Temporal;
    if (setup.transient || explicit_temporal ||
        (declared(setup.solver_mode) && mode != "steady") ||
        (declared(setup.numerics.transient_scheme) && temporal != "steady"))
        throw std::invalid_argument(
            executable + ": transient evolution is not implemented by this steady executable; "
            "the declared time integration cannot be silently ignored");
}

}  // namespace cfdx::physics
