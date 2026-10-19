#pragma once

#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/physics/turbulence_transport.h"

#include <cmath>
#include <cctype>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// Inlet turbulence intensity is a fraction (0 < I < 1), not a percentage.
// The length scale and velocity use the case's consistent length/velocity units.
struct SSTInletTurbulenceState {
    double k = 0.0;      // m^2/s^2
    double omega = 0.0;  // 1/s
};

// Derive the common isotropic inlet estimate k = 3/2 (|U| I)^2 and
// omega = sqrt(k) / (beta_star^(1/4) L). This is an inlet estimate, not a
// wall function or a claim of physical qualification.
inline SSTInletTurbulenceState derive_sst_inlet_turbulence_state(
    const cfdx::io::BoundarySpec& boundary,
    const TurbulenceTransportControls& controls)
{
    if (boundary.type != cfdx::io::BCType::INLET)
        throw std::invalid_argument("SST inlet turbulence state requires an inlet boundary");
    if (controls.model != TurbulenceModel::SST)
        throw std::invalid_argument("SST inlet turbulence state requires the SST model");
    validate_turbulence_controls(controls);
    if (!std::isfinite(controls.beta_star) || controls.beta_star <= 0.0)
        throw std::invalid_argument("SST beta_star must be finite and positive");
    if (!std::isfinite(boundary.turbulence_intensity) ||
        boundary.turbulence_intensity <= 0.0 ||
        boundary.turbulence_intensity >= 1.0)
        throw std::invalid_argument("SST inlet turbulence intensity must be a fraction in (0, 1)");
    if (!std::isfinite(boundary.turbulence_length_scale) ||
        boundary.turbulence_length_scale <= 0.0)
        throw std::invalid_argument("SST inlet turbulence length scale must be finite and positive");

    double speed = boundary.velocity_magnitude;
    if (!boundary.velocity_vector.empty()) {
        if (boundary.velocity_vector.size() != 3)
            throw std::invalid_argument("SST inlet velocity vector must have three components");
        double speed_squared = 0.0;
        for (const double component : boundary.velocity_vector) {
            if (!std::isfinite(component))
                throw std::invalid_argument("SST inlet velocity vector must be finite");
            speed_squared += component * component;
        }
        speed = std::sqrt(speed_squared);
    }
    if (!std::isfinite(speed) || speed <= 0.0)
        throw std::invalid_argument("SST inlet velocity magnitude must be finite and positive");

    const double speed_intensity = speed * boundary.turbulence_intensity;
    const double k = 1.5 * speed_intensity * speed_intensity;
    const double omega = std::sqrt(k) /
        (std::pow(controls.beta_star, 0.25) * boundary.turbulence_length_scale);
    if (!std::isfinite(k) || k <= 0.0 || !std::isfinite(omega) || omega <= 0.0)
        throw std::domain_error("SST inlet turbulence estimate is non-finite or non-positive");
    return {k, omega};
}

// Validate only the currently representable SST initialization and inlet
// prerequisites. CaseSetup does not yet express explicit k/omega wall and outlet
// scalar boundary contracts, so this function must not be treated as permission
// to enable SST in the production selector.
inline void validate_sst_case_prerequisites(
    const cfdx::io::CaseSetup& setup,
    const TurbulenceTransportControls& controls)
{
    if (controls.model != TurbulenceModel::SST)
        throw std::invalid_argument("SST case prerequisites require the SST model");
    validate_turbulence_controls(controls);

    std::string model;
    for (const unsigned char ch : setup.turbulence_model) {
        if (std::isalnum(ch))
            model.push_back(static_cast<char>(std::tolower(ch)));
    }
    if (model != "sst" && model != "komegasst")
        throw std::invalid_argument("case must explicitly select k-omega SST");

    if (!setup.has_initial_condition)
        throw std::invalid_argument("SST case must explicitly declare initial conditions");
    const auto k_it = setup.initial_condition.scalar_fields.find("k");
    const auto omega_it = setup.initial_condition.scalar_fields.find("omega");
    if (k_it == setup.initial_condition.scalar_fields.end() ||
        omega_it == setup.initial_condition.scalar_fields.end())
        throw std::invalid_argument("SST initial conditions must explicitly provide scalar fields 'k' and 'omega'");

    if (!std::isfinite(k_it->second) || k_it->second < controls.k_min)
        throw std::invalid_argument("SST initial k must be finite and at least k_min");
    if (!std::isfinite(omega_it->second) || omega_it->second < controls.omega_min)
        throw std::invalid_argument("SST initial omega must be finite and at least omega_min");

    std::map<std::string, bool> patch_names;
    for (const auto& boundary : setup.boundary_conditions) {
        if (boundary.patch_name.empty())
            throw std::invalid_argument("SST boundary patch names must not be empty");
        if (!patch_names.emplace(boundary.patch_name, true).second)
            throw std::invalid_argument("SST case contains duplicate boundary patch '" +
                                        boundary.patch_name + "'");
        if (boundary.type == cfdx::io::BCType::UNKNOWN ||
            boundary.type == cfdx::io::BCType::INTERNAL)
            throw std::invalid_argument("SST case contains an unsupported boundary type on patch '" +
                                        boundary.patch_name + "'");
        if (boundary.type == cfdx::io::BCType::INLET) {
            // Reuse the exact same validation and derivation used by the solver
            // when it constructs inlet scalar values.
            (void)derive_sst_inlet_turbulence_state(boundary, controls);
        }
    }
}

} // namespace cfdx::physics
