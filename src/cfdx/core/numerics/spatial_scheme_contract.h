// N1 — common spatial discretisation contracts.
//
// This header separates the identity of a spatial scheme from its concrete
// implementation. It is deliberately small: N1 defines the contract and
// compatibility rules; N2-N4 provide the actual gradient/reconstruction,
// diffusion and convection implementations.
#pragma once

#include "cfdx/core/numerics/numerical_method_contract.h"

#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::core {

enum class SpatialSchemeRole {
    Gradient,
    Reconstruction,
    Interpolation,
    Convection,
    Diffusion,
    Flux
};

inline const char* to_string(SpatialSchemeRole role) {
    switch (role) {
        case SpatialSchemeRole::Gradient: return "gradient";
        case SpatialSchemeRole::Reconstruction: return "reconstruction";
        case SpatialSchemeRole::Interpolation: return "interpolation";
        case SpatialSchemeRole::Convection: return "convection";
        case SpatialSchemeRole::Diffusion: return "diffusion";
        case SpatialSchemeRole::Flux: return "flux";
    }
    return "unknown";
}

// A role-level contract used by all spatial discretisations. The registry key
// remains the single selection identity; these fields describe what the
// selected operator is allowed to claim and are independently testable.
struct SpatialSchemeContract {
    SpatialSchemeRole role = SpatialSchemeRole::Gradient;
    std::string configuration_key;
    ConservationContract conservation = ConservationContract::NotApplicable;
    bool bounded = false;
    bool monotone = false;
    bool linear_exact = false;
    int formal_order = 0;
    std::string formulation;
    std::vector<std::string> verification_requirements;
};

inline void validate_spatial_scheme_contract(const SpatialSchemeContract& c) {
    if (c.configuration_key.empty())
        throw std::invalid_argument("spatial scheme contract requires configuration key");
    if (c.formulation.empty())
        throw std::invalid_argument("spatial scheme contract requires formulation");
    if (c.formal_order < 0)
        throw std::invalid_argument("spatial scheme formal order must be non-negative");
    if (c.bounded && c.role == SpatialSchemeRole::Convection &&
        c.conservation == ConservationContract::NotApplicable) {
        throw std::invalid_argument(
            "bounded convection scheme must declare a conservation contract");
    }
}

// Canonical case-level spatial selection. Concrete implementations resolve the
// key through numerical_method_selection.h; this type prevents callers from
// inventing a second scheme-selection vocabulary for spatial operators.
struct SpatialSchemeSelection {
    SpatialSchemeRole role = SpatialSchemeRole::Gradient;
    NumericalSchemeSelection method;
};

inline SpatialSchemeSelection make_spatial_scheme_selection(
    SpatialSchemeRole role,
    const NumericalMethodFamily family,
    const std::string& configuration_key)
{
    const auto method = make_scheme_selection(family, configuration_key);
    const bool compatible =
        (role == SpatialSchemeRole::Gradient && family == NumericalMethodFamily::Gradient) ||
        (role == SpatialSchemeRole::Reconstruction && family == NumericalMethodFamily::Reconstruction) ||
        (role == SpatialSchemeRole::Interpolation && family == NumericalMethodFamily::Interpolation) ||
        (role == SpatialSchemeRole::Convection && family == NumericalMethodFamily::Convection) ||
        (role == SpatialSchemeRole::Diffusion && family == NumericalMethodFamily::Diffusion) ||
        (role == SpatialSchemeRole::Flux && family == NumericalMethodFamily::Conservation);
    if (!compatible) {
        throw std::invalid_argument(
            "spatial scheme role '" + std::string(to_string(role)) +
            "' is incompatible with numerical family '" +
            std::string(to_string(family)) + "'");
    }
    return SpatialSchemeSelection{role, method};
}

} // namespace cfdx::core
