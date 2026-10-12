// Numerical Methods Maturity — Issue #461
// Common machine-checkable contract metadata for CFDX numerical methods.
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::core {

enum class NumericalMethodFamily {
    Gradient, Interpolation, Convection, Diffusion, Source, Temporal,
    TimeStep, Nonlinear, LinearSolver, Preconditioner, PressureVelocity,
    Conservation, Reconstruction, Initialization, Schur
};

enum class VerificationStatus { Planned, Implemented, Verified, Validated };
enum class ConservationContract {
    NotApplicable, LocalFaceConservative, GloballyConservative
};

struct NumericalMethodContract {
    std::string id;
    std::string name;
    NumericalMethodFamily family = NumericalMethodFamily::Convection;
    VerificationStatus status = VerificationStatus::Planned;
    ConservationContract conservation = ConservationContract::NotApplicable;
    bool bounded = false;
    bool monotone = false;
    bool linear_exact = false;
    int formal_spatial_order = 0;
    int formal_temporal_order = 0;
    std::string mathematical_formulation;
    std::string configuration_key;
    std::vector<std::string> verification_requirements;
};

inline const char* to_string(NumericalMethodFamily family) {
    switch (family) {
        case NumericalMethodFamily::Gradient: return "gradient";
        case NumericalMethodFamily::Interpolation: return "interpolation";
        case NumericalMethodFamily::Convection: return "convection";
        case NumericalMethodFamily::Diffusion: return "diffusion";
        case NumericalMethodFamily::Source: return "source";
        case NumericalMethodFamily::Temporal: return "temporal";
        case NumericalMethodFamily::TimeStep: return "time_step";
        case NumericalMethodFamily::Nonlinear: return "nonlinear";
        case NumericalMethodFamily::LinearSolver: return "linear_solver";
        case NumericalMethodFamily::Preconditioner: return "preconditioner";
        case NumericalMethodFamily::PressureVelocity: return "pressure_velocity";
        case NumericalMethodFamily::Conservation: return "conservation";
        case NumericalMethodFamily::Reconstruction: return "reconstruction";
        case NumericalMethodFamily::Initialization: return "initialization";
        case NumericalMethodFamily::Schur: return "schur";
    }
    return "unknown";
}

inline const char* to_string(VerificationStatus status) {
    switch (status) {
        case VerificationStatus::Planned: return "planned";
        case VerificationStatus::Implemented: return "implemented";
        case VerificationStatus::Verified: return "verified";
        case VerificationStatus::Validated: return "validated";
    }
    return "unknown";
}

inline const char* to_string(ConservationContract contract) {
    switch (contract) {
        case ConservationContract::NotApplicable: return "not_applicable";
        case ConservationContract::LocalFaceConservative: return "local_face_conservative";
        case ConservationContract::GloballyConservative: return "globally_conservative";
    }
    return "unknown";
}

inline void validate_numerical_method_contract(const NumericalMethodContract& c) {
    if (c.id.empty() || c.name.empty())
        throw std::invalid_argument("numerical method contract requires id and name");
    if (c.id.find(' ') != std::string::npos)
        throw std::invalid_argument("numerical method contract id must not contain spaces");
    if (c.formal_spatial_order < 0 || c.formal_temporal_order < 0)
        throw std::invalid_argument("numerical method formal orders must be non-negative");
    if (c.status >= VerificationStatus::Implemented &&
        (c.mathematical_formulation.empty() || c.configuration_key.empty()))
        throw std::invalid_argument("implemented numerical method requires formulation and configuration key");
    if (c.status == VerificationStatus::Validated && c.verification_requirements.empty())
        throw std::invalid_argument("validated numerical method requires verification requirements");
    if (c.bounded && !c.monotone && c.family == NumericalMethodFamily::Convection &&
        c.conservation == ConservationContract::NotApplicable)
        throw std::invalid_argument(
            "bounded convection contract must declare a conservation contract");
}

} // namespace cfdx::core
