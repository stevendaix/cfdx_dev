#pragma once

// NACA test-driver selection, NOT a replacement for the production COUPLED
// solver. The six existing Schur approximations must remain explicitly named
// in every experimental NACA COUPLED run to avoid silent fallback to PCD.
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include <stdexcept>
#include <string>

namespace cfdx::validation {

struct NacaCouplingSelection {
    cfdx::physics::PressureVelocityAlgorithm algorithm;
    cfdx::physics::CoupledSchurModel schur_model;
};

inline NacaCouplingSelection parse_naca_coupling_selection(
    const char* algorithm, const char* coupled_schur)
{
    using namespace cfdx::physics;
    const std::string method = algorithm ? algorithm : "SIMPLE";
    PressureVelocityAlgorithm selected;
    if (method == "SIMPLE") selected = PressureVelocityAlgorithm::SIMPLE;
    else if (method == "SIMPLEC") selected = PressureVelocityAlgorithm::SIMPLEC;
    else if (method == "COUPLED") selected = PressureVelocityAlgorithm::COUPLED;
    else throw std::invalid_argument(
        "CFDX_NACA_ALGORITHM expects SIMPLE, SIMPLEC or COUPLED");

    CoupledSchurModel schur = CoupledSchurModel::PCD;
    if (selected == PressureVelocityAlgorithm::COUPLED) {
        if (!coupled_schur || !parse_coupled_schur_model(coupled_schur, schur))
            throw std::invalid_argument(
                "CFDX_NACA_ALGORITHM=COUPLED requires explicit "
                "CFDX_NACA_COUPLED_SCHUR=block_local|pcd|lsc|bfbt|simple|simplec");
    } else if (coupled_schur) {
        throw std::invalid_argument(
            "CFDX_NACA_COUPLED_SCHUR is only available with COUPLED");
    }
    return {selected, schur};
}

inline const char* naca_algorithm_name(
    cfdx::physics::PressureVelocityAlgorithm algorithm) noexcept
{
    using cfdx::physics::PressureVelocityAlgorithm;
    switch (algorithm) {
        case PressureVelocityAlgorithm::SIMPLE: return "SIMPLE";
        case PressureVelocityAlgorithm::SIMPLEC: return "SIMPLEC";
        case PressureVelocityAlgorithm::COUPLED: return "COUPLED";
        default: return "UNSUPPORTED";
    }
}

} // namespace cfdx::validation