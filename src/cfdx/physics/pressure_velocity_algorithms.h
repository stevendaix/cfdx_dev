#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/flux.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

enum class PressureVelocityAlgorithm {
    SIMPLE,
    SIMPLEC,
    PISO,
    PIMPLE,
    // Fractional-step projection: momentum predictor followed by one or
    // more pressure projections. This is intentionally distinct from PISO:
    // each projection operates on the current conservative face flux.
    FRACTIONAL_STEP,
    // Fully coupled pressure-based solve. Momentum and continuity are
    // assembled in one block system with an implicit Rhie-Chow pressure block.
    COUPLED
};

struct CouplingControls {
    double alpha_u = 0.7;
    double alpha_p = 0.3;
    int n_outer_correctors = 1;
    int n_pressure_correctors = 2;
    int n_fractional_steps = 2;
    // Coupled linear solve controls. The coupled matrix is indefinite but the
    // Rhie-Chow pressure block makes the collocated formulation nonsingular
    // after pressure gauge fixing.
    std::size_t coupled_max_iterations = 2000;
    double coupled_linear_tolerance = 1e-10;
};

inline void validate_coupling_controls(const CouplingControls& c)
{
    if (!(c.alpha_u > 0.0 && c.alpha_u <= 1.0) ||
        !(c.alpha_p > 0.0 && c.alpha_p <= 1.0) ||
        c.n_outer_correctors < 1 || c.n_pressure_correctors < 1 ||
        c.n_fractional_steps < 1 || c.coupled_max_iterations == 0 ||
        !(c.coupled_linear_tolerance > 0.0))
        throw std::invalid_argument("pressure-velocity controls: invalid relaxation/corrector count");
}

inline double relaxed_value(double old_value, double computed_value, double alpha)
{
    if (!std::isfinite(alpha) || !std::isfinite(old_value) || !std::isfinite(computed_value) ||
        alpha <= 0.0 || alpha > 1.0)
        throw std::invalid_argument("relaxed_value: alpha must be in (0,1]");
    return old_value + alpha * (computed_value - old_value);
}

inline double rhie_chow_face_flux(double interpolated_flux,
                                  double pressure_owner,
                                  double pressure_neighbour,
                                  double d_f,
                                  double aP_owner,
                                  double aP_neighbour,
                                  double area)
{
    if (!std::isfinite(d_f) || !std::isfinite(aP_owner) ||
        !std::isfinite(aP_neighbour) || !std::isfinite(area) ||
        d_f <= 0.0 || aP_owner <= 0.0 || aP_neighbour <= 0.0 || area < 0.0)
        throw std::invalid_argument("rhie_chow_face_flux: invalid geometric/momentum coefficient");
    const double d = 0.5 * (1.0 / aP_owner + 1.0 / aP_neighbour);
    const double correction = d * (pressure_owner - pressure_neighbour) / d_f * area;
    return interpolated_flux - correction;
}

inline double piso_correction_gain(double diagonal, double neighbor_sum)
{
    if (!std::isfinite(diagonal) || !std::isfinite(neighbor_sum) || diagonal <= 0.0)
        throw std::invalid_argument("piso_correction_gain: diagonal must be positive");
    return 1.0 / std::max(diagonal - neighbor_sum, diagonal * 1e-12);
}


// The fractional-step pressure projection is deliberately explicit about its
// contract. It is a projection of the conservative face flux, not an alias for
// PISO. Keeping this helper here also makes the algorithm choice visible to
// callers without embedding policy in the finite-volume transport layer.
inline bool is_fractional_step_algorithm(PressureVelocityAlgorithm a)
{
    return a == PressureVelocityAlgorithm::FRACTIONAL_STEP;
}

inline bool is_coupled_algorithm(PressureVelocityAlgorithm a)
{
    return a == PressureVelocityAlgorithm::COUPLED;
}

}  // namespace cfdx::physics
