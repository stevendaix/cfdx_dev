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
    PIMPLE
};

struct CouplingControls {
    double alpha_u = 0.7;
    double alpha_p = 0.3;
    int n_outer_correctors = 1;
    int n_pressure_correctors = 2;
};

inline void validate_coupling_controls(const CouplingControls& c)
{
    if (!(c.alpha_u > 0.0 && c.alpha_u <= 1.0) ||
        !(c.alpha_p > 0.0 && c.alpha_p <= 1.0) ||
        c.n_outer_correctors < 1 || c.n_pressure_correctors < 1)
        throw std::invalid_argument("pressure-velocity controls: invalid relaxation/corrector count");
}

inline double relaxed_value(double old_value, double computed_value, double alpha)
{
    if (alpha <= 0.0 || alpha > 1.0)
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
    if (d_f <= 0.0 || aP_owner <= 0.0 || aP_neighbour <= 0.0 || area < 0.0)
        throw std::invalid_argument("rhie_chow_face_flux: invalid geometric/momentum coefficient");
    const double d = 0.5 * (1.0 / aP_owner + 1.0 / aP_neighbour);
    const double correction = d * (pressure_owner - pressure_neighbour) / d_f * area;
    return interpolated_flux - correction;
}

inline double piso_correction_gain(double diagonal, double neighbor_sum)
{
    if (diagonal <= 0.0)
        throw std::invalid_argument("piso_correction_gain: diagonal must be positive");
    return 1.0 / std::max(diagonal - neighbor_sum, diagonal * 1e-12);
}

}  // namespace cfdx::physics
