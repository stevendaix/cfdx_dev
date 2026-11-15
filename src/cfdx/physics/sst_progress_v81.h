#pragma once
// V81 forensic-only, not a solver or a convergence criterion.
// All quantities retain the units and thresholds of the production gates.
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {
struct SSTProgressDistanceV81 {
    double k_row_factor = 0.0;
    double omega_row_factor = 0.0;
    double dk_factor = 0.0;
    double domega_factor = 0.0;
    double worst_factor = 0.0;
    bool all_gates_met = false;
};
inline SSTProgressDistanceV81 sst_progress_distance_v81(
    double k_row, double omega_row, double dk, double domega,
    double tolerance)
{
    if (!std::isfinite(k_row) || !std::isfinite(omega_row) ||
        !std::isfinite(dk) || !std::isfinite(domega) ||
        !std::isfinite(tolerance) || !(tolerance > 0.) ||
        k_row < 0. || omega_row < 0. || dk < 0. || domega < 0.)
        throw std::invalid_argument("V81 progress: nonfinite or negative gate input");
    SSTProgressDistanceV81 out;
    out.k_row_factor=k_row/tolerance;
    out.omega_row_factor=omega_row/tolerance;
    out.dk_factor=dk/tolerance;
    out.domega_factor=domega/tolerance;
    out.worst_factor=std::max({out.k_row_factor,out.omega_row_factor,
                               out.dk_factor,out.domega_factor});
    out.all_gates_met=out.worst_factor<=1.;
    if (!std::isfinite(out.worst_factor))
        throw std::overflow_error("V81 progress: gate ratio overflow");
    return out;
}
// A dimensional CHECK of the viscous omega wall value, never a substitute
// for actual OpenFOAM omegaWallFunction owner-row or G matrix treatment.
inline double sst_viscous_wall_omega_v81(
    double nu, double beta1, double wall_distance, double factor)
{
    if (!std::isfinite(nu) || !std::isfinite(beta1) ||
        !std::isfinite(wall_distance) || !std::isfinite(factor) ||
        !(nu > 0.) || !(beta1 > 0.) || !(wall_distance > 0.) ||
        !(factor > 0.))
        throw std::invalid_argument("V81 viscous omega: invalid physical input");
    const double value=(factor*nu/beta1)/wall_distance/wall_distance;
    if (!std::isfinite(value) || !(value > 0.))
        throw std::overflow_error("V81 viscous omega: overflow");
    return value;
}
}
