#pragma once

#include "cfdx/core/field/field.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

struct TurbulenceCoefficients {
    double C_mu = 0.09;
    double C1 = 1.44;
    double C2 = 1.92;
    double sigma_k = 1.0;
    double sigma_epsilon = 1.3;
    double beta_star = 0.09;
    double beta = 0.075;
    double alpha = 0.52;
};

inline double turbulent_kinematic_viscosity_kepsilon(double k, double epsilon,
                                                      double C_mu = 0.09)
{
    if (k < 0.0 || epsilon <= 0.0 || C_mu <= 0.0)
        throw std::invalid_argument("k-epsilon: expected k>=0, epsilon>0 and C_mu>0");
    return C_mu * k * k / epsilon;
}

inline double turbulent_kinematic_viscosity_komega_sst(
    double k, double omega, double F2, double strain_rate, double a1 = 0.31)
{
    if (k < 0.0 || omega <= 0.0 || F2 < 0.0 || F2 > 1.0 ||
        strain_rate < 0.0 || !std::isfinite(strain_rate) || a1 <= 0.0)
        throw std::invalid_argument("k-omega SST: invalid k, omega, F2 or a1");
    return a1 * k / std::max(a1 * omega, strain_rate * F2);
}

// Backward-compatible overload retained for existing unit tests/callers.
// A zero strain-rate limiter recovers the pre-existing algebraic form.
inline double turbulent_kinematic_viscosity_komega_sst(
    double k, double omega, double F2, double a1 = 0.31)
{
    return turbulent_kinematic_viscosity_komega_sst(k, omega, F2, 0.0, a1);
}

inline double strain_rate_magnitude(const cfdx::core::Vec3& s)
{
    return std::sqrt(2.0 * (s.x*s.x + s.y*s.y + s.z*s.z));
}

inline double smagorinsky_eddy_viscosity(double delta, double strain,
                                         double Cs = 0.17)
{
    if (delta <= 0.0 || strain < 0.0 || Cs < 0.0)
        throw std::invalid_argument("Smagorinsky: invalid delta, strain or Cs");
    return (Cs * delta) * (Cs * delta) * strain;
}


// Algebraic DES helper: this is an eddy-viscosity length-scale model,
// not a full hybrid RANS-LES transport-equation DES implementation.
inline double des_eddy_viscosity(double delta, double wall_distance,
                                 double strain, double Cs = 0.17,
                                 double Cdes = 0.65)
{
    if (delta <= 0.0 || wall_distance <= 0.0 || strain < 0.0 ||
        Cs < 0.0 || Cdes <= 0.0)
        throw std::invalid_argument("DES: invalid model parameters");
    const double length = std::min(wall_distance, Cdes * delta);
    return (Cs * length) * (Cs * length) * strain;
}

inline double turbulence_production(const cfdx::core::Vec3& velocity_gradient_times_symmetric_part,
                                    double nu_t)
{
    if (nu_t < 0.0)
        throw std::invalid_argument("turbulence production: nu_t must be non-negative");
    return 2.0 * nu_t * (velocity_gradient_times_symmetric_part.x +
                         velocity_gradient_times_symmetric_part.y +
                         velocity_gradient_times_symmetric_part.z);
}

inline void validate_turbulence_field(const cfdx::core::Field<double, cfdx::core::Location::CELL>& f,
                                      std::size_t n, const char* name)
{
    if (f.size() != n || f.dimension() != 1)
        throw std::runtime_error(std::string(name) + ": invalid cell scalar field");
}

}  // namespace cfdx::physics
