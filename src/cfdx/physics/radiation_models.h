#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

enum class RadiationApproximation { Rosseland, P1, DOM };

struct RadiationModelSelector {
    double optical_thick = 3.0;
    double optical_thin = 0.1;

    RadiationApproximation select(double absorption, double scattering, double length) const {
        if (!std::isfinite(absorption) || !std::isfinite(scattering) ||
            !std::isfinite(length) || !std::isfinite(optical_thick) ||
            !std::isfinite(optical_thin) ||
            absorption < 0.0 || scattering < 0.0 || length <= 0.0 ||
            optical_thick <= optical_thin || optical_thin < 0.0) {
            throw std::invalid_argument("invalid radiation regime");
        }
        const double tau = (absorption + scattering) * length;
        if (tau >= optical_thick) return RadiationApproximation::Rosseland;
        if (tau >= optical_thin) return RadiationApproximation::P1;
        return RadiationApproximation::DOM;
    }
};

// Rosseland diffusion is an optically-thick approximation. For a gray,
// isotropically scattering medium the transport opacity is represented here
// by the extinction coefficient absorption+scattering.
inline double rosseland_conductivity(double temperature, double absorption)
{
    if (!std::isfinite(temperature) || !std::isfinite(absorption) ||
        temperature <= 0.0 || absorption <= 0.0)
        throw std::invalid_argument("invalid Rosseland state");
    constexpr double sigma = 5.670374419e-8;
    return 16.0 * sigma * std::pow(temperature, 3) / (3.0 * absorption);
}

inline double rosseland_conductivity(double temperature,
                                     double absorption,
                                     double scattering)
{
    if (!std::isfinite(scattering) || scattering < 0.0)
        throw std::invalid_argument("invalid Rosseland scattering");
    return rosseland_conductivity(temperature, absorption + scattering);
}

// P1 diffusion coefficient D = 1/[3 (kappa_a + sigma_s)] for isotropic
// scattering. This is the coefficient multiplying grad(G) in the P1 equation.
inline double p1_diffusion_coefficient(double absorption, double scattering)
{
    if (!std::isfinite(absorption) || !std::isfinite(scattering) ||
        absorption < 0.0 || scattering < 0.0 ||
        absorption + scattering <= 0.0)
        throw std::invalid_argument("invalid P1 coefficients");
    return 1.0 / (3.0 * (absorption + scattering));
}

// Kept as a compact reaction coefficient helper for callers that need the
// absorption/extinction product in reduced P1 algebra.
inline double p1_absorption_coefficient(double absorption, double scattering)
{
    if (!std::isfinite(absorption) || !std::isfinite(scattering) ||
        absorption < 0.0 || scattering < 0.0)
        throw std::invalid_argument("invalid P1 coefficients");
    return 3.0 * absorption * (absorption + scattering);
}

} // namespace cfdx::physics
