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
        if (absorption < 0.0 || scattering < 0.0 || length <= 0.0 ||
            optical_thick <= optical_thin || optical_thin < 0.0) {
            throw std::invalid_argument("invalid radiation regime");
        }
        const double tau = (absorption + scattering) * length;
        if (tau >= optical_thick) return RadiationApproximation::Rosseland;
        if (tau >= optical_thin) return RadiationApproximation::P1;
        return RadiationApproximation::DOM;
    }
};

inline double rosseland_conductivity(double temperature, double absorption) {
    if (temperature <= 0.0 || absorption <= 0.0)
        throw std::invalid_argument("invalid Rosseland state");
    constexpr double sigma = 5.670374419e-8;
    return 16.0 * sigma * std::pow(temperature, 3) / (3.0 * absorption);
}

inline double p1_absorption_coefficient(double absorption, double scattering) {
    if (absorption < 0.0 || scattering < 0.0)
        throw std::invalid_argument("invalid P1 coefficients");
    return 3.0 * absorption * (absorption + scattering);
}

} // namespace cfdx::physics
