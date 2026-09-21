#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;

inline double blackbody_emissive_power(double temperature)
{
    if (temperature < 0.0)
        throw std::invalid_argument("blackbody_emissive_power: temperature must be non-negative");
    return STEFAN_BOLTZMANN * std::pow(temperature, 4);
}

// Blackbody radiance/intensity integrated over a differential solid angle.
// E_b = pi I_b for a Lambertian blackbody.
inline double blackbody_intensity(double temperature)
{
    return blackbody_emissive_power(temperature) / M_PI;
}

inline double gray_surface_emissivity_flux(double emissivity,
                                           double temperature,
                                           double irradiation)
{
    if (emissivity < 0.0 || emissivity > 1.0 || temperature < 0.0 ||
        irradiation < 0.0)
        throw std::invalid_argument("gray_surface_emissivity_flux: invalid input");
    const double Eb = blackbody_emissive_power(temperature);
    return emissivity * (Eb - irradiation);
}

inline double two_surface_net_exchange(double emissivity1, double emissivity2,
                                       double T1, double T2, double view_factor_12)
{
    if (emissivity1 <= 0.0 || emissivity1 > 1.0 ||
        emissivity2 <= 0.0 || emissivity2 > 1.0 ||
        T1 < 0.0 || T2 < 0.0 || view_factor_12 < 0.0 || view_factor_12 > 1.0)
        throw std::invalid_argument("two_surface_net_exchange: invalid input");
    if (view_factor_12 == 0.0)
        return 0.0;
    const double resistance =
        (1.0 - emissivity1) / emissivity1 +
        1.0 / view_factor_12 +
        (1.0 - emissivity2) / emissivity2;
    if (!(resistance > 0.0))
        throw std::invalid_argument("two_surface_net_exchange: invalid view factor");
    return STEFAN_BOLTZMANN * (std::pow(T1, 4) - std::pow(T2, 4)) / resistance;
}

inline void validate_view_factor_matrix(const std::vector<double>& F,
                                        std::size_t n, double tol = 1e-10)
{
    if (F.size() != n * n)
        throw std::invalid_argument("view factors: matrix size mismatch");
    for (std::size_t i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double f = F[i*n+j];
            if (f < -tol || f > 1.0 + tol)
                throw std::invalid_argument("view factors: entry outside [0,1]");
            row_sum += f;
        }
        if (std::abs(row_sum - 1.0) > tol)
            throw std::invalid_argument("view factors: enclosure row does not sum to one");
    }
}

inline double p1_radiative_source(double absorption, double mean_intensity,
                                  double temperature)
{
    if (absorption < 0.0 || mean_intensity < 0.0 || temperature < 0.0)
        throw std::invalid_argument("p1_radiative_source: invalid input");
    // S_rad = 4 kappa (I_b - J), with I_b = sigma T^4 / pi.
    return 4.0 * absorption *
           (blackbody_intensity(temperature) - mean_intensity);
}

struct DiscreteDirection {
    double dx, dy, dz, weight;
};

inline void validate_discrete_directions(const std::vector<DiscreteDirection>& directions,
                                         double tol = 1e-10)
{
    if (directions.empty())
        throw std::invalid_argument("DOM: at least one direction is required");
    double weight_sum = 0.0;
    for (const auto& d : directions) {
        const double norm = std::sqrt(d.dx*d.dx + d.dy*d.dy + d.dz*d.dz);
        if (std::abs(norm - 1.0) > tol || d.weight < 0.0)
            throw std::invalid_argument("DOM: direction must be normalized and weight non-negative");
        weight_sum += d.weight;
    }
    if (std::abs(weight_sum - 4.0 * M_PI) > 1e-8)
        throw std::invalid_argument("DOM: quadrature weights must sum to 4*pi");
}

}  // namespace cfdx::physics
