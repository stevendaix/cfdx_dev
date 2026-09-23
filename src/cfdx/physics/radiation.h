#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;

inline void validate_finite(double value, const char* name)
{
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be finite");
}

inline double blackbody_emissive_power(double temperature)
{
    validate_finite(temperature, "temperature");
    if (temperature < 0.0)
        throw std::invalid_argument("blackbody_emissive_power: temperature must be non-negative");
    return STEFAN_BOLTZMANN * std::pow(temperature, 4);
}

inline double blackbody_intensity(double temperature)
{
    return blackbody_emissive_power(temperature) / M_PI;
}

inline double gray_surface_emissivity_flux(double emissivity,
                                           double temperature,
                                           double irradiation)
{
    validate_finite(emissivity, "emissivity");
    validate_finite(temperature, "temperature");
    validate_finite(irradiation, "irradiation");
    if (emissivity < 0.0 || emissivity > 1.0 || temperature < 0.0 ||
        irradiation < 0.0)
        throw std::invalid_argument("gray_surface_emissivity_flux: invalid input");
    const double Eb = blackbody_emissive_power(temperature);
    return emissivity * (Eb - irradiation);
}

// Diffuse-gray wall intensity emitted/reflected into the domain.
// I_w = epsilon*I_b + (1-epsilon)*G/(4*pi).
inline double gray_diffuse_wall_intensity(double emissivity,
                                          double temperature,
                                          double irradiation)
{
    validate_finite(emissivity, "emissivity");
    validate_finite(temperature, "temperature");
    validate_finite(irradiation, "irradiation");
    if (emissivity < 0.0 || emissivity > 1.0 || temperature < 0.0 ||
        irradiation < 0.0)
        throw std::invalid_argument("gray_diffuse_wall_intensity: invalid input");
    return emissivity * blackbody_intensity(temperature) +
           (1.0 - emissivity) * irradiation / (4.0 * M_PI);
}

// Net exchange per unit area of surface 1. The legacy overload assumes
// A1=A2; the area-aware overload must be used for general enclosures.
inline double two_surface_net_exchange(double emissivity1, double emissivity2,
                                       double T1, double T2, double view_factor_12,
                                       double area1, double area2)
{
    for (const auto& v : {emissivity1, emissivity2, T1, T2,
                          view_factor_12, area1, area2})
        validate_finite(v, "two_surface_net_exchange input");
    if (emissivity1 <= 0.0 || emissivity1 > 1.0 ||
        emissivity2 <= 0.0 || emissivity2 > 1.0 ||
        T1 < 0.0 || T2 < 0.0 || view_factor_12 < 0.0 || view_factor_12 > 1.0 ||
        !(area1 > 0.0) || !(area2 > 0.0))
        throw std::invalid_argument("two_surface_net_exchange: invalid input");
    if (view_factor_12 == 0.0)
        return 0.0;
    const double resistance =
        (1.0 - emissivity1) / emissivity1 +
        1.0 / view_factor_12 +
        (1.0 - emissivity2) / emissivity2 * (area1 / area2);
    if (!(resistance > 0.0))
        throw std::invalid_argument("two_surface_net_exchange: invalid view factor");
    return STEFAN_BOLTZMANN * (std::pow(T1, 4) - std::pow(T2, 4)) / resistance;
}

inline double two_surface_net_exchange(double emissivity1, double emissivity2,
                                       double T1, double T2, double view_factor_12)
{
    return two_surface_net_exchange(emissivity1, emissivity2, T1, T2,
                                    view_factor_12, 1.0, 1.0);
}

inline void validate_view_factor_matrix(const std::vector<double>& F,
                                        std::size_t n, double tol = 1e-10)
{
    if (F.size() != n * n || n == 0)
        throw std::invalid_argument("view factors: matrix size mismatch");
    if (!(tol > 0.0) || !std::isfinite(tol))
        throw std::invalid_argument("view factors: invalid tolerance");
    for (std::size_t i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double f = F[i*n+j];
            if (!std::isfinite(f) || f < -tol || f > 1.0 + tol)
                throw std::invalid_argument("view factors: entry outside [0,1]");
            row_sum += f;
        }
        if (std::abs(row_sum - 1.0) > tol)
            throw std::invalid_argument("view factors: enclosure row does not sum to one");
    }
}

inline void validate_view_factor_matrix(const std::vector<double>& F,
                                        std::size_t n,
                                        const std::vector<double>& areas,
                                        double tol = 1e-10)
{
    validate_view_factor_matrix(F, n, tol);
    if (areas.size() != n)
        throw std::invalid_argument("view factors: area size mismatch");
    for (double a : areas) {
        if (!std::isfinite(a) || !(a > 0.0))
            throw std::invalid_argument("view factors: areas must be finite and positive");
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            if (std::abs(areas[i] * F[i*n+j] - areas[j] * F[j*n+i]) > tol *
                std::max({1.0, std::abs(areas[i] * F[i*n+j]),
                          std::abs(areas[j] * F[j*n+i])}))
                throw std::invalid_argument("view factors: reciprocity violated");
}

inline double p1_radiative_source(double absorption, double mean_intensity,
                                  double temperature)
{
    validate_finite(absorption, "absorption");
    validate_finite(mean_intensity, "mean_intensity");
    validate_finite(temperature, "temperature");
    if (absorption < 0.0 || mean_intensity < 0.0 || temperature < 0.0)
        throw std::invalid_argument("p1_radiative_source: invalid input");
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
    if (!(tol > 0.0) || !std::isfinite(tol))
        throw std::invalid_argument("DOM: invalid tolerance");
    double weight_sum = 0.0;
    double mx = 0.0, my = 0.0, mz = 0.0;
    double mxx = 0.0, myy = 0.0, mzz = 0.0;
    double mxy = 0.0, mxz = 0.0, myz = 0.0;
    for (const auto& d : directions) {
        const double norm = std::sqrt(d.dx*d.dx + d.dy*d.dy + d.dz*d.dz);
        if (!std::isfinite(norm) || std::abs(norm - 1.0) > tol ||
            !std::isfinite(d.weight) || d.weight < 0.0)
            throw std::invalid_argument("DOM: direction must be normalized and weight non-negative");
        weight_sum += d.weight;
        mx += d.weight*d.dx; my += d.weight*d.dy; mz += d.weight*d.dz;
        mxx += d.weight*d.dx*d.dx; myy += d.weight*d.dy*d.dy;
        mzz += d.weight*d.dz*d.dz;
        mxy += d.weight*d.dx*d.dy; mxz += d.weight*d.dx*d.dz;
        myz += d.weight*d.dy*d.dz;
    }
    const double scale = 4.0 * M_PI;
    if (std::abs(weight_sum - scale) > 1e-8 ||
        std::abs(mx) > tol*scale || std::abs(my) > tol*scale ||
        std::abs(mz) > tol*scale ||
        std::abs(mxx - scale/3.0) > tol*scale ||
        std::abs(myy - scale/3.0) > tol*scale ||
        std::abs(mzz - scale/3.0) > tol*scale ||
        std::abs(mxy) > tol*scale || std::abs(mxz) > tol*scale ||
        std::abs(myz) > tol*scale)
        throw std::invalid_argument("DOM: quadrature does not reproduce isotropic first/second moments");
}

}  // namespace cfdx::physics
