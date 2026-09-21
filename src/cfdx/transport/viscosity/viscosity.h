// Transport Models - Viscosity
// 
// Architecture réorganisée P2:
//   viscosity.h - Dynamic viscosity models

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

namespace cfdx {
namespace transport {
namespace viscosity {

struct SutherlandParams {
    double mu0 = 1.716e-5;    // Reference dynamic viscosity [Pa.s] at T0
    double T0 = 273.15;       // Reference temperature [K]
    double S = 110.4;         // Sutherland constant [K]
};

inline double sutherland(double T, const SutherlandParams& params) {
    if (!std::isfinite(T) || T <= 0.0) {
        throw std::invalid_argument("sutherland: temperature must be finite and strictly positive");
    }
    if (!std::isfinite(params.mu0) || params.mu0 <= 0.0 ||
        !std::isfinite(params.T0) || params.T0 <= 0.0 ||
        !std::isfinite(params.S) || params.S <= -params.T0) {
        throw std::invalid_argument("sutherland: invalid model parameters");
    }
    const double mu = params.mu0 * std::pow(T / params.T0, 1.5) *
                      (params.T0 + params.S) / (T + params.S);
    if (!std::isfinite(mu) || mu <= 0.0) {
        throw std::runtime_error("sutherland: non-finite or non-positive viscosity");
    }
    return mu;
}

inline double sutherland(double T, double mu0_ref = 1.716e-5, double T_ref = 273.15, double S = 110.4) {
    return sutherland(T, SutherlandParams{mu0_ref, T_ref, S});
}

inline double power_law(double T, double mu0, double T0, double n) {
    if (!std::isfinite(T) || T <= 0.0 || !std::isfinite(mu0) || mu0 <= 0.0 ||
        !std::isfinite(T0) || T0 <= 0.0 || !std::isfinite(n)) {
        throw std::invalid_argument("power_law: invalid temperature or model parameters");
    }
    const double mu = mu0 * std::pow(T / T0, n);
    if (!std::isfinite(mu) || mu <= 0.0) {
        throw std::runtime_error("power_law: non-finite or non-positive viscosity");
    }
    return mu;
}

inline double constant(double mu) {
    return mu;
}

}  // namespace viscosity
}  // namespace transport
}  // namespace cfdx