// Transport Models - Viscosity
// 
// Architecture réorganisée P2:
//   viscosity.h - Dynamic viscosity models

#pragma once

#include <cmath>

namespace cfdx {
namespace transport {
namespace viscosity {

struct SutherlandParams {
    double mu0 = 1.716e-5;    // Reference dynamic viscosity [Pa.s] at T0
    double T0 = 273.15;       // Reference temperature [K]
    double S = 110.4;         // Sutherland constant [K]
};

inline double sutherland(double T, const SutherlandParams& params) {
    return params.mu0 * std::pow(T / params.T0, 1.5) * (params.T0 + params.S) / (T + params.S);
}

inline double sutherland(double T, double mu0_ref = 1.716e-5, double T_ref = 273.15, double S = 110.4) {
    return mu0_ref * std::pow(T / T_ref, 1.5) * (T_ref + S) / (T + S);
}

inline double power_law(double T, double mu0, double T0, double n) {
    return mu0 * std::pow(T / T0, n);
}

inline double constant(double mu) {
    return mu;
}

}  // namespace viscosity
}  // namespace transport
}  // namespace cfdx