// M0.14-T02 — Transport Models
//
// Viscosity, thermal conductivity, diffusivity models

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <cmath>

namespace cfdx {
namespace physics {

// Constantes pour l'air (Sutherland)
struct SutherlandParams {
    double mu0 = 1.716e-5;   // Viscosité de référence [Pa.s] à T0
    double T0 = 273.15;      // Température de référence [K]
    double S = 110.4;        // Constante de Sutherland [K]
};

// Viscosité constante
inline double constant_viscosity(double mu) {
    return mu;
}

// Loi de Sutherland pour viscosité
inline double sutherland_viscosity(double T, const SutherlandParams& params) {
    return params.mu0 * std::pow(T / params.T0, 1.5) * (params.T0 + params.S) / (T + params.S);
}

// Loi de puissance pour viscosité
inline double power_law_viscosity(double T, double mu0, double T0, double n) {
    return mu0 * std::pow(T / T0, n);
}

// Conductivité thermique constante
inline double constant_conductivity(double k) {
    return k;
}

// Conductivité thermique via nombre de Prandtl : k = mu * Cp / Pr
inline double prandtl_conductivity(double mu, double Cp, double Pr) {
    return mu * Cp / Pr;
}

// Diffusivité constante
inline double constant_diffusivity(double D) {
    return D;
}

// Diffusivité via nombre de Schmidt : D = mu / (rho * Sc)
inline double schmidt_diffusivity(double mu, double rho, double Sc) {
    return mu / (rho * Sc);
}

// Propriétés de transport regroupées
struct TransportProperties {
    double mu = 1.8e-5;      // Viscosité dynamique [Pa.s]
    double k = 0.026;        // Conductivité thermique [W/m/K]
    double D = 2.0e-5;       // Diffusivité massique [m²/s]
    double Cp = 1004.5;      // Chaleur spécifique [J/kg/K]
    double Pr = 0.71;        // Nombre de Prandtl
    double Sc = 0.7;         // Nombre de Schmidt
};

// Compute toutes les propriétés de transport à partir de T
inline TransportProperties compute_transport(double T,
                                              const SutherlandParams& suth = {},
                                              double Pr = 0.71,
                                              double Sc = 0.7,
                                              double Cp = 1004.5)
{
    TransportProperties tp;
    tp.mu = sutherland_viscosity(T, suth);
    tp.k = prandtl_conductivity(tp.mu, Cp, Pr);
    tp.D = schmidt_diffusivity(tp.mu, 1.2, Sc);  // rho ≈ 1.2 kg/m³ pour l'air
    tp.Cp = Cp;
    tp.Pr = Pr;
    tp.Sc = Sc;
    return tp;
}

// Compute sur un champ de température
inline void compute_transport_fields(const Field<double, Location::CELL>& T,
                                      Field<double, Location::CELL>& mu,
                                      Field<double, Location::CELL>& k,
                                      Field<double, Location::CELL>& D,
                                      const SutherlandParams& suth = {},
                                      double Pr = 0.71,
                                      double Sc = 0.7,
                                      double Cp = 1004.5)
{
    const std::size_t n = T.size();
    for (std::size_t c = 0; c < n; ++c) {
        TransportProperties tp = compute_transport(T(c), suth, Pr, Sc, Cp);
        mu(c) = tp.mu;
        k(c) = tp.k;
        D(c) = tp.D;
    }
}

}  // namespace physics
}  // namespace cfdx