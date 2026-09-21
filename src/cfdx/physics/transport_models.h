// M0.14-T02 — Transport Models
//
// Viscosité, conductivité thermique, diffusivité
// Basé sur les modèles Sutherland pour l'air

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/field/storage.h"
#include "equation_of_state.h"
#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx {
namespace physics {

using cfdx::core::Location;
using cfdx::core::Field;

constexpr double SUTHERTY_MU0 = 1.716e-5;
constexpr double SUTHERTY_T0 = 273.15;
constexpr double SUTHERTY_S = 110.4;
constexpr double AIR_R = 287.04749097718457;

struct TransportProperties {
    double mu = 1.8e-5;
    double k = 0.026;
    double D = 2.0e-5;
    double rho = 1.2;
    double Cp = 1004.5;
    double Pr = 0.71;
    double Sc = 0.7;
};

inline double sutherland_viscosity(double T, double mu0_ref = SUTHERTY_MU0, double T_ref = SUTHERTY_T0, double S = SUTHERTY_S) {
    if (!std::isfinite(T) || T <= 0.0) throw std::invalid_argument("sutherland_viscosity: temperature must be finite and strictly positive");
    if (!std::isfinite(mu0_ref) || mu0_ref <= 0.0 || !std::isfinite(T_ref) || T_ref <= 0.0 ||
        !std::isfinite(S) || S <= -T_ref) throw std::invalid_argument("sutherland_viscosity: invalid parameters");
    const double mu = mu0_ref * std::pow(T / T_ref, 1.5) * (T_ref + S) / (T + S);
    if (!std::isfinite(mu) || mu <= 0.0) throw std::runtime_error("sutherland_viscosity: non-finite or non-positive result");
    return mu;
}

inline double power_law_viscosity(double T, double mu0, double T0, double n) {
    if (!std::isfinite(T) || T <= 0.0 || !std::isfinite(mu0) || mu0 <= 0.0 ||
        !std::isfinite(T0) || T0 <= 0.0 || !std::isfinite(n))
        throw std::invalid_argument("power_law_viscosity: invalid parameters");
    const double mu = mu0 * std::pow(T / T0, n);
    if (!std::isfinite(mu) || mu <= 0.0) throw std::runtime_error("power_law_viscosity: invalid result");
    return mu;
}

inline double constant_viscosity(double mu) {
    return mu;
}

inline double prandtl_conductivity(double mu, double Cp, double Pr) {
    if (!std::isfinite(mu) || mu < 0.0 || !std::isfinite(Cp) || Cp <= 0.0 || !std::isfinite(Pr) || Pr <= 0.0)
        throw std::invalid_argument("prandtl_conductivity: invalid transport parameters");
    return mu * Cp / Pr;
}

inline double schmidt_diffusivity(double mu, double rho, double Sc) {
    if (!std::isfinite(mu) || mu < 0.0 || !std::isfinite(rho) || rho <= 0.0 ||
        !std::isfinite(Sc) || Sc <= 0.0) throw std::invalid_argument("schmidt_diffusivity: invalid transport parameters");
    return mu / (rho * Sc);
}

inline TransportProperties compute_transport(double T, double rho,
                                             double mu0_ref = SUTHERTY_MU0,
                                             double T_ref = SUTHERTY_T0,
                                             double S = SUTHERTY_S,
                                             double Pr = 0.71,
                                             double Sc = 0.7,
                                             double Cp = AIR_R * 1005.0 / 1.4) {
    TransportProperties tp;
    tp.mu = sutherland_viscosity(T, mu0_ref, T_ref, S);
    tp.k = prandtl_conductivity(tp.mu, Cp, Pr);
    tp.D = schmidt_diffusivity(tp.mu, rho, Sc);
    tp.rho = rho;
    tp.Cp = Cp;
    tp.Pr = Pr;
    tp.Sc = Sc;
    return tp;
}

inline TransportProperties compute_transport(double T,
                                              const IdealGasEOS& eos,
                                              double p,
                                              double Pr = 0.71,
                                              double Sc = 0.7) {
    double rho = eos.density(p, T);
    double mu = sutherland_viscosity(T);
    double Cp = eos.cp(p, T);
    TransportProperties tp;
    tp.mu = mu;
    tp.k = prandtl_conductivity(mu, Cp, Pr);
    tp.D = schmidt_diffusivity(mu, rho, Sc);
    tp.rho = rho;
    tp.Cp = Cp;
    tp.Pr = Pr;
    tp.Sc = Sc;
    return tp;
}

inline void compute_transport_fields(const Field<double, Location::CELL>& T,
                                      Field<double, Location::CELL>& mu,
                                      Field<double, Location::CELL>& k,
                                      Field<double, Location::CELL>& D,
                                      const IdealGasEOS& eos,
                                      const Field<double, Location::CELL>& p) {
    const std::size_t n = T.size();
    for (std::size_t c = 0; c < n; ++c) {
        TransportProperties tp = compute_transport(T(c), eos, p(c));
        mu(c) = tp.mu;
        k(c) = tp.k;
        D(c) = tp.D;
    }
}

inline void compute_transport_fields(const Field<double, Location::CELL>& T,
                                      Field<double, Location::CELL>& mu,
                                      Field<double, Location::CELL>& k,
                                      Field<double, Location::CELL>& D,
                                      double rho) {
    const std::size_t n = T.size();
    for (std::size_t c = 0; c < n; ++c) {
        TransportProperties tp = compute_transport(T(c), rho);
        mu(c) = tp.mu;
        k(c) = tp.k;
        D(c) = tp.D;
    }
}

}  // namespace physics
}  // namespace cfdx