#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/turbulence.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>

namespace cfdx::physics {

enum class TurbulenceModel {
    LAMINAR,
    KEPSILON,
    RNG_KEPSILON,
    REALIZABLE_KEPSILON,
    KOMEGA,
    SST,
    SPALART_ALLMARAS,
    SMAGORINSKY,
    WALE,
    DES,
    DDES,
    IDDES
};

struct TurbulenceTransportControls {
    TurbulenceModel model = TurbulenceModel::LAMINAR;
    double density = 1.0;
    double molecular_viscosity = 1e-3;
    double turbulent_prandtl = 0.9;
    double k_min = 1e-12;
    double epsilon_min = 1e-12;
    double omega_min = 1e-12;
    double sigma_k = 1.0;
    double sigma_epsilon = 1.3;
    double C_mu = 0.09;
    double C1 = 1.44;
    double C2 = 1.92;
    double beta_star = 0.09;
    double beta1 = 0.075;
    double beta2 = 0.0828;
    double gamma1 = 5.0/9.0;
    double gamma2 = 0.44;
    double a1 = 0.31;
    double rng_C_mu = 0.0845;
    double rng_C1 = 1.42;
    double rng_C2 = 1.68;
    double rng_eta0 = 4.38;
    double rng_beta = 0.012;
    double realizable_A0 = 4.04;
    double realizable_As = 0.0;
    double realizable_C1 = 1.44;
    double realizable_C2 = 1.9;
    double sigma_omega1 = 0.5;
    double sigma_omega2 = 0.856;
    double sigma_d = 0.0;
    double C_w1 = 0.492;
    double C_w2 = 0.0;
    double C_w3 = 2.0;
    double sigma_sa = 2.0/3.0;
    double Cb1 = 0.1355;
    double Cb2 = 0.622;
    double kappa = 0.41;
    double C_DES = 0.65;
    double C_DDES = 0.65;
    double C_IDDES = 0.65;
    double smagorinsky_constant = 0.17;
    double wale_constant = 0.5;
    double filter_width = 1.0;
};

inline void enforce_turbulence_bounds(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& second,
    const TurbulenceTransportControls& c)
{
    if (k.size()!=second.size()) throw std::invalid_argument("turbulence field size mismatch");
    for (std::size_t i=0;i<k.size();++i) {
        k(i)=std::max(k(i),c.k_min);
        second(i)=std::max(second(i),
            c.model==TurbulenceModel::SST ? c.omega_min : c.epsilon_min);
    }
}

inline double turbulence_nu_t(
    double k, double second, double strain, double wall_distance,
    const TurbulenceTransportControls& c)
{
    k=std::max(k,c.k_min);
    second=std::max(second,
        c.model==TurbulenceModel::SST ? c.omega_min : c.epsilon_min);
    switch(c.model) {
    case TurbulenceModel::LAMINAR: return 0.0;
    case TurbulenceModel::KEPSILON:
        return c.C_mu*k*k/second;
    case TurbulenceModel::RNG_KEPSILON: {
        const double eta = std::max(0.0, strain*k/second);
        const double fmu = std::exp(-0.0165*eta*eta);
        return c.rng_C_mu*fmu*k*k/second;
    }
    case TurbulenceModel::REALIZABLE_KEPSILON: {
        const double eta = std::max(0.0, strain*k/second);
        const double A = c.realizable_A0 + std::sqrt(6.0)*std::max(0.0,std::abs(strain));
        const double cmu = 1.0/std::max(A,1e-12);
        return cmu*k*k/second;
    }
    case TurbulenceModel::KOMEGA: {
        const double F2 = 1.0;
        return c.a1*k/std::max(c.a1*second,strain*F2);
    }
    case TurbulenceModel::SST: {
        const double F2=1.0;
        return c.a1*k/std::max(c.a1*second,strain*F2);
    }
    case TurbulenceModel::SPALART_ALLMARAS:
        return std::max(0.0, second);
    case TurbulenceModel::SMAGORINSKY: {
        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
        const double delta=std::max(c.filter_width,1e-12);
        return smagorinsky_eddy_viscosity(c.smagorinsky_constant*delta,strain);
    }
    case TurbulenceModel::WALE: {
        const double delta=std::max(c.filter_width,1e-12);
        const double s2=strain*strain;
        return std::pow(c.wale_constant*delta,2.0)*std::sqrt(std::max(s2,0.0));
    }
    case TurbulenceModel::DES:
    case TurbulenceModel::DDES:
    case TurbulenceModel::IDDES: {
        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
        const double shielding = (c.model==TurbulenceModel::DES) ? 1.0 : std::min(1.0,std::max(0.0,wall_distance/std::max(c.filter_width,1e-12)));
        const double cd = c.model==TurbulenceModel::IDDES ? c.C_IDDES : (c.model==TurbulenceModel::DDES ? c.C_DDES : c.C_DES);
        const double length=std::min(wall_distance,cd*std::cbrt(std::max(c.filter_width,1e-12)))*shielding;
        return length*length*strain;
    }
    }
    throw std::invalid_argument("unknown turbulence model");
}

inline void validate_turbulence_controls(const TurbulenceTransportControls& c)
{
    if(c.density<=0.0 || c.molecular_viscosity<0.0 ||
       c.turbulent_prandtl<=0.0 || c.k_min<=0.0 ||
       c.epsilon_min<=0.0 || c.omega_min<=0.0)
        throw std::invalid_argument("invalid turbulence controls");
    if(c.C_mu<=0.0 || c.C1<0.0 || c.C2<0.0 ||
       c.sigma_k<=0.0 || c.sigma_epsilon<=0.0)
        throw std::invalid_argument("invalid k-epsilon coefficients");
}

} // namespace cfdx::physics
