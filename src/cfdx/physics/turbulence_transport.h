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

enum class TurbulenceModel { LAMINAR, KEPSILON, KOMEGA, SST, SPALART_ALLMARAS, SMAGORINSKY, DES };

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
    // Spalart-Allmaras controls (SI kinematic-viscosity formulation).
    double sa_cb1 = 0.1355;
    double sa_cb2 = 0.622;
    double sa_sigma = 2.0/3.0;
    double sa_kappa = 0.41;
    double sa_cw2 = 0.3;
    double sa_cw3 = 2.0;
    double sa_cv1 = 7.1;
    double sa_ct3 = 1.2;
    double sa_ct4 = 0.5;
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
            (c.model==TurbulenceModel::SST || c.model==TurbulenceModel::KOMEGA)
                ? c.omega_min : (c.model==TurbulenceModel::SPALART_ALLMARAS ? 0.0 : c.epsilon_min));
    }
}

inline double turbulence_nu_t(
    double k, double second, double strain, double wall_distance,
    const TurbulenceTransportControls& c)
{
    k=std::max(k,c.k_min);
    second=std::max(second,
        (c.model==TurbulenceModel::SST || c.model==TurbulenceModel::KOMEGA)
            ? c.omega_min : (c.model==TurbulenceModel::SPALART_ALLMARAS ? 0.0 : c.epsilon_min));
    switch(c.model) {
    case TurbulenceModel::LAMINAR: return 0.0;
    case TurbulenceModel::KEPSILON:
        return c.C_mu*k*k/second;
    case TurbulenceModel::KOMEGA:
        return c.a1*k/second;
    case TurbulenceModel::SST: {
        const double F2=1.0;
        return c.a1*k/std::max(c.a1*second,strain*F2);
    }
    case TurbulenceModel::SPALART_ALLMARAS:
        return second;
    case TurbulenceModel::SMAGORINSKY: {
        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
        const double delta=std::cbrt(1.0);
        return smagorinsky_eddy_viscosity(delta,strain);
    }
    case TurbulenceModel::DES: {
        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
        const double length=std::min(wall_distance,0.65*std::cbrt(1.0));
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
    const double values[] = {c.C_mu,c.C1,c.C2,c.beta_star,c.beta1,c.beta2,
        c.gamma1,c.gamma2,c.a1,c.sa_cb1,c.sa_cb2,c.sa_sigma,c.sa_kappa,
        c.sa_cw2,c.sa_cw3,c.sa_cv1,c.sa_ct3,c.sa_ct4};
    for(double v : values) if(!std::isfinite(v))
        throw std::invalid_argument("non-finite turbulence coefficient");
    if(c.C_mu<=0.0 || c.C1<0.0 || c.C2<0.0 || c.beta_star<=0.0 ||
       c.beta1<=0.0 || c.beta2<=0.0 || c.gamma1<=0.0 || c.gamma2<=0.0 ||
       c.a1<=0.0 || c.sigma_k<=0.0 || c.sigma_epsilon<=0.0 ||
       c.sa_cb1<=0.0 || c.sa_cb2<0.0 || c.sa_sigma<=0.0 || c.sa_kappa<=0.0 ||
       c.sa_cw2<0.0 || c.sa_cw3<=0.0 || c.sa_cv1<=0.0 || c.sa_ct3<0.0 || c.sa_ct4<0.0)
        throw std::invalid_argument("invalid turbulence coefficients");
}

} // namespace cfdx::physics
