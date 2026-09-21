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

enum class TurbulenceModel { LAMINAR, KEPSILON, SST, SMAGORINSKY, DES };

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
    case TurbulenceModel::SST: {
        const double F2=1.0;
        return c.a1*k/std::max(c.a1*second,strain*F2);
    }
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
    if(c.C_mu<=0.0 || c.C1<0.0 || c.C2<0.0 ||
       c.sigma_k<=0.0 || c.sigma_epsilon<=0.0)
        throw std::invalid_argument("invalid k-epsilon coefficients");
}

} // namespace cfdx::physics
