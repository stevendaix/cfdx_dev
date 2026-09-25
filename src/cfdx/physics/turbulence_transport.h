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

enum class TurbulenceModel { LAMINAR, KEPSILON, RNG_KEPSILON, KOMEGA, SST, SPALART_ALLMARAS, SMAGORINSKY, DES };

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
    // RNG k-epsilon constants (Yakhot et al.; OpenFOAM v13 defaults).
    double rng_C_mu = 0.0845;
    double rng_C1 = 1.42;
    double rng_C2 = 1.68;
    double rng_sigma_k = 0.71942;
    double rng_sigma_epsilon = 0.71942;
    double rng_eta0 = 4.38;
    double rng_beta = 0.012;
    // Spalart-Allmaras constants, using kinematic nu-tilde.
    double sa_cb1 = 0.1355, sa_cb2 = 0.622, sa_sigma = 2.0/3.0;
    double sa_kappa = 0.41, sa_cw2 = 0.3, sa_cw3 = 2.0, sa_cv1 = 7.1;
    double sa_ct3 = 1.2, sa_ct4 = 0.5;
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
            (c.model==TurbulenceModel::SST || c.model==TurbulenceModel::KOMEGA) ? c.omega_min : (c.model==TurbulenceModel::SPALART_ALLMARAS ? 0.0 : c.epsilon_min));
    }
}

inline double turbulence_nu_t(
    double k, double second, double strain, double wall_distance,
    const TurbulenceTransportControls& c, double cell_volume = -1.0, double F2 = 1.0)
{
    k=std::max(k,c.k_min);
    const bool omega_based =
        c.model==TurbulenceModel::SST || c.model==TurbulenceModel::KOMEGA;
    if(c.model!=TurbulenceModel::SPALART_ALLMARAS)
        second=std::max(second, omega_based ? c.omega_min : c.epsilon_min);
    switch(c.model) {
    case TurbulenceModel::LAMINAR: return 0.0;
    case TurbulenceModel::KEPSILON: return c.C_mu*k*k/second;
    case TurbulenceModel::RNG_KEPSILON: return c.rng_C_mu*k*k/second;
    case TurbulenceModel::KOMEGA: {
        const double omega_tilde=std::max(
            second, c.komega_clim*std::max(strain,0.0)/std::sqrt(c.beta_star));
        return k/omega_tilde;
    }
    case TurbulenceModel::SPALART_ALLMARAS: {
        const double nt=std::max(second,0.0);
        if(!(c.molecular_viscosity>0.0))
            throw std::invalid_argument("SA eddy viscosity requires molecular viscosity");
        const double chi3=std::pow(nt/c.molecular_viscosity,3.0);
        return nt*chi3/(chi3+std::pow(c.sa_cv1,3.0));
    }
    case TurbulenceModel::SST:
        if(F2<0.0 || F2>1.0) throw std::invalid_argument("SST F2 must lie in [0,1]");
        return c.a1*k/std::max(c.a1*second,std::max(strain,0.0)*F2);
    case TurbulenceModel::SMAGORINSKY:
        if(!(cell_volume>0.0) || !std::isfinite(cell_volume))
            throw std::invalid_argument("Smagorinsky requires positive cell volume");
        return smagorinsky_eddy_viscosity(std::cbrt(cell_volume),strain,c.smagorinsky_Cs);
    case TurbulenceModel::DES:
        if(!(cell_volume>0.0) || !std::isfinite(cell_volume) || !(wall_distance>0.0))
            throw std::invalid_argument("DES requires positive cell volume and wall distance");
        return des_eddy_viscosity(std::cbrt(cell_volume),wall_distance,strain,
                                  c.smagorinsky_Cs,c.des_Cdes);
    }
    throw std::invalid_argument("unknown turbulence model");
}

inline void validate_turbulence_controls(const TurbulenceTransportControls& c)
{
    if(!std::isfinite(c.density) || !std::isfinite(c.molecular_viscosity) ||
       !std::isfinite(c.turbulent_prandtl) || !std::isfinite(c.k_min) ||
       !std::isfinite(c.epsilon_min) || !std::isfinite(c.omega_min) ||
       c.density<=0.0 || c.molecular_viscosity<0.0 ||
       c.turbulent_prandtl<=0.0 || c.k_min<=0.0 ||
       c.epsilon_min<=0.0 || c.omega_min<=0.0)
        throw std::invalid_argument("invalid turbulence controls");
    const double values[] = {c.C_mu,c.C1,c.C2,c.beta_star,c.beta1,c.beta2,
        c.gamma1,c.gamma2,c.a1,c.rng_C_mu,c.rng_C1,c.rng_C2,c.rng_sigma_k,
        c.rng_sigma_epsilon,c.rng_eta0,c.rng_beta,c.sa_cb1,c.sa_cb2,c.sa_sigma,
        c.sa_kappa,c.sa_cw2,c.sa_cw3,c.sa_cv1,c.sa_ct3,c.sa_ct4};
    for (double v : values)
        if (!std::isfinite(v)) throw std::invalid_argument("non-finite turbulence coefficient");
    if(c.C_mu<=0.0 || c.C1<0.0 || c.C2<0.0 || c.beta_star<=0.0 ||
       c.beta1<=0.0 || c.beta2<=0.0 || c.gamma1<=0.0 || c.gamma2<=0.0 || c.a1<=0.0 ||
       c.sigma_k<=0.0 || c.sigma_epsilon<=0.0 || c.rng_C_mu<=0.0 || c.rng_C1<=0.0 ||
       c.rng_C2<=0.0 || c.rng_sigma_k<=0.0 || c.rng_sigma_epsilon<=0.0 ||
       c.rng_eta0<=0.0 || c.rng_beta<=0.0 || c.sa_cb1<=0.0 || c.sa_cb2<0.0 ||
       c.sa_sigma<=0.0 || c.sa_kappa<=0.0 || c.sa_cw2<0.0 || c.sa_cw3<=0.0 ||
       c.sa_cv1<=0.0 || c.sa_ct3<0.0 || c.sa_ct4<0.0)
        throw std::invalid_argument("invalid turbulence coefficients");
}

} // namespace cfdx::physics
