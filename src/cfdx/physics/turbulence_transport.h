#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/turbulence.h"
#include "cfdx/physics/turbulence_models.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>

namespace cfdx::physics {

using TurbulenceModel = AdvancedTurbulenceModel;

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
    // Realizable k-epsilon defaults used by OpenFOAM's realizableKE model.
    double realizable_A0 = 4.0;
    double realizable_C2 = 1.9;
    double realizable_sigma_k = 1.0;
    double realizable_sigma_epsilon = 1.2;
    // Spalart-Allmaras constants, using kinematic nu-tilde.
    double sa_cb1 = 0.1355, sa_cb2 = 0.622, sa_sigma = 2.0/3.0;
    double sa_kappa = 0.41, sa_cw2 = 0.3, sa_cw3 = 2.0, sa_cv1 = 7.1;
    double sa_ct3 = 1.2, sa_ct4 = 0.5;
    double komega_alpha = 13.0/25.0, komega_beta0 = 0.0708;
    double komega_sigma_k = 0.6, komega_sigma_w = 0.5;
    double komega_sigma_d0 = 1.0/8.0, komega_clim = 7.0/8.0;
    double sst_sigma_k1 = 0.85, sst_sigma_k2 = 1.0;
    double sst_sigma_w1 = 0.5, sst_sigma_w2 = 0.856;
    double sst_production_limiter = 10.0;
    double smagorinsky_Cs = 0.17, wale_Cw = 0.325, des_Cdes = 0.65;
    TurbulenceCorrectionControls corrections{};
};

struct TurbulenceBoundDiagnostics {
    std::size_t first_floor_activations = 0;
    std::size_t second_floor_activations = 0;
    std::size_t first_nonfinite = 0;
    std::size_t second_nonfinite = 0;
    double first_total_correction = 0.0;
    double second_total_correction = 0.0;
    double first_max_correction = 0.0;
    double second_max_correction = 0.0;

    std::size_t nonfinite_values() const noexcept {
        return first_nonfinite + second_nonfinite;
    }
};

// Preserve non-finite values so invalid states remain visible to acceptance.
inline TurbulenceBoundDiagnostics enforce_turbulence_bounds(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& second,
    const TurbulenceTransportControls& c)
{
    if (k.size()!=second.size())
        throw std::invalid_argument("turbulence field size mismatch");
    TurbulenceBoundDiagnostics out;
    const double second_floor =
        (c.model==TurbulenceModel::SST || c.model==TurbulenceModel::KOMEGA)
            ? c.omega_min
            : (c.model==TurbulenceModel::SPALART_ALLMARAS ? 0.0 : c.epsilon_min);
    for (std::size_t i=0;i<k.size();++i) {
        const double ki = k(i);
        if (!std::isfinite(ki)) {
            ++out.first_nonfinite;
        } else if (ki < c.k_min) {
            const double correction = c.k_min - ki;
            k(i) = c.k_min;
            ++out.first_floor_activations;
            out.first_total_correction += correction;
            out.first_max_correction = std::max(out.first_max_correction, correction);
        }
        const double si = second(i);
        if (!std::isfinite(si)) {
            ++out.second_nonfinite;
        } else if (si < second_floor) {
            const double correction = second_floor - si;
            second(i) = second_floor;
            ++out.second_floor_activations;
            out.second_total_correction += correction;
            out.second_max_correction = std::max(out.second_max_correction, correction);
        }
    }
    return out;
}

inline double turbulence_nu_t(
    double k, double second, double strain, double wall_distance,
    const TurbulenceTransportControls& c, double cell_volume = -1.0, double F2 = 1.0,
    double rotation = 0.0, double third_invariant = 0.0, double ddes_r = 0.0,
    double iddes_stress_blend = 1.0, double wale_S2 = -1.0, double wale_Sd2 = -1.0)
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
    case TurbulenceModel::REALIZABLE_KEPSILON: {
        const RealizableKEpsilonInvariants invariants{
            std::max(strain, 0.0), std::max(rotation, 0.0), third_invariant};
        const double cmu = realizable_kepsilon_cmu_from_invariants(
            invariants, k, second, c.realizable_A0);
        return cmu*k*k/second;
    }
    case TurbulenceModel::KOMEGA: {
        return k/second;
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
    case TurbulenceModel::WALE:
        if(!(cell_volume>0.0) || !std::isfinite(cell_volume))
            throw std::invalid_argument("WALE requires positive cell volume");
        if(!(wale_S2>=0.0) || !(wale_Sd2>=0.0))
            throw std::invalid_argument("WALE requires tensor invariants S2 and Sd2");
        return wale_eddy_viscosity({wale_S2,wale_Sd2,std::cbrt(cell_volume)},c.wale_Cw);
    case TurbulenceModel::DES:
    case TurbulenceModel::DDES:
    case TurbulenceModel::IDDES: {
        if(!(cell_volume>0.0) || !std::isfinite(cell_volume) || !(wall_distance>0.0))
            throw std::invalid_argument("DES family requires positive cell volume and wall distance");
        const double delta=std::cbrt(cell_volume);
        if(c.model==TurbulenceModel::DES)
            return des_hybrid_eddy_viscosity(wall_distance,delta,strain,c.smagorinsky_Cs,c.des_Cdes);
        if(c.model==TurbulenceModel::DDES)
            return ddes_hybrid_eddy_viscosity(wall_distance,delta,strain,ddes_r,c.smagorinsky_Cs,c.des_Cdes);
        return iddes_hybrid_eddy_viscosity(wall_distance,delta,strain,ddes_r,iddes_stress_blend,
                                    c.smagorinsky_Cs,c.des_Cdes);
    }
    case TurbulenceModel::DYNAMIC_KEQN:
        if(!(cell_volume>0.0) || !std::isfinite(cell_volume))
            throw std::invalid_argument("dynamic LES requires positive cell volume");
        return dynamic_one_equation_eddy_viscosity(second,std::cbrt(cell_volume),0.1);
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
        c.rng_sigma_epsilon,c.rng_eta0,c.rng_beta,c.realizable_A0,
        c.realizable_C2,c.realizable_sigma_k,c.realizable_sigma_epsilon,
        c.sa_cb1,c.sa_cb2,c.sa_sigma,
        c.sa_kappa,c.sa_cw2,c.sa_cw3,c.sa_cv1,c.sa_ct3,c.sa_ct4,
        c.komega_alpha,c.komega_beta0,c.komega_sigma_k,c.komega_sigma_w,
        c.komega_sigma_d0,c.komega_clim,c.sst_sigma_k1,c.sst_sigma_k2,c.sst_sigma_w1,c.sst_sigma_w2,
        c.sst_production_limiter,c.smagorinsky_Cs,c.wale_Cw,c.des_Cdes};
    for (double v : values)
        if (!std::isfinite(v)) throw std::invalid_argument("non-finite turbulence coefficient");
    if(c.C_mu<=0.0 || c.C1<0.0 || c.C2<0.0 || c.beta_star<=0.0 ||
       c.beta1<=0.0 || c.beta2<=0.0 || c.gamma1<=0.0 || c.gamma2<=0.0 || c.a1<=0.0 ||
       c.sigma_k<=0.0 || c.sigma_epsilon<=0.0 || c.rng_C_mu<=0.0 || c.rng_C1<=0.0 ||
       c.rng_C2<=0.0 || c.rng_sigma_k<=0.0 || c.rng_sigma_epsilon<=0.0 ||
       c.rng_eta0<=0.0 || c.rng_beta<=0.0 || c.realizable_A0<=0.0 ||
       c.realizable_C2<=0.0 || c.realizable_sigma_k<=0.0 ||
       c.realizable_sigma_epsilon<=0.0 || c.sa_cb1<=0.0 || c.sa_cb2<0.0 ||
       c.sa_sigma<=0.0 || c.sa_kappa<=0.0 || c.sa_cw2<0.0 || c.sa_cw3<=0.0 ||
       c.sa_cv1<=0.0 || c.sa_ct3<0.0 || c.sa_ct4<0.0 || c.komega_alpha<=0.0 ||
       c.komega_beta0<=0.0 || c.komega_sigma_k<=0.0 || c.komega_sigma_w<=0.0 ||
       c.komega_sigma_d0<0.0 || c.komega_clim<0.0 || c.sst_sigma_k1<=0.0 ||
       c.sst_sigma_k2<=0.0 || c.sst_sigma_w1<=0.0 || c.sst_sigma_w2<=0.0 ||
       c.sst_production_limiter<=0.0 || c.smagorinsky_Cs<0.0 || c.des_Cdes<=0.0)
        throw std::invalid_argument("invalid turbulence coefficients");
    validate_turbulence_correction_controls(c.corrections);
}

} // namespace cfdx::physics
