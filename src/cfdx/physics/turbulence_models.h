#pragma once
#include "cfdx/physics/thermophysical_models.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfdx::physics {
enum class AdvancedTurbulenceModel { LAMINAR, KEPSILON, RNG_KEPSILON, REALIZABLE_KEPSILON, KOMEGA, SST, SPALART_ALLMARAS, SMAGORINSKY, WALE, DYNAMIC_KEQN, DES, DDES, IDDES };
struct TurbulenceModelCoefficients {
 double C_mu=.09,C1=1.44,C2=1.92,sigma_k=1.0,sigma_epsilon=1.3,sigma_omega=.5,beta_star=.09,beta1=.075,beta2=.0828,gamma1=5.0/9.0,gamma2=.44,a1=.31,sigma_nu=2.0/3.0,Cb1=.1355,Cb2=.622,sigma_s=.3,Cs=.17,Cw=.3,CDES=.65,Cddes=.65;
};
inline void validate_turbulence_model_coefficients(const TurbulenceModelCoefficients& c){
 if(!std::isfinite(c.C_mu)||c.C_mu<=0||!std::isfinite(c.beta_star)||c.beta_star<=0||!std::isfinite(c.Cs)||c.Cs<0||!std::isfinite(c.CDES)||c.CDES<=0) throw std::invalid_argument("invalid turbulence model coefficients");
}
inline double k_epsilon_eddy_viscosity(double k,double epsilon,double Cmu=.09){if(k<0||epsilon<=0||Cmu<=0||!std::isfinite(k)||!std::isfinite(epsilon))throw std::invalid_argument("k-epsilon invalid inputs");return Cmu*k*k/epsilon;}
inline double rng_kepsilon_eddy_viscosity(double k,double epsilon,double Cmu=.0845){return k_epsilon_eddy_viscosity(k,epsilon,Cmu);}
inline double realizable_kepsilon_eddy_viscosity(double k,double epsilon,double Cmu=.09){return k_epsilon_eddy_viscosity(k,epsilon,Cmu);}
inline double komega_eddy_viscosity(double k,double omega,double betaStar=.09){if(k<0||omega<=0||betaStar<=0||!std::isfinite(k)||!std::isfinite(omega))throw std::invalid_argument("k-omega invalid inputs");return k/omega;}
inline double sst_eddy_viscosity(double k,double omega,double strain,double a1=.31){if(k<0||omega<=0||strain<0||a1<=0)throw std::invalid_argument("SST invalid inputs");return a1*k/std::max(a1*omega,strain);}
inline double spalart_allmaras_nu_t(double nu_tilde,double distance,double molecular_nu){if(nu_tilde<0||distance<=0||molecular_nu<=0)throw std::invalid_argument("Spalart-Allmaras invalid inputs");const double chi=nu_tilde/molecular_nu,cv1=7.1,chi3=chi*chi*chi,fv1=chi3/(chi3+cv1*cv1*cv1);return fv1*nu_tilde;}
inline double smagorinsky_nu_t(double Cs,double delta,double strain){if(Cs<0||delta<=0||strain<0)throw std::invalid_argument("Smagorinsky invalid inputs");return (Cs*delta)*(Cs*delta)*strain;}
inline double wale_nu_t(double Cw,double delta,double strain){if(Cw<0||delta<=0||strain<0)throw std::invalid_argument("WALE invalid inputs");return (Cw*delta)*(Cw*delta)*strain;}
inline double des_length_scale(double wall_distance,double delta,double Cdes=.65){if(wall_distance<=0||delta<=0||Cdes<=0)throw std::invalid_argument("DES invalid inputs");return std::min(wall_distance,Cdes*delta);}
inline double ddes_length_scale(double wall_distance,double delta,double Cdes=.65,double shielding=1){if(shielding<0||shielding>1)throw std::invalid_argument("DDES shielding must be in [0,1]");return std::min(wall_distance,Cdes*delta*std::max(shielding,1e-12));}
inline double iddes_length_scale(double wall_distance,double delta,double Cdes=.65,double shielding=1,double stress_blend=1){if(shielding<0||shielding>1||stress_blend<0||stress_blend>1)throw std::invalid_argument("IDDES blend values must be in [0,1]");const double les=Cdes*delta,fd=std::clamp(shielding,0.0,1.0),fb=std::clamp(stress_blend,0.0,1.0);return (1-fb)*wall_distance+fb*std::min(wall_distance,fd*les);}
}