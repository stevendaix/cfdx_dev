#pragma once
#include "cfdx/physics/thermophysical_models.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
namespace cfdx::physics {
enum class AdvancedTurbulenceModel { LAMINAR, KEPSILON, RNG_KEPSILON, REALIZABLE_KEPSILON, KOMEGA, SST, SPALART_ALLMARAS, SMAGORINSKY, WALE, DYNAMIC_KEQN, DES, DDES, IDDES };
enum class TurbulenceImplementationKind { CLOSURE, TRANSPORT_MODEL };
struct RealizableKEpsilonInvariants {
    double strain_magnitude=0.0;
    double rotation_magnitude=0.0;
    double third_invariant=0.0;
};

inline double realizable_kepsilon_cmu_from_invariants(
    const RealizableKEpsilonInvariants& inv,
    double k, double epsilon, double A0=4.04)
{
    if(!std::isfinite(inv.strain_magnitude) || !std::isfinite(inv.rotation_magnitude) ||
       !std::isfinite(inv.third_invariant) || inv.strain_magnitude<0.0 ||
       inv.rotation_magnitude<0.0 || k<0.0 || epsilon<=0.0 || A0<=0.0)
        throw std::invalid_argument("invalid realizable k-epsilon inputs");
    const double S2=std::max(inv.strain_magnitude*inv.strain_magnitude,1e-30);
    const double W=std::clamp(inv.third_invariant,-1.0/std::sqrt(6.0),1.0/std::sqrt(6.0));
    const double phi=std::acos(std::clamp(std::sqrt(6.0)*W,-1.0,1.0))/3.0;
    const double As=std::sqrt(6.0)*std::cos(phi);
    const double Ustar=std::sqrt(S2+inv.rotation_magnitude*inv.rotation_magnitude);
    const double denom=A0+As*Ustar*k/std::max(epsilon,1e-300);
    if(!(denom>0.0) || !std::isfinite(denom))
        throw std::domain_error("realizable k-epsilon Cmu denominator is invalid");
    return 1.0/denom;
}

inline double rng_kepsilon_c1_star(
    double eta, double C1=1.42, double eta0=4.38, double beta=0.012)
{
    if(!std::isfinite(eta) || eta<0.0 || C1<=0.0 || eta0<=0.0 || beta<=0.0)
        throw std::invalid_argument("invalid RNG k-epsilon inputs");
    return C1-eta*(1.0-eta/eta0)/(1.0+beta*eta*eta*eta);
}

struct SSTBlendedCoefficients {
    double sigma_k=0.0;
    double sigma_omega=0.0;
    double beta=0.0;
    double gamma=0.0;
};

inline SSTBlendedCoefficients sst_blended_coefficients(
    double F1,
    double sigma_k1, double sigma_k2,
    double sigma_w1, double sigma_w2,
    double beta1, double beta2,
    double gamma1, double gamma2)
{
    if(!std::isfinite(F1) || F1<0.0 || F1>1.0)
        throw std::invalid_argument("SST F1 must lie in [0,1]");
    const double f=std::clamp(F1,0.0,1.0);
    return {
        f*sigma_k1+(1.0-f)*sigma_k2,
        f*sigma_w1+(1.0-f)*sigma_w2,
        f*beta1+(1.0-f)*beta2,
        f*gamma1+(1.0-f)*gamma2
    };
}

enum class TurbulenceWallRegime { VISCOSITY_AFFECTED, BUFFER, LOG_LAYER };

inline TurbulenceWallRegime classify_wall_y_plus(double y_plus)
{
    if(!std::isfinite(y_plus) || y_plus<0.0)
        throw std::invalid_argument("wall y+ must be finite and non-negative");
    if(y_plus<=5.0) return TurbulenceWallRegime::VISCOSITY_AFFECTED;
    if(y_plus<30.0) return TurbulenceWallRegime::BUFFER;
    return TurbulenceWallRegime::LOG_LAYER;
}

struct TurbulenceModelDescriptor {
 AdvancedTurbulenceModel model=AdvancedTurbulenceModel::LAMINAR;
 TurbulenceImplementationKind implementation=TurbulenceImplementationKind::CLOSURE;
};
inline constexpr TurbulenceImplementationKind implementation_kind(AdvancedTurbulenceModel model){
 switch(model) {
 case AdvancedTurbulenceModel::KEPSILON:
 case AdvancedTurbulenceModel::RNG_KEPSILON:
 case AdvancedTurbulenceModel::KOMEGA:
 case AdvancedTurbulenceModel::SST:
 case AdvancedTurbulenceModel::SPALART_ALLMARAS:
     return TurbulenceImplementationKind::TRANSPORT_MODEL;
 case AdvancedTurbulenceModel::REALIZABLE_KEPSILON:
     // The Realizable k-epsilon kernel currently supplies algebraic Cmu only.
     // Do not advertise a transport equation until the k/epsilon solver is wired.
     return TurbulenceImplementationKind::CLOSURE;
 default:
     return TurbulenceImplementationKind::CLOSURE;
 }
}

inline constexpr bool has_transport_equation(AdvancedTurbulenceModel model)
{
 return implementation_kind(model)==TurbulenceImplementationKind::TRANSPORT_MODEL;
}
struct TurbulenceModelCoefficients {
 double C_mu=.09,C1=1.44,C2=1.92,sigma_k=1.0,sigma_epsilon=1.3,sigma_omega=.5,beta_star=.09,beta1=.075,beta2=.0828,gamma1=5.0/9.0,gamma2=.44,a1=.31,sigma_nu=2.0/3.0,Cb1=.1355,Cb2=.622,sigma_s=.3,Cs=.17,Cw=.3,CDES=.65,Cddes=.65;
};
inline void validate_turbulence_model_coefficients(const TurbulenceModelCoefficients& c){
 if(!std::isfinite(c.C_mu)||c.C_mu<=0||!std::isfinite(c.beta_star)||c.beta_star<=0||!std::isfinite(c.Cs)||c.Cs<0||!std::isfinite(c.CDES)||c.CDES<=0) throw std::invalid_argument("invalid turbulence model coefficients");
}
inline double k_epsilon_eddy_viscosity(double k,double epsilon,double Cmu=.09){if(k<0||epsilon<=0||Cmu<=0||!std::isfinite(k)||!std::isfinite(epsilon))throw std::invalid_argument("k-epsilon invalid inputs");return Cmu*k*k/epsilon;}
inline double rng_kepsilon_eddy_viscosity(double k,double epsilon,double Cmu=.0845){return k_epsilon_eddy_viscosity(k,epsilon,Cmu);}
inline double realizable_kepsilon_eddy_viscosity(double k,double epsilon,double Cmu=.09){
 if(k<0||epsilon<=0||Cmu<=0||!std::isfinite(k)||!std::isfinite(epsilon)) throw std::invalid_argument("realizable k-epsilon invalid inputs");
 return Cmu*k*k/epsilon;
}
inline double komega_eddy_viscosity(double k,double omega,double betaStar=.09){if(k<0||omega<=0||betaStar<=0||!std::isfinite(k)||!std::isfinite(omega))throw std::invalid_argument("k-omega invalid inputs");return k/omega;}
inline double sst_eddy_viscosity(double k,double omega,double strain,double a1=.31,double F2=1.0){
 if(k<0||omega<=0||strain<0||a1<=0||!std::isfinite(F2)||F2<0.0||F2>1.0)
     throw std::invalid_argument("SST invalid inputs");
 // Menter SST: the strain-rate branch is limited by F2 in the eddy-viscosity
 // denominator. Keeping F2 explicit avoids silently using the outer-layer form
 // in the near-wall branch.
 return a1*k/std::max(a1*omega,strain*F2);
}
enum class NegativeNuTildePolicy { CLAMP_ZERO, REJECT };
inline double spalart_allmaras_nu_t(double nu_tilde,double distance,double molecular_nu,
                                     NegativeNuTildePolicy policy=NegativeNuTildePolicy::CLAMP_ZERO){
 if(!std::isfinite(nu_tilde)||distance<=0||molecular_nu<=0) throw std::invalid_argument("Spalart-Allmaras invalid inputs");
 if(nu_tilde<0 && policy==NegativeNuTildePolicy::REJECT) throw std::domain_error("Spalart-Allmaras negative nu_tilde");
 const double nu_eff=std::max(0.0,nu_tilde);
 const double chi=nu_eff/molecular_nu,cv1=7.1,chi3=chi*chi*chi,fv1=chi3/(chi3+cv1*cv1*cv1);
 return fv1*nu_eff;
}
inline double smagorinsky_nu_t(double Cs,double delta,double strain){if(Cs<0||delta<=0||strain<0)throw std::invalid_argument("Smagorinsky invalid inputs");return (Cs*delta)*(Cs*delta)*strain;}
inline double wale_nu_t(double Cw,double delta,double strain){if(Cw<0||delta<=0||strain<0)throw std::invalid_argument("WALE invalid inputs");return (Cw*delta)*(Cw*delta)*strain;}
inline double des_length_scale(double wall_distance,double delta,double Cdes=.65){if(wall_distance<=0||delta<=0||Cdes<=0)throw std::invalid_argument("DES invalid inputs");return std::min(wall_distance,Cdes*delta);}
inline double ddes_shielding(double r_d){
 if(!std::isfinite(r_d)||r_d<0) throw std::invalid_argument("DDES shielding parameter must be non-negative");
 return 1.0-std::tanh(std::pow(8.0*r_d,3));
}
inline double ddes_length_scale_from_rd(double wall_distance,double delta,double Cdes,double r_d){
 if(wall_distance<=0||delta<=0||Cdes<=0||!std::isfinite(r_d)||r_d<0)
   throw std::invalid_argument("DDES length-scale inputs are invalid");
 const double fd=ddes_shielding(r_d);
 return wall_distance-fd*std::max(0.0,wall_distance-Cdes*delta);
}
inline double ddes_length_scale(double wall_distance,double delta,double Cdes=.65,double shielding=1){
 if(wall_distance<=0||delta<=0||Cdes<=0||shielding<0||shielding>1) throw std::invalid_argument("DDES invalid inputs");
 // shielding=1 means fully shielded (RANS); shielding=0 means unshielded LES branch.
 return wall_distance-shielding*std::max(0.0,wall_distance-Cdes*delta);
}
inline double iddes_shielding(double r_d,double stress_blend){
 if(!std::isfinite(r_d)||r_d<0||!std::isfinite(stress_blend)||stress_blend<0||stress_blend>1)
   throw std::invalid_argument("IDDES shielding inputs are invalid");
 const double fd=ddes_shielding(r_d);
 return std::clamp((1.0-stress_blend)+stress_blend*fd,0.0,1.0);
}
inline double iddes_length_scale_from_rd(double wall_distance,double delta,double Cdes,double r_d,double stress_blend){
 if(wall_distance<=0||delta<=0||Cdes<=0) throw std::invalid_argument("IDDES length-scale inputs are invalid");
 const double shielding=iddes_shielding(r_d,stress_blend);
 return wall_distance-shielding*std::max(0.0,wall_distance-Cdes*delta);
}
inline double iddes_length_scale(double wall_distance,double delta,double Cdes=.65,double shielding=1,double stress_blend=1){
 if(wall_distance<=0||delta<=0||Cdes<=0||shielding<0||shielding>1||stress_blend<0||stress_blend>1)
   throw std::invalid_argument("IDDES blend values must be in [0,1]");
 const double les=Cdes*delta;
 const double fd=std::clamp(shielding,0.0,1.0);
 const double fb=std::clamp(stress_blend,0.0,1.0);
 return (1-fb)*wall_distance+fb*std::min(wall_distance,fd*les);
}
}