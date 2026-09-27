#pragma once
#include "cfdx/physics/thermophysical_models.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <limits>
namespace cfdx::physics {
enum class AdvancedTurbulenceModel { LAMINAR, KEPSILON, RNG_KEPSILON, REALIZABLE_KEPSILON, KOMEGA, SST, SPALART_ALLMARAS, SMAGORINSKY, WALE, DYNAMIC_KEQN, DES, DDES, IDDES };
enum class TurbulenceImplementationKind { CLOSURE, TRANSPORT_MODEL };
enum class TurbulenceImplementationStatus { SOLVER_READY, KERNEL_ONLY, PLANNED };
enum class TurbulenceFamily { LAMINAR, RANS, LES, HYBRID };
enum class TurbulenceWallTreatment { RESOLVED, WALL_FUNCTION, ALL_Y_PLUS, NONE };
enum class TurbulenceCorrection {
    NONE,
    ROTATION_CURVATURE,
    COMPRESSIBILITY,
    ROUGHNESS,
    PRODUCTION_LIMITER,
    KATO_LAUNDER,
    QCR
};

struct TurbulenceCorrectionControls {
    bool rotation_curvature=false;
    bool compressibility=false;
    bool roughness=false;
    bool production_limiter=true;
    bool kato_launder=false;
    bool qcr=false;
    double rotation_coefficient=1.0;
    double curvature_coefficient=1.0;
    double compressibility_coefficient=1.0;
    double roughness_height=0.0;
    double roughness_coefficient=0.0;
    double production_limit=10.0;
    double kato_coefficient=1.0;
    double qcr_coefficient=0.3;
};

inline void validate_turbulence_correction_controls(const TurbulenceCorrectionControls& c)
{
    const double v[] = {c.rotation_coefficient,c.curvature_coefficient,
        c.compressibility_coefficient,c.roughness_height,c.roughness_coefficient,
        c.production_limit,c.kato_coefficient,c.qcr_coefficient};
    for(double x : v) if(!std::isfinite(x)) throw std::invalid_argument("non-finite turbulence correction");
    if(c.rotation_coefficient<0 || c.curvature_coefficient<0 ||
       c.compressibility_coefficient<0 || c.roughness_height<0 ||
       c.roughness_coefficient<0 || c.production_limit<=0 ||
       c.kato_coefficient<0 || c.qcr_coefficient<0)
        throw std::invalid_argument("invalid turbulence correction controls");
}

struct TurbulenceInvariants {
    double strain=0.0;
    double rotation=0.0;
    double curvature=0.0;
    double divergence=0.0;
    double wall_distance=0.0;
};

inline bool valid_turbulence_invariants(const TurbulenceInvariants& x)
{
    return std::isfinite(x.strain) && std::isfinite(x.rotation) &&
           std::isfinite(x.curvature) && std::isfinite(x.divergence) &&
           std::isfinite(x.wall_distance) &&
           x.strain>=0 && x.rotation>=0 && x.wall_distance>=0;
}

inline double rotation_curvature_factor(const TurbulenceInvariants& x,
                                        const TurbulenceCorrectionControls& c)
{
    if(!valid_turbulence_invariants(x)) throw std::invalid_argument("invalid rotation/curvature invariants");
    const double denom=std::max(x.strain*x.strain,1e-30);
    const double ratio=(x.rotation*x.rotation + c.curvature_coefficient*x.curvature*x.curvature)/denom;
    return std::clamp(1.0 + c.rotation_coefficient*ratio,0.0,10.0);
}

inline double compressibility_factor(const TurbulenceInvariants& x,
                                     const TurbulenceCorrectionControls& c,
                                     double speed_of_sound)
{
    if(!valid_turbulence_invariants(x) || speed_of_sound<=0 || !std::isfinite(speed_of_sound))
        throw std::invalid_argument("invalid compressibility inputs");
    const double mach=x.strain/std::max(speed_of_sound,1e-30);
    return 1.0/(1.0 + c.compressibility_coefficient*mach*mach);
}

inline double roughness_factor(double roughness_height,double wall_distance,
                               const TurbulenceCorrectionControls& c)
{
    if(roughness_height<0 || wall_distance<=0 || !std::isfinite(roughness_height) ||
       !std::isfinite(wall_distance))
        throw std::invalid_argument("invalid roughness inputs");
    const double ks=roughness_height/std::max(wall_distance,1e-30);
    return 1.0 + c.roughness_coefficient*ks;
}

inline double kato_launder_production_factor(const TurbulenceInvariants& x,
                                             const TurbulenceCorrectionControls& c)
{
    if(!valid_turbulence_invariants(x)) throw std::invalid_argument("invalid Kato-Launder invariants");
    return c.kato_coefficient * x.strain * x.rotation;
}

inline double qcr_stress_factor(double strain,double rotation,
                                const TurbulenceCorrectionControls& c)
{
    if(strain<0 || rotation<0 || !std::isfinite(strain) || !std::isfinite(rotation))
        throw std::invalid_argument("invalid QCR invariants");
    return std::clamp(1.0 + c.qcr_coefficient*rotation/std::max(strain,1e-30),0.0,10.0);
}

inline double corrected_turbulence_production(double production,
                                              const TurbulenceInvariants& x,
                                              const TurbulenceCorrectionControls& c,
                                              double speed_of_sound=1.0)
{
    if(!std::isfinite(production) || production<0 || !valid_turbulence_invariants(x))
        throw std::invalid_argument("invalid turbulence production");
    double p=production;
    if(c.rotation_curvature) p*=rotation_curvature_factor(x,c);
    if(c.compressibility) p*=compressibility_factor(x,c,speed_of_sound);
    if(c.roughness) p*=roughness_factor(c.roughness_height,x.wall_distance,c);
    if(c.kato_launder) p+=kato_launder_production_factor(x,c);
    if(c.production_limiter) p=std::min(p,c.production_limit);
    return std::max(0.0,p);
}

struct WALEInvariants {
    double S2=0.0;
    double Sd2=0.0;
    double Delta=0.0;
};

inline double wale_eddy_viscosity(const WALEInvariants& x,double Cw=0.325)
{
    if(!std::isfinite(x.S2)||!std::isfinite(x.Sd2)||x.S2<0||x.Sd2<0||x.Delta<=0||Cw<0)
        throw std::invalid_argument("invalid WALE invariants");
    const double denom=std::pow(x.S2,2.5)+std::pow(x.Sd2,1.25);
    if(denom<=0) return 0.0;
    return std::pow(Cw*x.Delta,2.0)*std::pow(x.Sd2,1.5)/denom;
}

inline double dynamic_les_coefficient(double resolved_stress,double test_stress,
                                      double denominator,double cmin=0.0,double cmax=0.23)
{
    if(!std::isfinite(resolved_stress)||!std::isfinite(test_stress)||
       !std::isfinite(denominator)||denominator<=0||cmin<0||cmax<cmin)
        throw std::invalid_argument("invalid dynamic LES inputs");
    return std::clamp((resolved_stress-test_stress)/denominator,cmin,cmax);
}

struct TurbulenceCorrectionCapabilities {
    bool rotation_curvature=false;
    bool compressibility=false;
    bool roughness=false;
    bool production_limiter=false;
    bool kato_launder=false;
    bool qcr=false;
};

struct TurbulenceCapability {
    AdvancedTurbulenceModel model=AdvancedTurbulenceModel::LAMINAR;
    const char* key="LAMINAR";
    const char* label="Laminar";
    TurbulenceFamily family=TurbulenceFamily::LAMINAR;
    TurbulenceImplementationStatus status=TurbulenceImplementationStatus::SOLVER_READY;
    TurbulenceImplementationKind implementation=TurbulenceImplementationKind::CLOSURE;
    int transported_equations=0;
    bool supports_steady=true;
    bool supports_transient=true;
    bool supports_2d=true;
    bool supports_3d=true;
    const char* required_fields="";
    const char* missing="";
    TurbulenceCorrectionCapabilities corrections{};
};
inline constexpr std::array<TurbulenceCapability,13> turbulence_capabilities{{
    {AdvancedTurbulenceModel::LAMINAR,"LAMINAR","Laminar",TurbulenceFamily::LAMINAR,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::CLOSURE,0,true,true,true,true,"U","turbulence transport, wall treatment and turbulence V&V are not applicable"},
    {AdvancedTurbulenceModel::KEPSILON,"KEPSILON","k-epsilon",TurbulenceFamily::RANS,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::TRANSPORT_MODEL,2,true,true,true,true,"k,epsilon","model-specific wall treatment, turbulence inlet specification and quantitative benchmark campaign"},
    {AdvancedTurbulenceModel::RNG_KEPSILON,"RNG_KEPSILON","RNG k-epsilon",TurbulenceFamily::RANS,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::TRANSPORT_MODEL,2,true,true,true,true,"k,epsilon","wall treatment contract and dedicated quantitative RNG benchmark campaign"},
    {AdvancedTurbulenceModel::REALIZABLE_KEPSILON,"REALIZABLE_KEPSILON","Realizable k-epsilon",TurbulenceFamily::RANS,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::TRANSPORT_MODEL,2,true,true,true,true,"k,epsilon","full mesh/y+ sensitivity and quantitative benchmark qualification"},
    {AdvancedTurbulenceModel::KOMEGA,"KOMEGA","k-omega",TurbulenceFamily::RANS,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::TRANSPORT_MODEL,2,true,true,true,true,"k,omega","wall treatment contract and quantitative benchmark qualification"},
    {AdvancedTurbulenceModel::SST,"SST","k-omega SST",TurbulenceFamily::RANS,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::TRANSPORT_MODEL,2,true,true,true,true,"k,omega","transition/corrections are not yet exposed; complete wall-treatment and quantitative benchmark qualification"},
    {AdvancedTurbulenceModel::SPALART_ALLMARAS,"SPALART_ALLMARAS","Spalart-Allmaras",TurbulenceFamily::RANS,TurbulenceImplementationStatus::SOLVER_READY,TurbulenceImplementationKind::TRANSPORT_MODEL,1,true,true,true,true,"nu_tilde","negative/rotation/compressibility/QCR variants and quantitative benchmark qualification"},
    {AdvancedTurbulenceModel::SMAGORINSKY,"SMAGORINSKY","Smagorinsky LES",TurbulenceFamily::LES,TurbulenceImplementationStatus::KERNEL_ONLY,TurbulenceImplementationKind::CLOSURE,0,false,true,false,true,"velocity gradient,Delta","LES transient driver, SGS BC contract, wall treatment and V&V"},
    {AdvancedTurbulenceModel::WALE,"WALE","WALE LES",TurbulenceFamily::LES,TurbulenceImplementationStatus::KERNEL_ONLY,TurbulenceImplementationKind::CLOSURE,0,false,true,false,true,"velocity-gradient tensor,Delta","full tensor WALE kernel integration, LES transient driver, wall treatment and V&V"},
    {AdvancedTurbulenceModel::DYNAMIC_KEQN,"DYNAMIC_KEQN","Dynamic one-equation LES",TurbulenceFamily::LES,TurbulenceImplementationStatus::PLANNED,TurbulenceImplementationKind::CLOSURE,1,false,true,false,true,"k_sgs,velocity gradient,Delta","complete dynamic SGS transport, test filtering, clipping, LES driver and V&V"},
    {AdvancedTurbulenceModel::DES,"DES","DES",TurbulenceFamily::HYBRID,TurbulenceImplementationStatus::KERNEL_ONLY,TurbulenceImplementationKind::CLOSURE,0,false,true,false,true,"base RANS fields,Delta,wall distance","RANS base-model coupling, shielding/length-scale integration, transient 3-D driver and V&V"},
    {AdvancedTurbulenceModel::DDES,"DDES","DDES",TurbulenceFamily::HYBRID,TurbulenceImplementationStatus::KERNEL_ONLY,TurbulenceImplementationKind::CLOSURE,0,false,true,false,true,"base RANS fields,Delta,wall distance,shielding","base-model coupling, complete shielding integration, wall treatment, transient 3-D driver and V&V"},
    {AdvancedTurbulenceModel::IDDES,"IDDES","IDDES",TurbulenceFamily::HYBRID,TurbulenceImplementationStatus::KERNEL_ONLY,TurbulenceImplementationKind::CLOSURE,0,false,true,false,true,"base RANS fields,Delta,wall distance,shielding","IDDES stress/wake shielding and blending integration, wall treatment, transient 3-D driver and V&V"}
}};
inline constexpr const TurbulenceCapability& turbulence_capability(AdvancedTurbulenceModel model) {
    for (const auto& c : turbulence_capabilities) if (c.model==model) return c;
    return turbulence_capabilities[0];
}
inline constexpr bool turbulence_correction_supported(
    AdvancedTurbulenceModel model, TurbulenceCorrection correction)
{
    const auto& c=turbulence_capability(model).corrections;
    switch(correction) {
    case TurbulenceCorrection::NONE: return true;
    case TurbulenceCorrection::ROTATION_CURVATURE: return c.rotation_curvature;
    case TurbulenceCorrection::COMPRESSIBILITY: return c.compressibility;
    case TurbulenceCorrection::ROUGHNESS: return c.roughness;
    case TurbulenceCorrection::PRODUCTION_LIMITER: return c.production_limiter;
    case TurbulenceCorrection::KATO_LAUNDER: return c.kato_launder;
    case TurbulenceCorrection::QCR: return c.qcr;
    }
    return false;
}

inline constexpr bool turbulence_solver_ready(AdvancedTurbulenceModel model) {
    return turbulence_capability(model).status==TurbulenceImplementationStatus::SOLVER_READY;
}
struct RealizableKEpsilonInvariants {
    double strain_magnitude=0.0;
    double rotation_magnitude=0.0;
    double third_invariant=0.0;
};

inline bool valid_realizable_kepsilon_inputs(
    const RealizableKEpsilonInvariants& inv,
    double k, double epsilon, double A0=4.04)
{
    return std::isfinite(inv.strain_magnitude) &&
           std::isfinite(inv.rotation_magnitude) &&
           std::isfinite(inv.third_invariant) &&
           std::isfinite(k) && std::isfinite(epsilon) && std::isfinite(A0) &&
           inv.strain_magnitude>=0.0 && inv.rotation_magnitude>=0.0 &&
           k>=0.0 && epsilon>0.0 && A0>0.0;
}

namespace detail {
inline double realizable_kepsilon_cmu_kernel(
    const RealizableKEpsilonInvariants& inv,
    double k, double epsilon, double A0)
{
    const double S2=std::max(inv.strain_magnitude*inv.strain_magnitude,1e-30);
    const double W=std::clamp(inv.third_invariant,-1.0/std::sqrt(6.0),1.0/std::sqrt(6.0));
    const double phi=std::acos(std::sqrt(6.0)*W)/3.0;
    const double As=std::sqrt(6.0)*std::cos(phi);
    const double Ustar=std::sqrt(S2+inv.rotation_magnitude*inv.rotation_magnitude);
    return A0+As*Ustar*k/std::max(epsilon,1e-300);
}
} // namespace detail

inline double realizable_kepsilon_cmu_from_invariants_unchecked(
    const RealizableKEpsilonInvariants& inv,
    double k, double epsilon, double A0=4.04)
{
    assert(valid_realizable_kepsilon_inputs(inv,k,epsilon,A0));
    const double denom=detail::realizable_kepsilon_cmu_kernel(inv,k,epsilon,A0);
    assert(denom>0.0 && std::isfinite(denom));
    return (denom>0.0 && std::isfinite(denom)) ? 1.0/denom : 0.0;
}

inline double realizable_kepsilon_cmu_from_invariants(
    const RealizableKEpsilonInvariants& inv,
    double k, double epsilon, double A0=4.04)
{
    if(!valid_realizable_kepsilon_inputs(inv,k,epsilon,A0))
        throw std::invalid_argument("invalid realizable k-epsilon inputs");
    const double denom=detail::realizable_kepsilon_cmu_kernel(inv,k,epsilon,A0);
    if(!(denom>0.0) || !std::isfinite(denom))
        throw std::domain_error("realizable k-epsilon Cmu denominator is invalid");
    return 1.0/denom;
}

inline double rng_kepsilon_c1_star(
    double eta, double C1=1.42, double eta0=4.38, double beta=0.012)
{
    if(!std::isfinite(eta) || !std::isfinite(C1) || !std::isfinite(eta0) ||
       !std::isfinite(beta) || eta<0.0 || C1<=0.0 || eta0<=0.0 || beta<=0.0)
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
    const std::array<double,8> coefficients{
        sigma_k1,sigma_k2,sigma_w1,sigma_w2,beta1,beta2,gamma1,gamma2};
    if(!std::isfinite(F1) || F1<0.0 || F1>1.0 ||
       !std::all_of(coefficients.begin(),coefficients.end(),
                    [](double value){return std::isfinite(value) && value>0.0;}))
        throw std::invalid_argument("SST blending inputs must be finite and positive");
    const double f=std::clamp(F1,0.0,1.0);
    return {
        f*sigma_k1+(1.0-f)*sigma_k2,
        f*sigma_w1+(1.0-f)*sigma_w2,
        f*beta1+(1.0-f)*beta2,
        f*gamma1+(1.0-f)*gamma2
    };
}

enum class TurbulenceWallRegime { VISCOSITY_AFFECTED, BUFFER, LOG_LAYER };

inline bool valid_wall_y_plus(double y_plus)
{
    return std::isfinite(y_plus) && y_plus>=0.0;
}

inline TurbulenceWallRegime classify_wall_y_plus_unchecked(double y_plus)
{
    assert(valid_wall_y_plus(y_plus));
    if(y_plus<=5.0) return TurbulenceWallRegime::VISCOSITY_AFFECTED;
    if(y_plus<30.0) return TurbulenceWallRegime::BUFFER;
    return TurbulenceWallRegime::LOG_LAYER;
}

inline TurbulenceWallRegime classify_wall_y_plus(double y_plus)
{
    if(!valid_wall_y_plus(y_plus))
        throw std::invalid_argument("wall y+ must be finite and non-negative");
    return classify_wall_y_plus_unchecked(y_plus);
}

struct TurbulenceModelDescriptor {
 AdvancedTurbulenceModel model=AdvancedTurbulenceModel::LAMINAR;
 TurbulenceImplementationKind implementation=TurbulenceImplementationKind::CLOSURE;
 TurbulenceImplementationStatus status=TurbulenceImplementationStatus::SOLVER_READY;
};
inline constexpr TurbulenceModelDescriptor turbulence_model_descriptor(AdvancedTurbulenceModel model) {
 const auto& c=turbulence_capability(model);
 return {model,c.implementation,c.status};
}
inline constexpr TurbulenceImplementationKind implementation_kind(AdvancedTurbulenceModel model){
 switch(model) {
 case AdvancedTurbulenceModel::KEPSILON:
 case AdvancedTurbulenceModel::RNG_KEPSILON:
 case AdvancedTurbulenceModel::REALIZABLE_KEPSILON:
 case AdvancedTurbulenceModel::KOMEGA:
 case AdvancedTurbulenceModel::SST:
 case AdvancedTurbulenceModel::SPALART_ALLMARAS:
     return TurbulenceImplementationKind::TRANSPORT_MODEL;
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
inline double k_epsilon_eddy_viscosity(double k,double epsilon,double Cmu=.09){if(k<0||epsilon<=0||Cmu<=0||!std::isfinite(k)||!std::isfinite(epsilon)||!std::isfinite(Cmu))throw std::invalid_argument("k-epsilon invalid inputs");return Cmu*k*k/epsilon;}
inline double rng_kepsilon_eddy_viscosity(double k,double epsilon,double Cmu=.0845){return k_epsilon_eddy_viscosity(k,epsilon,Cmu);}
inline double realizable_kepsilon_eddy_viscosity(double k,double epsilon,double Cmu){
 if(k<0||epsilon<=0||Cmu<=0||!std::isfinite(k)||!std::isfinite(epsilon)||!std::isfinite(Cmu)) throw std::invalid_argument("realizable k-epsilon invalid inputs");
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
