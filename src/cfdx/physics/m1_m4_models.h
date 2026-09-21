#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics::m1m4 {

constexpr double kappa = 0.41;
constexpr double sigma_sb = 5.670374419e-8;

inline void require_finite(double x, const char* name)
{
    if (!std::isfinite(x)) throw std::invalid_argument(std::string(name) + " must be finite");
}

inline double safe_positive(double x, double floor, const char* name)
{
    require_finite(x, name);
    if (x <= floor) throw std::invalid_argument(std::string(name) + " must be > floor");
    return x;
}

// ---------------- M1: pressure-velocity coupling ----------------

struct CouplingOptions {
    int outer_correctors = 1;
    int pressure_correctors = 2;
    int non_orthogonal_correctors = 0;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    bool consistent = false;
};

inline void validate_coupling(const CouplingOptions& c)
{
    if (c.outer_correctors < 1 || c.pressure_correctors < 1 ||
        c.non_orthogonal_correctors < 0)
        throw std::invalid_argument("invalid pressure-velocity corrector counts");
    if (!(c.velocity_relaxation > 0.0 && c.velocity_relaxation <= 1.0) ||
        !(c.pressure_relaxation > 0.0 && c.pressure_relaxation <= 1.0))
        throw std::invalid_argument("relaxation factors must be in (0,1]");
}

inline double under_relax(double old_value, double new_value, double alpha)
{
    if (!(alpha > 0.0 && alpha <= 1.0)) throw std::invalid_argument("alpha must be in (0,1]");
    return std::fma(alpha, new_value - old_value, old_value);
}

// Momentum predictor coefficient d = V/aP.
// The pressure correction is U' = -d grad(p') for the conventional
// discretisation A_P U_P = H_P - grad(p)_P V.
inline double pressure_velocity_coefficient(double volume, double aP)
{
    if (!(volume > 0.0) || !(aP > 0.0)) throw std::invalid_argument("volume and aP must be positive");
    return volume / aP;
}

// Rhie-Chow: pressure-gradient correction of the face mass flux.
// This is written in mass-flux form and therefore receives density explicitly.
inline double rhie_chow_mass_flux(double rho,
                                  double interpolated_velocity_flux,
                                  double pressure_owner,
                                  double pressure_neighbour,
                                  double grad_p_owner,
                                  double grad_p_neighbour,
                                  double d_f,
                                  double area)
{
    if (!(rho > 0.0) || !(d_f > 0.0) || !(area >= 0.0))
        throw std::invalid_argument("invalid Rhie-Chow input");
    const double grad_p_face = 0.5 * (grad_p_owner + grad_p_neighbour);
    const double dp_over_d = (pressure_neighbour - pressure_owner) / d_f;
    const double correction = rho * d_f * (dp_over_d - grad_p_face) * area;
    return interpolated_velocity_flux - correction;
}

// ---------------- M2: turbulence ----------------

struct KEpsilonCoefficients {
    double Cmu = 0.09;
    double C1 = 1.44;
    double C2 = 1.92;
    double sigma_k = 1.0;
    double sigma_epsilon = 1.3;
};

inline double k_epsilon_nut(double k, double epsilon, const KEpsilonCoefficients& c = {})
{
    if (k < 0.0 || epsilon <= 0.0 || c.Cmu <= 0.0)
        throw std::invalid_argument("k-epsilon requires k>=0 and epsilon>0");
    return c.Cmu * k * k / epsilon;
}

inline double k_epsilon_production(double nut, double strain_rate)
{
    if (nut < 0.0 || strain_rate < 0.0) throw std::invalid_argument("invalid k-epsilon production input");
    return 2.0 * nut * strain_rate * strain_rate;
}

inline double k_epsilon_dissipation_source(double rho, double C2, double epsilon, double k)
{
    if (rho <= 0.0 || C2 < 0.0 || epsilon < 0.0 || k <= 0.0)
        throw std::invalid_argument("invalid epsilon dissipation input");
    return rho * C2 * epsilon * epsilon / k;
}

struct SSTCoefficients {
    double betaStar = 0.09;
    double beta1 = 0.075;
    double beta2 = 0.0828;
    double gamma1 = 5.0 / 9.0;
    double gamma2 = 0.44;
    double a1 = 0.31;
};

inline double sst_limiter(double k, double omega, double strain_rate,
                          double a1 = 0.31)
{
    if (k < 0.0 || omega <= 0.0 || strain_rate < 0.0 || a1 <= 0.0)
        throw std::invalid_argument("invalid SST limiter input");
    return a1 * omega / std::max(a1 * omega, strain_rate);
}

inline double sst_nut(double k, double omega, double strain_rate,
                      double F2, const SSTCoefficients& c = {})
{
    if (k < 0.0 || omega <= 0.0 || strain_rate < 0.0 ||
        F2 < 0.0 || F2 > 1.0)
        throw std::invalid_argument("invalid SST input");
    const double denominator = std::max(c.a1 * omega, strain_rate * F2);
    return c.a1 * k / denominator;
}

inline double sst_blended(double F1, double value1, double value2)
{
    if (F1 < 0.0 || F1 > 1.0) throw std::invalid_argument("F1 outside [0,1]");
    return F1 * value1 + (1.0 - F1) * value2;
}

inline double van_driest_damping(double y_plus, double A_plus = 26.0)
{
    if (y_plus < 0.0 || A_plus <= 0.0) throw std::invalid_argument("invalid wall damping input");
    return 1.0 - std::exp(-y_plus / A_plus);
}

inline double wall_function_u_plus(double y_plus, double E = 9.793)
{
    if (y_plus <= 0.0 || E <= 0.0) throw std::invalid_argument("invalid wall-function input");
    if (y_plus < 11.0) return y_plus;
    return std::log(E * y_plus) / kappa;
}

inline double smagorinsky_nut(double delta, double strain_rate, double Cs = 0.17)
{
    if (delta <= 0.0 || strain_rate < 0.0 || Cs < 0.0)
        throw std::invalid_argument("invalid Smagorinsky input");
    const double l = Cs * delta;
    return l * l * strain_rate;
}

inline double des_length_scale(double delta, double wall_distance, double Cdes = 0.65)
{
    if (delta <= 0.0 || wall_distance <= 0.0 || Cdes <= 0.0)
        throw std::invalid_argument("invalid DES length-scale input");
    return std::min(wall_distance, Cdes * delta);
}

// ---------------- M3: thermal / CHT ----------------

inline double thermal_diffusivity(double conductivity, double rho, double cp)
{
    if (conductivity < 0.0 || rho <= 0.0 || cp <= 0.0)
        throw std::invalid_argument("invalid thermal properties");
    return conductivity / (rho * cp);
}

inline double conductive_flux(double conductivity, double T_owner,
                              double T_neighbour, double distance)
{
    if (conductivity < 0.0 || distance <= 0.0)
        throw std::invalid_argument("invalid conductive-flux input");
    return -conductivity * (T_neighbour - T_owner) / distance;
}

inline double interface_conductance(double k1, double k2, double d1, double d2,
                                    double area)
{
    if (k1 <= 0.0 || k2 <= 0.0 || d1 <= 0.0 || d2 <= 0.0 || area <= 0.0)
        throw std::invalid_argument("invalid interface conductance input");
    return area / (d1 / k1 + d2 / k2);
}

inline double interface_heat_flux(double conductance, double T1, double T2)
{
    if (conductance < 0.0) throw std::invalid_argument("conductance must be non-negative");
    return conductance * (T1 - T2);
}

inline double energy_source_linearization(double source, double dsource_dT,
                                          double T0, double& Su, double& Sp)
{
    require_finite(source, "source");
    require_finite(dsource_dT, "dsource_dT");
    require_finite(T0, "T0");
    Sp = std::min(0.0, dsource_dT);
    Su = source - Sp * T0;
    return Su + Sp * T0;
}

// ---------------- M4: radiation ----------------

inline double blackbody(double T)
{
    if (T < 0.0) throw std::invalid_argument("temperature must be non-negative");
    return sigma_sb * std::pow(T, 4);
}

inline double gray_emission(double emissivity, double T)
{
    if (emissivity < 0.0 || emissivity > 1.0 || T < 0.0)
        throw std::invalid_argument("invalid gray-body input");
    return emissivity * blackbody(T);
}

inline double two_surface_exchange(double e1, double e2, double T1, double T2,
                                   double F12)
{
    if (e1 <= 0.0 || e1 > 1.0 || e2 <= 0.0 || e2 > 1.0 ||
        T1 < 0.0 || T2 < 0.0 || F12 < 0.0 || F12 > 1.0)
        throw std::invalid_argument("invalid radiation exchange input");
    if (F12 == 0.0) return 0.0;
    const double resistance = (1.0-e1)/e1 + 1.0/F12 + (1.0-e2)/e2;
    return sigma_sb * (std::pow(T1,4)-std::pow(T2,4)) / resistance;
}

inline void validate_view_factors(const std::vector<double>& F, std::size_t n,
                                  double tol = 1e-10)
{
    if (F.size() != n*n) throw std::invalid_argument("view-factor matrix size mismatch");
    for (std::size_t i=0; i<n; ++i) {
        double sum = 0.0;
        for (std::size_t j=0; j<n; ++j) {
            const double f = F[i*n+j];
            if (f < -tol || f > 1.0+tol) throw std::invalid_argument("view factor outside [0,1]");
            sum += f;
        }
        if (std::abs(sum-1.0) > tol) throw std::invalid_argument("view-factor row does not close");
    }
}

inline void validate_reciprocity(const std::vector<double>& F,
                                 const std::vector<double>& areas,
                                 std::size_t n, double tol = 1e-10)
{
    if (F.size()!=n*n || areas.size()!=n) throw std::invalid_argument("view-factor dimensions mismatch");
    for (std::size_t i=0;i<n;++i)
        for (std::size_t j=0;j<n;++j)
            if (std::abs(areas[i]*F[i*n+j]-areas[j]*F[j*n+i]) > tol)
                throw std::invalid_argument("view-factor reciprocity violation");
}

inline double blackbody_intensity(double T)
{
    if (T < 0.0) throw std::invalid_argument("temperature must be non-negative");
    return blackbody(T) / M_PI;
}

inline double p1_source(double absorption, double mean_intensity, double T)
{
    if (absorption < 0.0 || mean_intensity < 0.0 || T < 0.0)
        throw std::invalid_argument("invalid P1 source input");
    // P1 source convention: S_r = 4*kappa*(I_b - J).
    return 4.0 * absorption * (blackbody_intensity(T) - mean_intensity);
}

struct Direction {
    double x, y, z, weight;
};

inline void validate_dom(const std::vector<Direction>& dirs, double tol = 1e-10)
{
    if (dirs.empty()) throw std::invalid_argument("DOM requires directions");
    double weight_sum=0.0;
    for (const auto& d: dirs) {
        const double n=std::sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
        if (std::abs(n-1.0)>tol || d.weight<0.0)
            throw std::invalid_argument("invalid DOM direction");
        weight_sum += d.weight;
    }
    if (std::abs(weight_sum-4.0*M_PI)>1e-8)
        throw std::invalid_argument("DOM weights must integrate 4*pi");
}

// ---------------- Analytical validation oracles ----------------

inline double plane_poiseuille_velocity(double y, double H, double dpdx, double mu)
{
    if (H <= 0.0 || mu <= 0.0) throw std::invalid_argument("invalid Poiseuille parameters");
    return -dpdx * y * (H-y) / (2.0*mu);
}

inline double plane_poiseuille_flow_rate_per_width(double H, double dpdx, double mu)
{
    if (H <= 0.0 || mu <= 0.0) throw std::invalid_argument("invalid Poiseuille parameters");
    return -dpdx * H*H*H / (12.0*mu);
}

inline double fully_developed_temperature(double x, double Tin, double qdot,
                                          double rho, double cp, double U, double A,
                                          double perimeter)
{
    if (rho<=0.0 || cp<=0.0 || U<=0.0 || A<=0.0 || perimeter<0.0)
        throw std::invalid_argument("invalid thermal validation parameters");
    return Tin + qdot * perimeter * x / (rho * U * A * cp);
}

} // namespace cfdx::physics::m1m4
