#pragma once

// CFDX wall-function reference kernels are backend-independent.

#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cfdx::physics::wall {

struct WallFunctionConstants {
    double kappa = 0.41;
    double E = 9.8;
};

inline double y_plus(double wall_distance,double friction_velocity,double nu) {
    if(wall_distance<=0.0 || friction_velocity<0.0 || nu<=0.0) throw std::invalid_argument("y_plus: invalid inputs");
    return wall_distance*friction_velocity/nu;
}

inline double u_plus_log(double yplus,double kappa=0.41,double E=9.8) {
    if(yplus<=0.0 || kappa<=0.0 || E<=0.0) throw std::invalid_argument("u_plus_log: invalid inputs");
    return std::log(E*yplus)/kappa;
}

inline double log_layer_crossover(double kappa=0.41,double E=9.8) {
    if(kappa<=0.0 || E<=1.0) throw std::invalid_argument("log_layer_crossover: invalid inputs");
    double y=11.0;
    for(int it=0;it<200;++it) {
        const double next=std::log(E*y)/kappa;
        if(std::abs(next-y)<1e-13*y) return next;
        y=next;
    }
    throw std::runtime_error("log_layer_crossover: no convergence");
}

inline double turbulent_viscosity_log(double wall_distance,double friction_velocity,double nu,double kappa=0.41,double E=9.8) {
    const double yp=y_plus(wall_distance,friction_velocity,nu);
    if(yp<=log_layer_crossover(kappa,E)) return 0.0;
    const double up=u_plus_log(yp,kappa,E);
    if(up<=0.0) return 0.0;
    return std::max(0.0, wall_distance*friction_velocity/up-nu);
}

inline double omega_log(double wall_distance,double friction_velocity,double beta_star=0.09,double kappa=0.41) {
    if(wall_distance<=0.0 || friction_velocity<=0.0 || beta_star<=0.0 || kappa<=0.0)
        throw std::invalid_argument("omega_log: invalid inputs");
    return friction_velocity/(std::sqrt(beta_star)*kappa*wall_distance);
}

} // namespace cfdx::physics::wall
