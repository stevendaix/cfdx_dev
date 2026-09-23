#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

struct SpalartAllmarasModel {
    double cb1=0.1355, cb2=0.622, sigma=2.0/3.0, kappa=0.41;
    double cw2=0.3, cw3=2.0, cv1=7.1;
    double ct3=1.2, ct4=0.5, cw1=0.0;

    SpalartAllmarasModel()
        : cw1(cb1/(kappa*kappa) + (1.0+cb2)/sigma) {}

    static double positive(double x,const char* name) {
        if(!std::isfinite(x)||x<0.0) throw std::invalid_argument(name);
        return x;
    }

    double fv1(double chi) const {
        chi=positive(chi,"SA chi must be finite and non-negative");
        const double chi3=chi*chi*chi;
        return chi3/(chi3+cv1*cv1*cv1);
    }

    double fv2(double chi) const {
        chi=positive(chi,"SA chi must be finite and non-negative");
        return 1.0-chi/(1.0+chi*fv1(chi));
    }

    double ft2(double chi) const {
        chi=positive(chi,"SA chi must be finite and non-negative");
        return ct3*std::exp(-ct4*chi*chi);
    }

    double turbulent_viscosity(double rho,double nu_tilde,double nu) const {
        if(!std::isfinite(rho)||!std::isfinite(nu)||rho<=0.0||nu<=0.0)
            throw std::invalid_argument("invalid SA fluid state");
        nu_tilde=positive(nu_tilde,"SA nu_tilde must be finite and non-negative");
        return rho*nu_tilde*fv1(nu_tilde/nu);
    }

    double wall_r(double nu_tilde,double wall_distance,double s_hat) const {
        if(!std::isfinite(nu_tilde)||nu_tilde<0.0 ||
           !std::isfinite(wall_distance)||wall_distance<=0.0 ||
           !std::isfinite(s_hat)||s_hat<=0.0)
            throw std::invalid_argument("invalid SA wall state");
        return nu_tilde/(s_hat*kappa*kappa*wall_distance*wall_distance);
    }

    double destruction_coefficient(double r) const {
        if(!std::isfinite(r)||r<0.0)
            throw std::invalid_argument("invalid SA destruction ratio");
        const double g=r+cw2*(std::pow(r,3.0)-r);
        const double g6=std::pow(g,6.0);
        return g*(1.0+std::pow(cw3,6.0))/(g6+std::pow(cw3,6.0));
    }

    double production_coefficient(double chi,double vorticity) const {
        if(!std::isfinite(vorticity)||vorticity<0.0)
            throw std::invalid_argument("invalid SA production state");
        return cb1*(1.0-ft2(chi))*vorticity;
    }
};

} // namespace cfdx::physics
