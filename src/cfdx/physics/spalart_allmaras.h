#pragma once
#include <cmath>
#include <stdexcept>
namespace cfdx::physics {
struct SpalartAllmarasModel{double cb1=0.1355,cb2=0.622,kappa=0.41,cw2=0.3,cw3=2.0,cv1=7.1;
double fv1(double chi)const{double c3=chi*chi*chi;return c3/(c3+cv1*cv1*cv1);}
double turbulent_viscosity(double rho,double nu_tilde,double nu)const{if(rho<=0.0||nu<=0.0||nu_tilde<0.0)throw std::invalid_argument("invalid SA state");return rho*nu_tilde*fv1(nu_tilde/nu);}
double destruction_coefficient(double r)const{double g=r+cw2*(r*r*r-r);return cw3*g/(std::pow(g,4.0)+std::pow(cw3,4.0));}};
}