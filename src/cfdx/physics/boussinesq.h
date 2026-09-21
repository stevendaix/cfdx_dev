#pragma once
#include <cmath>
#include <stdexcept>
namespace cfdx::physics {
struct BoussinesqModel{double rho_ref=1.0,beta=0.0,T_ref=300.0,gravity=9.81;
double density(double T)const{if(rho_ref<=0.0||beta<0.0)throw std::invalid_argument("invalid Boussinesq parameters");return rho_ref*(1.0-beta*(T-T_ref));}
double buoyancy_acceleration(double T)const{return gravity*beta*(T-T_ref);}};
}