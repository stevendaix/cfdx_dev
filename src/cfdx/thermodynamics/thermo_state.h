#pragma once
#include <cmath>
#include <stdexcept>
namespace cfdx::thermodynamics {
struct ThermoState{double rho=1.0,mu=1e-5,cp=1000.0,conductivity=0.02,enthalpy=0.0,speed_of_sound=0.0;bool is_valid()const{return std::isfinite(rho)&&std::isfinite(mu)&&std::isfinite(cp)&&std::isfinite(conductivity)&&rho>0.0&&mu>=0.0&&cp>0.0&&conductivity>=0.0;}};
struct IdealGasThermoModel{double R=287.05,gamma=1.4,cp=1004.5,mu_ref=1.716e-5,T_ref=273.15,S=110.4,k_ref=0.0241;
ThermoState state(double p,double T)const{if(p<=0.0||T<=0.0)throw std::invalid_argument("invalid ideal-gas state");ThermoState s;s.rho=p/(R*T);s.cp=cp;s.mu=mu_ref*std::pow(T/T_ref,1.5)*(T_ref+S)/(T+S);s.conductivity=k_ref*std::pow(T/T_ref,0.76);s.enthalpy=cp*T;s.speed_of_sound=std::sqrt(gamma*R*T);return s;}};
}