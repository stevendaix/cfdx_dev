#pragma once
#include "cfdx/thermodynamics/thermo_state.h"
#include <array>
#include <cmath>
#include <stdexcept>
namespace cfdx::physics {
struct CompressibleState{double rho=1.0,u=0.0,v=0.0,w=0.0,p=101325.0,T=300.0;};
struct ConservativeState{double rho=0.0,rhou=0.0,rhov=0.0,rhow=0.0,rhoE=0.0;};
inline ConservativeState primitive_to_conservative(const CompressibleState&s,const thermodynamics::IdealGasThermoModel&eos){
 if(s.rho<=0.0||s.p<=0.0||s.T<=0.0)throw std::invalid_argument("invalid compressible state");
 const double e=eos.cp*s.T-s.p/s.rho,kinetic=0.5*(s.u*s.u+s.v*s.v+s.w*s.w);
 return {s.rho,s.rho*s.u,s.rho*s.v,s.rho*s.w,s.rho*(e+kinetic)};
}
inline std::array<double,5> rusanov_flux(const CompressibleState&L,const CompressibleState&R,const std::array<double,3>&n,const thermodynamics::IdealGasThermoModel&eos){
 const auto UL=primitive_to_conservative(L,eos),UR=primitive_to_conservative(R,eos);
 const double unL=L.u*n[0]+L.v*n[1]+L.w*n[2],unR=R.u*n[0]+R.v*n[1]+R.w*n[2];
 const double aL=eos.state(L.p,L.T).speed_of_sound,aR=eos.state(R.p,R.T).speed_of_sound;
 const double a=std::max(std::abs(unL)+aL,std::abs(unR)+aR);
 const std::array<double,5> FL{L.rho*unL,L.rho*L.u*unL+L.p*n[0],L.rho*L.v*unL+L.p*n[1],L.rho*L.w*unL+L.p*n[2],(UL.rhoE+L.p)*unL};
 const std::array<double,5> FR{R.rho*unR,R.rho*R.u*unR+R.p*n[0],R.rho*R.v*unR+R.p*n[1],R.rho*R.w*unR+R.p*n[2],(UR.rhoE+R.p)*unR};
 const std::array<double,5> qL{UL.rho,UL.rhou,UL.rhov,UL.rhow,UL.rhoE},qR{UR.rho,UR.rhou,UR.rhov,UR.rhow,UR.rhoE};
 std::array<double,5> F{};for(std::size_t i=0;i<5;++i)F[i]=0.5*(FL[i]+FR[i])-0.5*a*(qR[i]-qL[i]);return F;
}
} // namespace cfdx::physics
