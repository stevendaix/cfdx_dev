#include "cfdx/physics/turbulence.h"
#include "cfdx/physics/radiation.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
struct RefPoint { double y_plus; double u_plus; };
constexpr RefPoint channel_dns[]={{0.0,0.0},{0.053648,0.053639},{0.21456,0.21443},{0.48263,0.48197},{0.85771,0.85555},{1.3396,1.3339},{1.9279,1.9148},{2.6224,2.5939},{3.4226,3.3632},{4.3280,4.2095},{5.3381,5.1133},{6.4523,6.0493},{7.6700,6.9892},{8.9902,7.9052},{10.412,8.7741}};
double flat_plate_laminar_cf(double Re_x){if(Re_x<=0)throw std::invalid_argument("Re_x");return 0.664/std::sqrt(Re_x);}
double flat_plate_turbulent_cf(double Re_x){if(Re_x<=0)throw std::invalid_argument("Re_x");return 0.0592/std::pow(Re_x,0.2);}
}
int main(){try{
 const double k=0.81,eps=0.27;
 if(std::abs(turbulent_kinematic_viscosity_kepsilon(k,eps)-0.09*k*k/eps)>1e-14)throw std::runtime_error("k-epsilon closure");
 if(std::abs(turbulent_kinematic_viscosity_komega_sst(k,2.0,0.4)-0.31*k/std::max(0.31*2.0,2.0*0.4))>1e-14)throw std::runtime_error("SST closure");
 const Vec3 s{2,-1,0.5}; if(!(smagorinsky_eddy_viscosity(0.2,strain_rate_magnitude(s))>0))throw std::runtime_error("Smagorinsky");
 if(!(des_eddy_viscosity(0.2,0.01,1.0)>0))throw std::runtime_error("DES");
 for(const auto&p:channel_dns)if(!std::isfinite(p.y_plus)||!std::isfinite(p.u_plus))throw std::runtime_error("DNS reference");
 if(std::abs(flat_plate_laminar_cf(1e6)-0.000664)>1e-15)throw std::runtime_error("laminar flat plate");
 if(std::abs(flat_plate_turbulent_cf(1e6)-0.00373515)>1e-6)throw std::runtime_error("turbulent flat plate");
 const double qb=STEFAN_BOLTZMANN*(std::pow(900.,4)-std::pow(500.,4));
 if(std::abs(qb-two_surface_net_exchange(1,1,900,500,1))>1e-10*qb)throw std::runtime_error("radiation limit");
 std::cout<<"LEVEL_B_REFERENCE_BENCHMARKS: PASS\n";return 0;
}catch(const std::exception&e){std::cerr<<"LEVEL_B_REFERENCE_BENCHMARKS: FAIL: "<<e.what()<<"\n";return 1;}}
