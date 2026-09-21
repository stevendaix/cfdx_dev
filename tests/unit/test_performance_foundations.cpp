#include "cfdx/core/linalg/solver_workspace.h"
#include "cfdx/core/linalg/mixed_precision.h"
#include "cfdx/physics/adaptive_cfl.h"
#include "cfdx/physics/boussinesq.h"
#include "cfdx/physics/spalart_allmaras.h"
#include "cfdx/thermodynamics/thermo_state.h"
#include <cmath>
#include <iostream>
int main(){using namespace cfdx;core::KrylovWorkspace k;k.resize(64);core::GmresWorkspace g;g.resize(64,8);if(k.r.size()!=64||g.basis.size()!=64*9)return 1;core::Vector a(3,1.0),b(3,2.0);if(std::abs(core::mixed_precision_dot(a,b)-6.0)>1e-12)return 2;if(!(physics::adaptive_time_step(1.0,2.0)<1.0))return 3;physics::BoussinesqModel bx{1.0,0.01,300.0,9.81};if(std::abs(bx.density(310.0)-0.9)>1e-12)return 4;physics::SpalartAllmarasModel sa;if(!(sa.turbulent_viscosity(1.0,1e-4,1e-5)>0.0))return 5;thermodynamics::IdealGasThermoModel gas;if(!gas.state(101325.0,300.0).is_valid())return 6;std::cout<<"performance/physics foundations: PASS\n";return 0;}