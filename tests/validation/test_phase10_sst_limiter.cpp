#include "cfdx/physics/turbulence.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>

int main() {
    try {
        const double k=0.1, omega=10.0, F2=0.8, strain=5.0;
        const double nut=cfdx::physics::turbulent_kinematic_viscosity_komega_sst(k,omega,F2,strain);
        const double expected=0.31*k/std::max(0.31*omega,strain*F2);
        if(std::abs(nut-expected)>1e-14) throw std::runtime_error("SST limiter mismatch");
        std::cout<<"PHASE10_SST_LIMITER: PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"PHASE10_SST_LIMITER: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
