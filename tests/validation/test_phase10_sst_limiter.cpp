#include "cfdx/physics/turbulence.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        const double k=0.1, omega=10.0, F2=0.8, strain=5.0, a1=0.31;
        const double nut=cfdx::physics::turbulent_kinematic_viscosity_komega_sst(
            k,omega,F2,strain,a1);
        const double expected=a1*k/std::max(a1*omega,strain*F2);
        if(std::abs(nut-expected)>1e-14) throw std::runtime_error("SST limiter mismatch");
        std::cout<<"PHASE10_SST_LIMITER: PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"PHASE10_SST_LIMITER: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
