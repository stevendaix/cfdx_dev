#include "cfdx/physics/sst_solver.h"
#include <cmath>
#include <iostream>
#include <limits>
int main() {
    using namespace cfdx::physics;
    using namespace cfdx::core;
    Field<double,Location::CELL> oldk(1,"k_old","m2/s2",1), nextk(1,"k_new","m2/s2",1);
    Field<double,Location::CELL> oldw(1,"w_old","1/s",1), nextw(1,"w_new","1/s",1);
    oldk(0)=1.e-3; nextk(0)=1.0000001e-3;
    oldw(0)=1.e5; nextw(0)=1.0000001e5;
    const auto kd=sst_field_update_diagnostic(oldk,nextk,1.e-12);
    const auto wd=sst_field_update_diagnostic(oldw,nextw,1.e-12);
    const auto pass=sst_relative_convergence_certificate(kd,wd,1.e-12,1.e-4,2.e-7);
    if (!pass.converged || pass.k_equation_relative > 2.e-7 || pass.omega_equation_relative > 2.e-7) return 1;
    const auto miss=sst_relative_convergence_certificate(kd,wd,1.e-7,1.e-4,2.e-7);
    if (miss.converged || miss.k_equation_relative < 1.e-5) return 2;
    // Changing time units for omega cannot change its relative diagnostics.
    oldw(0)*=60; nextw(0)*=60;
    const auto wmin=sst_field_update_diagnostic(oldw,nextw,60.e-12);
    const auto consistent=sst_relative_convergence_certificate(kd,wmin,1.e-12,60.e-4,2.e-7);
    if (!consistent.converged ||
        std::abs(consistent.omega_equation_relative - pass.omega_equation_relative)>1.e-15) return 3;
    std::cout << "PASS V11 SST relative update plus true equation defect certificate\n";
}
