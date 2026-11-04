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
    // Regression: a large upstream omega can mask a badly solved near-wall
    // equation if every cell is divided by the GLOBAL maximum omega.
    ScalarEquation eq;
    eq.matrix=SparseMatrix(3,3);
    for(std::size_t i=0;i<3;++i) eq.matrix.push_back(i,i,1.0);
    eq.matrix.finalize();
    eq.rhs=Vector(3,0.0);
    eq.rhs(0)=1.e6; eq.rhs(1)=1.0; eq.rhs(2)=1.0;
    Vector x(3,0.0);
    x(0)=1.e6; x(1)=1.1; x(2)=1.0;
    const double global_solution_defect=sst_equation_solution_defect(eq,x);
    // The V11 global certificate sees only ~1e-7 and would pass 1e-6.
    if (!(global_solution_defect / 1.e6 < 1.e-6)) return 4;
    const auto local=sst_local_equation_defect(eq,x,1.e-12);
    if (!(local.maximum_relative > 0.09 && local.maximum_relative < 0.091) ||
        local.worst_cell != 1) return 5;
    // Same physical system in units scaled by 60: the relative local defect
    // must not depend on the omega time-unit choice.
    for(std::size_t i=0;i<3;++i){eq.rhs(i)*=60.0; x(i)*=60.0;}
    const auto rescaled=sst_local_equation_defect(eq,x,60.e-12);
    if(std::abs(rescaled.maximum_relative-local.maximum_relative)>1.e-13) return 6;
    std::cout << "PASS V12 SST global and local residuals with dimensional invariance\n";
}
