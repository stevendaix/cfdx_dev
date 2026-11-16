#include "cfdx/physics/global_joint_jv_v40.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
using cfdx::physics::audit_joint_jv_v40;
int main() {
    const std::vector<double> x{1,2,3,4,5,6}; // Ux, Uy, Uz, p, k, omega
    const std::vector<double> v{1,2,1,0,2,3}; // last two: logarithmic direction
    const std::vector<double> scales(6,1.0);
    const auto physical=[](const std::vector<double>& z) {
        return std::vector<double>{z[0]+0.4*z[4],
                                   z[1]+0.1*z[5],z[2],
                                   z[3]+0.2*z[0],z[4]+2*z[5],z[5]+0.5*z[1]};
    };
    const auto r=audit_joint_jv_v40(x,v,scales,scales,1,physical,1e-12,1e-12,0);
    assert(r.smooth && r.residual_finite && r.function_evaluations==3);
    assert(std::abs(r.scaled_jv[0]-(1+0.4*5*2))<1e-7);
    assert(std::abs(r.scaled_jv[3]-0.2)<1e-7);
    assert(std::abs(r.scaled_jv[4]-(5*2+2*6*3))<1e-7);
    assert(std::abs(r.scaled_jv[5]-(6*3+0.5*2))<1e-7);
    const auto kink=[](const std::vector<double>& z) {
        return std::vector<double>{std::abs(z[0]-1.0),z[1],z[2],z[3],z[4],z[5]};
    };
    const auto singular=audit_joint_jv_v40(x,v,scales,scales,1,kink,1e-12,1e-12);
    assert(!singular.smooth && singular.max_relative_one_sided_gap>1.0);
    auto p_direction=v;p_direction[3]=1.0;
    try {(void)audit_joint_jv_v40(x,p_direction,scales,scales,1,physical,1e-12,1e-12,0);assert(false);}
    catch(const std::invalid_argument&) {}
    const auto incomplete=[](const std::vector<double>&) {return std::vector<double>(5);};
    try {(void)audit_joint_jv_v40(x,v,scales,scales,1,incomplete,1e-12,1e-12);assert(false);}
    catch(const std::runtime_error&) {}
    std::cout<<"V40_GLOBAL_JV_CONTRACT 7/7 PASS not_production_FV_qualified=1\n";
}
