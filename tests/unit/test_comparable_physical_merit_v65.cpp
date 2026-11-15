#include "cfdx/physics/comparable_physical_merit_v65.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool condition, const char* msg) {
    if(!condition) throw std::runtime_error(msg);
}
static void rejected(const std::function<void()>& op) {
    bool denied=false;
    try {op();}catch(const std::exception&){denied=true;}
    require(denied,"V65 negative test did not reject");
}
static PhysicalCandidateV51 state() {
    PhysicalCandidateV51 x; x.cells=2;
    x.fields={3,4,5,6,0,0,1,1,0.1,0.2,3,3};
    return x;
}
static PhysicalResidualV51 result(double momentum=0.10, double scale=0.1) {
    PhysicalResidualV51 r;
    r.signed_raw={momentum,0.,0.,0.,0.,0.,0.04,-0.02,0.02,-0.02,0.05,0.03};
    r.normalized=r.signed_raw;
    for(auto& x:r.normalized)x/=scale;
    r.face_flux={0.02,0.02,0.0};
    r.eddy_viscosity={0.001,0.001};
    r.flux_fixed_point_defect=1e-13;
    r.merit=10;
    r.flux_closed=true;
    return r;
}
int main() {
    const auto x=state();
    const auto old=result();
    auto lower=result(0.05);
    // The candidate may use a different own normalization, never compared by V65.
    for(auto& v:lower.normalized)v*=100.;
    auto e=evaluate_comparable_merit_v65(x,x,old,lower,1e-6,1e-7);
    require(e.accepted && e.merit_ratio<0.995,"fixed scale contraction");
    auto higher=result(0.2);
    e=evaluate_comparable_merit_v65(x,x,old,higher,1e-6,1e-7);
    require(!e.accepted && !e.joint_allowed,"growth must reject");
    e=evaluate_comparable_merit_v65(x,x,old,higher,1e-6,1e-7,true);
    require(e.accepted && e.baseline,"first-step baseline");
    // Zero-step means exactly the same F(X) (not a newly assembled rAU history).
    e=evaluate_comparable_merit_v65(x,x,old,old,1e-6,1e-7);
    require(e.merit_ratio==1. && !e.monotone_allowed,"identity zero-step");
    // Reject per-candidate normalization hiding actual physical growth.
    for(auto& v:higher.normalized)v*=1e-6;
    e=evaluate_comparable_merit_v65(x,x,old,higher,1e-6,1e-7);
    require(!e.accepted,"per-state normalization cannot hide growth");
    // Reconstructed Rhie-Chow must reproduce both the flux and raw continuity.
    const std::vector<std::size_t> owners={0,0,1};
    const std::vector<int> neighbours={-1,1,-1};
    const std::vector<double> phi={0.02,0.02,0.0};
    const auto flux=audit_comparable_flux_v65(x,old,phi,
        [&](std::size_t f){return owners[f];},
        [&](std::size_t f){return neighbours[f];},&phi);
    require(flux.continuity_linf<1e-15 && flux.reconstructed_relative_l1==0.,"RC closure");
    auto corrupt=phi;corrupt[1]+=1e-4;
    rejected([&]{audit_comparable_flux_v65(x,old,corrupt,
        [&](std::size_t f){return owners[f];},
        [&](std::size_t f){return neighbours[f];});});
    auto wrong=old;wrong.signed_raw[6]=0.1;
    rejected([&]{audit_comparable_flux_v65(x,wrong,phi,
        [&](std::size_t f){return owners[f];},
        [&](std::size_t f){return neighbours[f];});});
    auto no_closure=old;no_closure.flux_fixed_point_defect=2e-9;
    rejected([&]{evaluate_comparable_merit_v65(x,x,no_closure,lower,1e-6,1e-7);});
    rejected([&]{evaluate_comparable_merit_v65(x,x,old,lower,-1,1e-7);});
    std::cout<<"V65_SAME_OPERATOR_TEST_PASS\n";
}
