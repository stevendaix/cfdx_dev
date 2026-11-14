#include "cfdx/physics/state_only_six_residual_v50.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using cfdx::physics::make_native_six_residual_v50;
static void require(bool yes) {
    if (!yes) throw std::runtime_error("V50 six-field residual test failed");
}
int main() {
    const std::array<std::vector<double>,3> momentum{{
        {1.0,-2.0},{-3.0,4.0},{0.0,0.0}}};
    const std::array<double,3> scale{{2.0,4.0,1.0}};
    const std::vector<double> div{0.1,-0.1};
    const std::vector<double> sst{0.6,-0.6,3.0,-2.0};
    const std::vector<double> row{2.0,2.0,5.0,4.0};
    const std::vector<double> field{3.0,3.0,2.0,2.0};
    const auto a=make_native_six_residual_v50(momentum,scale,div,0.5,
                                              sst,row,field,true);
    const auto b=make_native_six_residual_v50(momentum,scale,div,0.5,
                                              sst,row,field,true);
    require(!cfdx::physics::native_six_final_stop_eligible_v50(a,0.0,0.5));
    require(cfdx::physics::native_six_final_stop_eligible_v50(a,0.0,1.0));
    require(!cfdx::physics::native_six_final_stop_eligible_v50(a,1.01,1.0));
    require(a.momentum_global_relative==1.0);
    require(a.cells==2 && a.normalized.size()==12 &&
            a.fingerprint==b.fingerprint && a.normalized==b.normalized);
    require(a.normalized[0]==0.5 && a.normalized[1]==-1.0);
    require(a.normalized[2]==-0.75 && a.normalized[3]==1.0);
    require(a.normalized[6]==0.2 && a.normalized[7]==-0.2);
    require(std::abs(a.continuity_l1-0.4)<1e-15 &&
            std::abs(a.linf[4]-0.1)<1e-15 &&
            std::abs(a.linf[5]-0.3)<1e-15 &&
            a.worst_cell[1]==1 && a.merit==1.0);
    auto fails=[&](auto operation) {
        bool thrown=false;
        try {operation();} catch(const std::exception&) {thrown=true;}
        require(thrown);
    };
    fails([&] {make_native_six_residual_v50(momentum,scale,div,0.5,sst,row,field,false);});
    fails([&] {make_native_six_residual_v50(momentum,scale,div,0.,sst,row,field,true);});
    fails([&] {make_native_six_residual_v50(momentum,scale,div,0.5,{},row,field,true);});
    auto badmom=momentum;badmom[0][0]=std::numeric_limits<double>::quiet_NaN();
    fails([&] {make_native_six_residual_v50(badmom,scale,div,0.5,sst,row,field,true);});
    auto badrow=row;badrow[0]=0.0;
    fails([&] {make_native_six_residual_v50(momentum,scale,div,0.5,sst,badrow,field,true);});
    auto badscale=scale;badscale[1]=0.0;
    fails([&] {make_native_six_residual_v50(momentum,badscale,div,0.5,sst,row,field,true);});
    auto badsst=sst;badsst[0]=std::numeric_limits<double>::infinity();
    fails([&] {make_native_six_residual_v50(momentum,scale,div,0.5,badsst,row,field,true);});
    fails([&] {cfdx::physics::native_six_final_stop_eligible_v50(a,0.0,0.0);});
    require(!cfdx::physics::native_six_final_stop_eligible_v50(
        a,std::numeric_limits<double>::quiet_NaN(),1.0));
    // The digest must detect perturbations in *any* of the six blocks.
    for (std::size_t d=0;d<6;++d) {
        auto mod=momentum;
        auto cd=div;
        auto sd=sst;
        if (d<3) mod[d][0]+=0.0001;
        else if(d==3) cd[0]+=0.0001;
        else sd[(d-4)*2]+=0.0001;
        const auto candidate=make_native_six_residual_v50(
            mod,scale,cd,0.5,sd,row,field,true);
        require(candidate.fingerprint!=a.fingerprint);
    }
}
