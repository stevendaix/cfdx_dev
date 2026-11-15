#include "cfdx/physics/bfs_sst_flux_picard_v69.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool ok,const char* error){if(!ok)throw std::runtime_error(error);}
template<class Fn>static void must_throw(Fn&& fn){bool failed=false;try{fn();}catch(const std::exception&){failed=true;}check(failed,"expected fail-closed rejection");}
static PhysicalResidualV51 sample(const PhysicalCandidateV51& x){
    PhysicalResidualV51 f;
    f.flux_closed=true;f.flux_fixed_point_defect=1e-14;
    f.signed_raw.assign(6*x.cells,0.);
    f.normalized.assign(6*x.cells,0.);
    f.signed_raw[3*x.cells]=.25;
    f.signed_raw[4*x.cells]=x.fields[4*x.cells]-1.;
    f.signed_raw[5*x.cells]=x.fields[5*x.cells]-2.;
    f.normalized=f.signed_raw;
    f.merit=std::max({std::abs(f.signed_raw[3*x.cells]),
                      std::abs(f.signed_raw[4*x.cells]),
                      std::abs(f.signed_raw[5*x.cells])});
    f.face_flux={x.fields[4*x.cells]-.1};
    f.eddy_viscosity={x.fields[5*x.cells]};
    return f;
}
int main(){
    PhysicalCandidateV51 state;
    state.cells=1;state.fields={3.,0.,0.,0.,2.,4.};
    const auto original=state;
    const auto physical=sample(state);
    std::size_t transports=0;
    const auto correct=[&](const PhysicalCandidateV51& x,const std::vector<double>& phi){
        check(phi.size()==1 && std::isfinite(phi[0]),"missing RC flux");
        ++transports;
        auto updated=x;
        updated.fields[4]=1.+.5*(x.fields[4]-1.);
        updated.fields[5]=2.+.5*(x.fields[5]-2.);
        return updated;
    };
    const auto run=probe_bfs_sst_flux_picard_v69(state,physical,sample,correct,0.,0.,2);
    check(transports==2 && run.size()==2,"not two fresh SST solves");
    check(run[0].omega_before==2. && run[0].omega_after==1.,"first omega contraction");
    check(run[1].omega_after==.5 && run[1].omega_ratio==.5,"second omega contraction");
    check(run[0].merit_ratio==.5 && run[0].rc_closed,"physical merit/RC monitor");
    check(state.fields==original.fields,"state leaked from read-only experiment");
    const auto gated=probe_bfs_sst_flux_picard_v69(
        state,physical,sample,correct,0.,0.,2,
        &state,&physical,1e-8,1e-8);
    check(gated[0].v65_would_pass && gated[0].v65_joint_ok &&
          gated[0].v65_monotone_ok && gated[0].v65_merit_ratio==.5,
          "V69 must reuse the original frozen-scale V65 gates");
    must_throw([&]{probe_bfs_sst_flux_picard_v69(
        state,physical,sample,correct,0.,0.,1,&state,nullptr);});
    must_throw([&]{probe_bfs_sst_flux_picard_v69(state,physical,sample,correct,0.,0.,0);});
    must_throw([&]{probe_bfs_sst_flux_picard_v69(state,physical,sample,
        [](auto x,const auto&){x.fields[0]+=1.;return x;},0.,0.,1);});
    must_throw([&]{probe_bfs_sst_flux_picard_v69(state,physical,sample,
        [](auto x,const auto&){x.fields[5]=-1.;return x;},0.,0.,1);});
    must_throw([&]{probe_bfs_sst_flux_picard_v69(state,physical,sample,
        [](auto x,const auto&){x.fields[4]=NAN;return x;},0.,0.,1);});
    must_throw([&]{probe_bfs_sst_flux_picard_v69(state,physical,
        [](const auto& x){auto r=sample(x);r.flux_closed=false;return r;},
        correct,0.,0.,1);});
    std::size_t evaluations=0;
    must_throw([&]{probe_bfs_sst_flux_picard_v69(state,physical,
        [&](const auto& value){auto f=sample(value);f.face_flux[0]+=.1*evaluations++;return f;},
        correct,0.,0.,1);});
    std::cout<<"V69_SST_RC_PICARD_PROBE_PASS cycles="<<run.size()
             <<" omega_after="<<run.back().omega_after<<"\n";
}