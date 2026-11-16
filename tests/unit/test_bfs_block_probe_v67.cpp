#include "cfdx/physics/bfs_block_probe_v67.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
template<class F>static void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,"V67 must fail closed");}
static PhysicalCandidateV51 make_state(){
    PhysicalCandidateV51 x;x.cells=2;
    x.fields={3,4,0,0,0,0,1,1,0.4,0.6,4,6};return x;
}
static PhysicalResidualV51 physical(const PhysicalCandidateV51& x){
    PhysicalResidualV51 r;
    r.signed_raw.assign(6*x.cells,0.0);
    const double omega=x.fields[5*x.cells];
    const double k=x.fields[4*x.cells];
    // A smooth, fully deterministic manufactured residual with omega
    // coupling to momentum and k: not a replacement for native SST.
    r.signed_raw[0]=1.0+0.5*(omega-4.0);
    r.signed_raw[4*x.cells]=0.1*(k-0.4);
    r.signed_raw[5*x.cells]=0.3*(omega-4.0);
    r.normalized=r.signed_raw;
    r.face_flux={0.01,0.01};
    r.eddy_viscosity.assign(x.cells,0.01);
    r.inverse_momentum_diagonal.assign(3*x.cells,0.2);
    r.flux_closed=true;r.flux_fixed_point_defect=0.;r.merit=1.0;
    return r;
}
int main(){
    const auto old=make_state();auto trial=old;
    trial.fields[0]=4.0;
    trial.fields[4*old.cells]=0.8;
    trial.fields[5*old.cells]=6.0;
    const auto frozen=mixed_sst_candidate_v67(old,trial,{0.,0.},0.,0.);
    for(std::size_t b=0;b<4;++b)
        check(frozen.fields[b*old.cells]==trial.fields[b*old.cells],"carrier changed");
    check(frozen.fields[4*old.cells]==old.fields[4*old.cells] &&
          frozen.fields[5*old.cells]==old.fields[5*old.cells],"SST freeze");
    const auto whole=mixed_sst_candidate_v67(old,trial,{1.,1.},0.,0.);
    check(whole.fields==trial.fields,"full SST trial must be identical");
    const auto half=mixed_sst_candidate_v67(old,trial,{0.5,0.25},0.,0.);
    check(std::abs(half.fields[4*old.cells]-0.6)<1e-14 &&
          std::abs(half.fields[5*old.cells]-4.5)<1e-14,"independent damping");
    const auto six=mixed_six_candidate_v67(old,trial,0.5,0.,0.);
    check(std::abs(six.fields[0]-(old.fields[0]+trial.fields[0])/2)<1e-14,
          "full-six carrier interpolation");
    check(std::abs(six.fields[5*old.cells]-5.0)<1e-14,
          "full-six omega interpolation");
    rejects([&]{mixed_six_candidate_v67(old,trial,0.,0.,0.);});
    rejects([&]{mixed_six_candidate_v67(old,trial,1.01,0.,0.);});
    const auto r0=physical(old);
    const auto outcome=probe_sst_block_mix_v67(old,trial,r0,{0.5,0.25},
          physical,1.e-5,1.e-5,0.,0.);
    check(outcome.merit_ratio>0. && std::isfinite(outcome.merit_ratio),"merit");
    const auto six_merit=probe_full_six_line_v67(
        old,trial,r0,0.5,physical,1.e-5,1.e-5,0.,0.);
    check(std::isfinite(six_merit.merit_ratio),"full-six physical merit");
    check(outcome.would_pass_original_guards ==
          (outcome.joint_allowed && outcome.monotone_allowed),"original gates");
    check(std::abs(v67_armijo_reference_ratio(1.)-.995)<1e-15,"full step match");
    check(std::abs(v67_armijo_reference_ratio(1./1024.)-
          (1.-.005/1024.))<1e-15,"short step slope");
    rejects([&]{v67_armijo_reference_ratio(0.);});
    rejects([&]{v67_armijo_reference_ratio(1.1);});
    rejects([&]{mixed_sst_candidate_v67(old,trial,{1.1,1.},0.,0.);});
    rejects([&]{mixed_sst_candidate_v67(old,trial,{NAN,1.},0.,0.);});
    rejects([&]{auto x=trial;x.cells=3;mixed_sst_candidate_v67(old,x,{1.,1.},0.,0.);});
    rejects([&]{auto x=trial;x.fields[5*old.cells]=-1.;mixed_sst_candidate_v67(old,x,{1.,1.},0.,0.);});
    rejects([&]{auto bad=r0;bad.flux_closed=false;probe_sst_block_mix_v67(
                  old,trial,bad,BfsBlockMixV67{1.,1.},physical,1.e-5,1.e-5,0.,0.);});
    std::cout<<"V67_BLOCK_COUNTERFACTUAL_PASS\n";
}