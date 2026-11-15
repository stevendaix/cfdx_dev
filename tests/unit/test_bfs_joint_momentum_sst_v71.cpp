// V71: falsify stale rAU, unclosed pressure, nonphysical SST changes,
// in-place mutation, and failure to recompute F_phys after each substep.
#include "cfdx/physics/bfs_joint_momentum_sst_v71.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>static void fails(F&& f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,"falsification did not reject");}
int main(){
    constexpr std::array<std::size_t,3> own={0,0,1};
    constexpr std::array<int,3> nb={-1,1,-1};
    const auto owner=[&](std::size_t f){return own.at(f);};
    const auto neighbour=[&](std::size_t f){return nb.at(f);};
    PhysicalCandidateV51 base;
    base.cells=2;base.fields.assign(12,0.);
    base.fields[0]=1.;base.fields[1]=1.;
    base.fields[8]=.3;base.fields[9]=.3;
    base.fields[10]=2.;base.fields[11]=2.;
    const auto flux=[](const PhysicalCandidateV51& s){
        return std::vector<double>{-1., .2 + (s.fields[6]-s.fields[7]) + .2*(s.fields[0]-s.fields[1]),
                                   .3+s.fields[7]};
    };
    std::size_t evaluations=0;
    const auto evaluate=[&](const PhysicalCandidateV51& s) {
        ++evaluations;
        PhysicalResidualV51 r;
        r.flux_closed=true;r.flux_fixed_point_defect=0.;
        r.face_flux=flux(s);r.signed_raw.assign(12,0.);r.normalized.assign(12,0.);
        r.eddy_viscosity.resize(2);r.inverse_momentum_diagonal.assign(6,.2);
        for(std::size_t i=0;i<2;++i){
            r.signed_raw[i]=s.fields[i]-(1.5-.1*s.fields[8]);
            r.signed_raw[8+i]=s.fields[8+i]-(.2+.02*s.fields[i]);
            r.signed_raw[10+i]=s.fields[10+i]-(2.+.1*s.fields[i]);
            r.eddy_viscosity[i]=s.fields[8+i]/s.fields[10+i];
        }
        r.signed_raw[6]=r.face_flux[0]+r.face_flux[1];
        r.signed_raw[7]=-r.face_flux[1]+r.face_flux[2];
        for(std::size_t i=0;i<12;++i){
            r.normalized[i]=r.signed_raw[i];
            r.merit=std::max(r.merit==std::numeric_limits<double>::infinity()?0.:r.merit,
                             std::abs(r.signed_raw[i]));
        }
        return r;
    };
    const auto start_f=evaluate(base);
    auto projected_calls=0,transport_calls=0;
    const auto project=[&](const PhysicalCandidateV51& candidate,
                           const PhysicalResidualV51& residual,double fraction){
        ++projected_calls;
        const auto response=[&](const std::vector<double>& delta){
            BfsPressureResponseV70 r;r.state=candidate;
            r.state.fields[6]+=delta[0];r.state.fields[7]+=delta[1];
            r.face_flux=flux(r.state);
            return r;
        };
        return bfs_project_pressure_v70(candidate,residual,response,
                                        owner,neighbour,true,0,1e-10,fraction);
    };
    const auto transport=[&](const PhysicalCandidateV51& candidate,
                             const std::vector<double>& f){
        ++transport_calls;check(f.size()==3,"transport got wrong flux");
        auto t=candidate;
        for(std::size_t i=0;i<2;++i){
            t.fields[8+i]=.2+.02*candidate.fields[i];
            t.fields[10+i]=2.+.1*candidate.fields[i];
        }
        return t;
    };
    const auto snapshot=base;
    const auto row=bfs_joint_private_cycle_v71(
        base,start_f,base,start_f,evaluate,project,transport,
        owner,neighbour,.25,1.,0.,0.,1e-8,1e-8);
    check(row.mass_after_pressure<1e-10,"pressure response did not conserve mass");
    check(row.final_mass_closed,"same-state mass after SST not conserved");
    check(row.momentum_update_linf>0.,"momentum predictor not exercised");
    check(projected_calls==1 && transport_calls==1 && evaluations>=4,
          "all three physical substeps must execute");
    check(row.omega_after<1e-12,"same flux SST transport did not solve omega");
    check(row.private_only && snapshot.fields==base.fields,"private cycle wrote to live state");
    auto stale=start_f;stale.inverse_momentum_diagonal.clear();
    fails([&]{(void)bfs_momentum_predict_v71(base,stale,.5);});
    fails([&]{(void)bfs_momentum_predict_v71(base,start_f,2.);});
    auto broken_projection=[&](const PhysicalCandidateV51& s,const PhysicalResidualV51&,double){
        BfsPressureProjectionV70 p;p.state=s;p.face_flux={1.,2.,3.};p.fully_projected=true;return p;
    };
    fails([&]{(void)bfs_joint_private_cycle_v71(base,start_f,base,start_f,
        evaluate,broken_projection,transport,owner,neighbour,.25,1.,0.,0.,1e-8,1e-8);});
    auto dishonest_project=[&](const PhysicalCandidateV51& s,
                                const PhysicalResidualV51& r,double frac){
        auto out=project(s,r,frac);
        out.state.fields[10]+=1.; // pressure response secretly changes omega
        return out;
    };
    fails([&]{(void)bfs_joint_private_cycle_v71(base,start_f,base,start_f,
        evaluate,dishonest_project,transport,owner,neighbour,.25,1.,0.,0.,1e-8,1e-8);});
    auto broken_transport=[&](const PhysicalCandidateV51& s,const std::vector<double>&){
        auto t=s;t.fields[0]+=1.;return t;
    };
    fails([&]{(void)bfs_joint_private_cycle_v71(base,start_f,base,start_f,
        evaluate,project,broken_transport,owner,neighbour,.25,1.,0.,0.,1e-8,1e-8);});
    std::cout<<"V71_JOINT_PRIVATE_PASS mass="<<row.mass_after_sst
             <<" omega="<<row.omega_after<<" momentum_step="<<row.momentum_update_linf
             <<" commit=0\n";
}
