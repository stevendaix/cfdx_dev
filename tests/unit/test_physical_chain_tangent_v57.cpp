#include "cfdx/physics/physical_chain_tangent_v57.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
using namespace cfdx::physics;
static void check(bool p,const char* msg){if(!p)throw std::runtime_error(msg);}
static PhysicalResidualV51 physical(const PhysicalCandidateV51& x){
    const auto n=x.cells;
    PhysicalResidualV51 r;
    r.merit=0.;
    r.signed_raw.assign(6*n,0.);
    r.normalized.assign(6*n,0.);
    r.face_flux={x.fields[0]+0.05*x.fields[4],x.fields[1]-0.03*x.fields[5]};
    r.eddy_viscosity.assign(n,0.);
    auto& t=r.stages_v57.emplace();
    t.inverse_momentum_diagonal.assign(3*n,0.);
    t.effective_dynamic_viscosity.assign(n,0.);
    t.strain_rate.assign(n,0.);
    t.sst_blend_f2.assign(n,0.);
    t.wall_momentum_face_mu={-1.,0.003};
    t.wall_production.assign(n,-1.);
    t.momentum_convection.assign(3*n,0.);
    for(std::size_t i=0;i<n;++i){
        const auto U=x.fields[i],k=x.fields[4*n+i],omega=x.fields[5*n+i];
        const auto mu=1.+0.01*k/omega;
        r.eddy_viscosity[i]=mu;
        t.effective_dynamic_viscosity[i]=mu;
        t.strain_rate[i]=std::sqrt(1.+U*U);
        t.sst_blend_f2[i]=1. / (1.+omega);
        t.wall_production[i]=k*U;
        for(std::size_t d=0;d<3;++d)
            t.inverse_momentum_diagonal[d*n+i]=1./(2.+0.1*mu+std::abs(r.face_flux[i]));
        for(std::size_t b=0;b<6;++b){
            const auto idx=b*n+i;
            r.signed_raw[idx]=x.fields[idx]+(b+1)*0.13*mu*U+
                    0.05*r.face_flux[i]*r.face_flux[i];
            r.normalized[idx]=r.signed_raw[idx];
            r.merit=std::max(r.merit,std::abs(r.normalized[idx]));
        }
    }
    r.flux_closed=true;r.flux_fixed_point_defect=1e-13;
    return r;
}
int main(){
    PhysicalCandidateV51 x; x.cells=2;
    x.fields={0.3,-0.4,0.1,0.2,0.,0.,0.3,0.5,1.3,1.6,2.2,2.8};
    std::vector<double> direction(x.fields.size(),0.);
    for(std::size_t i=0;i<direction.size();++i)
        direction[i]=0.05*std::sin(static_cast<double>(i+1));
    const auto a=audit_physical_chain_v57(x,direction,physical,1e-5);
    const auto repeat=audit_physical_chain_v57(x,direction,physical,1e-5);
    check(a.evaluations==5 && a.signed_raw.size()==12 && a.drau.size()==6 &&
          a.dphi.size()==2 && a.dmu.size()==2 && a.dwall_production.size()==2,
          "V57 missing candidate physical stage");
    check(a.signed_raw==repeat.signed_raw && a.dphi==repeat.dphi &&
          a.drau==repeat.drau,"V57 candidate re-evaluation is history-dependent");
    check(a.global_halving<1e-6 && v57_norm2(a.drau)>1e-6 &&
          v57_norm2(a.dmu)>1e-6 && v57_norm2(a.dwall_production)>1e-6,
          "V57 failed to differentiate RC/rAU/viscosity/SST/wall dependencies");
    check(a.wall_active_faces==1 && a.smooth_both_sides,
          "V57 malformed wall/branch diagnostic");
    bool rejected=false;
    try {
        auto no_trace=[](const PhysicalCandidateV51& x){auto r=physical(x);r.stages_v57.reset();return r;};
        (void)audit_physical_chain_v57(x,direction,no_trace,1e-5);
    }catch(const std::runtime_error&){rejected=true;}
    check(rejected,"V57 missing physical stage silently accepted");
    auto same=direction;for(auto& v:same)v=0.;
    rejected=false;
    try{(void)audit_physical_chain_v57(x,same,physical,1e-5);}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"V57 zero tangent silently accepted");
    auto kink=[](const PhysicalCandidateV51& x){
        auto r=physical(x);
        const double nonsmooth=std::abs(x.fields[0]-0.3);
        r.signed_raw[0]+=nonsmooth;
        r.normalized[0]=r.signed_raw[0];
        return r;
    };
    const auto kinked=audit_physical_chain_v57(x,direction,kink,1e-5);
    check(kinked.global_sided_nonsmoothness>0.01 &&
          !kinked.smooth_both_sides,
          "V57 non-smooth sided physics was wrongly certified by h-halving");
    std::cout<<"V57_CHAIN_PASS halving="<<a.global_halving
             <<" rc="<<v57_norm2(a.dphi)<<" rau="<<v57_norm2(a.drau)
             <<" mu="<<v57_norm2(a.dmu)<<" wall="<<v57_norm2(a.dwall_production)<<"\n";
}
