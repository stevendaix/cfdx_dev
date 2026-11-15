@@ -0,0 +1,78 @@
#include "cfdx/physics/physical_remainder_rates_v61.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool b,const char* msg){if(!b)throw std::runtime_error(msg);}
static PhysicalResidualV51 physical(const PhysicalCandidateV51& x,bool ties){
    PhysicalResidualV51 r;
    const auto n=x.cells;
    r.signed_raw.assign(6*n,0.);
    r.normalized.assign(6*n,0.);
    r.eddy_viscosity.assign(n,1.);
    r.face_flux.assign(n+1,0.);
    CandidatePhysicalStagesV57 t;
    t.momentum_convection.assign(3*n,0.);
    t.inverse_momentum_diagonal.assign(3*n,1.);
    t.effective_dynamic_viscosity.assign(n,1.);
    t.strain_rate.assign(n,2.);
    t.sst_blend_f2.assign(n,0.5);
    t.wall_momentum_face_mu.assign(n+1,1.);
    t.wall_production.assign(n,1.);
    for(std::size_t i=0;i<n;++i){
        const auto ux=x.fields[i],uy=x.fields[n+i],k=x.fields[4*n+i],w=x.fields[5*n+i];
        for(std::size_t b=0;b<3;++b){
            const double u=x.fields[b*n+i];
            t.momentum_convection[b*n+i]=0.1*u;
            r.signed_raw[b*n+i]=u+.2*k+.1*w+t.momentum_convection[b*n+i];
        }
        if(ties && i==0)
            r.signed_raw[0]+=4.*std::max({ux-.1,uy-.2,0.});
        r.signed_raw[3*n+i]=ux-.1*uy;
        r.signed_raw[4*n+i]=2.*k+.1*ux;
        r.signed_raw[5*n+i]=2.*w+.1*k;
        t.inverse_momentum_diagonal[i]+=.01*k;
        t.effective_dynamic_viscosity[i]+=.01*w;
        r.face_flux[i]=.1*ux+.03*x.fields[3*n+i];
    }
    r.face_flux.back()=0.1*x.fields[1];
    r.stages_v57=t;
    r.normalized=r.signed_raw;
    r.merit=0.;for(const auto a:r.normalized)r.merit=std::max(r.merit,std::abs(a));
    r.flux_closed=true;r.flux_fixed_point_defect=0.;
    return r;
}
int main(){
    PhysicalCandidateV51 x;
    x.cells=2;x.fields={.1,.1,.2,.2,.1,.1,.1,.1,.7,.7,2.,2.};
    std::vector<double> v(x.fields.size(),0.),w(v);
    v[0]=1.;w[x.cells]=1.;
    v[4*x.cells]=0.1;w[5*x.cells]=0.2;
    const auto smooth=[&](const PhysicalCandidateV51& z){return physical(z,false);};
    const auto s=audit_physical_remainder_rates_v61(x,v,w,1e-4,smooth);
    require(s.evaluations==18 && s.all_scales_additive &&
            !s.persistent_defect_suspected && !s.roundoff_suspected,
            "V61 smooth multi-rate oracle rejected a linear six-field residual");
    const auto kink=[&](const PhysicalCandidateV51& z){return physical(z,true);};
    const auto k=audit_physical_remainder_rates_v61(x,v,w,1e-4,kink);
    require(!k.all_scales_additive && k.persistent_defect_suspected &&
            k.worst_block==0 &&
            k.maximum_physical_component[2]>1e-6,
            "V61 failed to track persistent momentum limiter non-additivity");
    bool threw=false;
    try { (void)audit_physical_remainder_rates_v61(x,v,w,0.,smooth); }
    catch(const std::invalid_argument&){threw=true;}
    require(threw,"V61 accepted zero finite-difference step");
    const auto forged=[&](const PhysicalCandidateV51& z){
        auto r=physical(z,false);r.flux_closed=false;return r;
    };
    threw=false;
    try { (void)audit_physical_remainder_rates_v61(x,v,w,1e-5,forged); }
    catch(const std::runtime_error&){threw=true;}
    require(threw,"V61 accepted a stale/invalid physical face closure");
    std::cout<<"V61_PHYSICAL_RATES_PASS evaluations="<<k.evaluations
             <<" coarse="<<k.maximum_physical_component[0]
             <<" quarter="<<k.maximum_physical_component[2]
             <<" newton_committed=0\n";
}
