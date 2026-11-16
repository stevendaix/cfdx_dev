#include "cfdx/physics/physical_remainder_audit_v60.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
static PhysicalResidualV51 actual(const PhysicalCandidateV51& x,bool corner=false){
    const auto n=x.cells;
    PhysicalResidualV51 r;
    r.signed_raw.assign(6*n,0.);
    r.normalized.assign(6*n,0.);
    r.eddy_viscosity.assign(n,0.);
    r.face_flux.assign(n+1,0.);
    CandidatePhysicalStagesV57 tape;
    tape.inverse_momentum_diagonal.resize(3*n);
    tape.effective_dynamic_viscosity.resize(n);
    tape.strain_rate.resize(n);
    tape.sst_blend_f2.resize(n);
    tape.wall_production.resize(n);
    tape.wall_momentum_face_mu.assign(n+1,-1.);
    tape.momentum_convection.resize(3*n);
    for(std::size_t i=0;i<n;++i){
        const auto get=[&](std::size_t b){return x.fields[b*n+i];};
        const double ux=get(0),uy=get(1),p=get(3),k=get(4),w=get(5);
        for(std::size_t b=0;b<3;++b) {
            const double u=get(b);
            tape.momentum_convection[b*n+i]=0.04*u+0.015*k;
            tape.inverse_momentum_diagonal[b*n+i]=0.7+0.01*u+0.002*w;
            r.signed_raw[b*n+i]=tape.momentum_convection[b*n+i]+u+0.03*p+0.01*w;
        }
        // A nonsmooth point at the base (two active states + a zero arm).
        if(corner)r.signed_raw[0]+=4.*std::max({ux-.1,uy-.2,0.});
        r.signed_raw[3*n+i]=0.07*ux-0.02*uy+0.03*k+0.001*p;
        r.signed_raw[4*n+i]=3.*k+0.2*ux+0.1*w;
        r.signed_raw[5*n+i]=2.*w+0.3*k+0.03*uy;
        r.eddy_viscosity[i]=1.+0.02*k+0.01*w;
        tape.effective_dynamic_viscosity[i]=1.3+0.02*k;
        tape.strain_rate[i]=1.4+0.07*ux+0.03*uy;
        tape.sst_blend_f2[i]=0.6+0.05*k;
        tape.wall_production[i]=0.4+0.01*w;
        r.face_flux[i]=0.2*ux+0.1*p+0.02*k;
    }
    r.face_flux[n]=0.03*x.fields[0];
    tape.wall_momentum_face_mu[0]=1.2+0.01*x.fields[0];
    r.normalized=r.signed_raw;
    r.merit=0.;for(double a:r.normalized)r.merit=std::max(r.merit,std::abs(a));
    r.flux_fixed_point_defect=0.;r.flux_closed=true;
    r.stages_v57=std::move(tape);
    return r;
}
int main(){
    PhysicalCandidateV51 x;
    x.cells=4;x.fields.assign(6*x.cells,0.);
    for(std::size_t i=0;i<x.cells;++i){
        x.fields[i]=.1;x.fields[x.cells+i]=.2;
        x.fields[2*x.cells+i]=.3;
        x.fields[3*x.cells+i]=.1;
        x.fields[4*x.cells+i]=.7;
        x.fields[5*x.cells+i]=2.;
    }
    std::vector<double> v(x.fields.size(),0.),w(v);
    v[0]=1.;w[x.cells]=1.;
    v[4*x.cells+2]=0.04;w[5*x.cells+1]=0.06;
    const auto linear=[&](const PhysicalCandidateV51& y){return actual(y);};
    const auto passed=audit_physical_remainder_additivity_v60(x,v,w,1e-5,linear);
    require(passed.audited && passed.evaluations==6,"V60 did not reassemble 6 physical probes");
    require(passed.global_remainder<1e-9 && passed.rhie_chow_flux<1e-9 &&
            passed.inverse_diagonal<1e-9 && passed.sst_k<1e-9 &&
            passed.sst_omega<1e-9,"V60 smooth candidate should be additive");
    const auto sharp=[&](const PhysicalCandidateV51& y){return actual(y,true);};
    const auto blocked=audit_physical_remainder_additivity_v60(x,v,w,1e-5,sharp);
    require(blocked.momentum_remainder>0.15 && blocked.global_remainder>0.1,
            "V60 failed to detect nonsmooth remainder (third max branch)");
    require(blocked.convection<1e-9 && blocked.rhie_chow_flux<1e-9,
            "V60 incorrectly blamed smooth convection or Rhie-Chow");
    bool refused=false;
    try{
        const auto forged=[&](const PhysicalCandidateV51& y){
            auto r=actual(y);r.flux_closed=false;return r;};
        (void)audit_physical_remainder_additivity_v60(x,v,w,1e-5,forged);
    }catch(const std::runtime_error&){refused=true;}
    require(refused,"V60 accepted an unclosed physical flux");
    refused=false;
    try{
        const auto no_tape=[&](const PhysicalCandidateV51& y){
            auto r=actual(y);r.stages_v57.reset();return r;};
        (void)audit_physical_remainder_additivity_v60(x,v,w,1e-5,no_tape);
    }catch(const std::invalid_argument&){refused=true;}
    require(refused,"V60 accepted absent native stage tape");
    std::cout<<"V60_REMAINDER_AUDIT_PASS smooth="<<passed.global_remainder
             <<" nonsmooth="<<blocked.global_remainder
             <<" rc="<<blocked.rhie_chow_flux<<" evaluations="<<blocked.evaluations<<'\n';
}
