#include "cfdx/physics/hybrid_physical_jv_v57.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static void check(bool good,const char* what){if(!good)throw std::runtime_error(what);}
static PhysicalResidualV51 evaluate(const PhysicalCandidateV51& x){
    PhysicalResidualV51 r;
    r.merit=0.;r.flux_closed=true;r.flux_fixed_point_defect=1e-13;
    const double u=x.fields[0],p=x.fields[3],k=x.fields[4],w=x.fields[5];
    const double phi=u+0.1*k+0.02*p;
    r.face_flux={phi};r.eddy_viscosity={1.+0.01*k/w};
    r.signed_raw={phi*u+std::exp(0.1*p)-0.2*w,
        p+u*k,0.5*w,phi,k+0.3*u*u,w+0.4*k*p};
    r.normalized=r.signed_raw;
    for(double v:r.signed_raw)r.merit=std::max(r.merit,std::abs(v));
    auto& s=r.stages_v57.emplace();
    s.inverse_momentum_diagonal.assign(3,1.);
    s.effective_dynamic_viscosity=r.eddy_viscosity;
    s.strain_rate={std::sqrt(1.+u*u)};
    s.sst_blend_f2={0.6};
    s.wall_momentum_face_mu={};
    s.wall_production={-1.};
    s.momentum_convection={phi*u,0.,0.};
    return r;
}
int main(){
    PhysicalCandidateV51 x{1,{0.5,0.2,0.1,0.4,1.5,2.1}};
    const auto f0=evaluate(x);
    const std::vector<double> direction{0.2,0.1,0.03,0.04,0.12,0.1};
    const double h=1e-5;
    auto xp=x,xm=x;
    for(std::size_t i=0;i<6;++i){xp.fields[i]+=h*direction[i];xm.fields[i]-=h*direction[i];}
    const auto plus=evaluate(xp),minus=evaluate(xm);
    const std::array<double,6> ones{1.,1.,1.,1.,1.,1.};
    const auto conv=[](const PhysicalCandidateV51& x,const std::vector<double>& d,
        const std::vector<double>& phi,const std::vector<double>& dphi){
        return ConvectiveTangentV57{{dphi[0]*x.fields[0]+phi[0]*d[0],0.,0.},0};
    };
    const auto hybrid=combine_physical_jv_v57(x,f0,plus,minus,direction,
        ones,ones,h,conv,false,0);
    check(hybrid.action.size()==6 && hybrid.nonsmooth_faces==0 &&
          hybrid.convective_relative_discrepancy<1e-7,
          "V57 hybrid mismatched analytic/FD limited convection");
    const double phi=f0.face_flux[0];
    const double dphi=direction[0]+0.1*direction[4]+0.02*direction[3];
    const double truth=dphi*x.fields[0]+phi*direction[0]+
        0.1*std::exp(0.1*x.fields[3])*direction[3]-0.2*direction[5];
    check(std::abs(hybrid.action[0]-truth)<1e-8 &&
          std::abs(hybrid.action[3]-dphi)<1e-8 &&
          std::abs(hybrid.action[4]-(direction[4]+0.6*x.fields[0]*direction[0]))<1e-8,
          "V57 hybrid lost RC or momentum-SST cross term");
    const auto pin=combine_physical_jv_v57(x,f0,plus,minus,direction,
        ones,ones,h,conv,true,0);
    check(std::abs(pin.action[3]-direction[3])<1e-14,
          "V57 pressure gauge algebra improperly differentiated");
    bool rejected=false;
    try{
        const auto bad=[](const PhysicalCandidateV51&,
                const std::vector<double>&,const std::vector<double>&,
                const std::vector<double>&){return ConvectiveTangentV57{{0.,0.,0.},1};};
        auto unqualified=combine_physical_jv_v57(x,f0,plus,minus,
            direction,ones,ones,h,bad,false,0);
        check(unqualified.nonsmooth_faces>0 &&
              unqualified.convective_relative_discrepancy>0.1,
              "V57 nonsmooth/wrong analytic tangent was incorrectly certified");
    }catch(const std::runtime_error&){rejected=true;}
    check(!rejected,"V57 diagnostic unexpectedly threw on an explicitly bad tangent");
    std::cout<<"V57_HYBRID_PASS convective_gap="<<hybrid.convective_relative_discrepancy
             <<" dF="<<hybrid.action[0]<<"\n";
}
