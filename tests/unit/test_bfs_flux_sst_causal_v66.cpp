@@ -0,0 +1,48 @@
#include "cfdx/physics/bfs_flux_sst_causal_v66.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static void need(bool c,const char* why){if(!c)throw std::runtime_error(why);}
template<class F>static void rejects(F&& f){bool yes=false;try{f();}catch(const std::exception&){yes=true;} need(yes,"V66 negative test did not reject");}
static PhysicalCandidateV51 state(){PhysicalCandidateV51 s;s.cells=2;s.fields={3,4,5,6,0,0,1,1,0.1,0.2,3,3};return s;}
static PhysicalResidualV51 result(){
    PhysicalResidualV51 r;r.signed_raw={.1,0,0,0,0,0, .04,-.02,.02,-.02,.05,.03};
    r.normalized=r.signed_raw;r.face_flux={.02,.02,0};r.eddy_viscosity={.001,.002};
    r.inverse_momentum_diagonal={.1,.1,.2,.2,.3,.3};
    r.flux_closed=true;r.flux_fixed_point_defect=1e-13;r.merit=1.;return r;
}
int main(){
    const auto s=state();const auto r=result();
    auto t=s;auto tr=r;t.fields[5*2+1]+=1.;tr.signed_raw[5*2+1]=4.;tr.normalized=tr.signed_raw;
    tr.eddy_viscosity[1]=.003;tr.inverse_momentum_diagonal[1]=.12;
    auto owner=[](std::size_t f){return f==2?1u:0u;};
    auto neighbour=[](std::size_t f){return f==1?1:-1;};
    auto x=audit_bfs_flux_sst_v66(s,t,r,tr,r.face_flux,{.025,.018,0},owner,neighbour);
    need(x.valid && x.omega_worst_change.cell==1 && x.omega_worst_trial.cell==1,
         "omega hotspot must follow signed physical row");
    need(std::abs(x.omega_worst_change.delta_signed-3.97)<1e-13,
         "omega signed residual delta");
    need(x.trial_carrier.worst_face==0 && x.trial_carrier.worst_cell==0,
         "face and cell conservative mismatch");
    need(x.rau_worst_cell==1 && x.nut_worst_cell==1,
         "rAU / nut hotspots");
    rejects([&]{auto bad=tr;bad.inverse_momentum_diagonal.clear();audit_bfs_flux_sst_v66(s,t,r,bad,r.face_flux,r.face_flux,owner,neighbour);});
    rejects([&]{auto bad=r.face_flux;bad[0]=INFINITY;v66_carrier_gap(s,r,bad,owner,neighbour);});
    rejects([&]{auto bad=r.face_flux;bad.pop_back();v66_carrier_gap(s,r,bad,owner,neighbour);});
    rejects([&]{auto bad=r;bad.flux_closed=false;audit_bfs_flux_sst_v66(s,t,bad,tr,r.face_flux,r.face_flux,owner,neighbour);});
    const auto gap=v66_coefficient_gap({.1,.2}, {.1,.4},2);
    need(gap.worst_cell==1 && std::abs(gap.relative_linf-.5)<1e-15,
         "native-vs-carrier coefficient comparison");
    rejects([&]{v66_coefficient_gap({.1,.2},{.1,-.4},2);});
    V66OmegaBudgetParts old{1.,2.,-3.,4.,-.5,3.5,0.,0.,0.,false};
    auto next=old;next.convection=1.5;next.physical_production=-2.;next.full_residual=5.;
    auto d=v66_compare_omega_budgets(old,next);
    need(std::abs(d.sum-d.residual)<1e-14 && std::abs(d.production-1.)<1e-14,
         "omega budget signed conservation");
    rejects([&]{auto bad=next;bad.full_residual+=.1;v66_compare_omega_budgets(old,bad);});
    rejects([&]{auto bad=next;bad.constraint=true;v66_compare_omega_budgets(old,bad);});
    rejects([&]{auto bad=next;bad.source_identity_error=1.;v66_compare_omega_budgets(old,bad);});
    std::cout<<"V66_BFS_FLUX_SST_CAUSAL_PASS\n";
}
