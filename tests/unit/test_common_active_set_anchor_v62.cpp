#include "cfdx/physics/global_physical_newton_v52.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void must(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
static PhysicalResidualV51 oracle(const PhysicalCandidateV51& x){
    PhysicalResidualV51 r;
    const std::size_t n=x.fields.size();
    r.signed_raw.resize(n);r.normalized.resize(n);
    r.eddy_viscosity.assign(x.cells,1.0);
    r.face_flux.assign(x.cells+1,0.0);
    for(std::size_t i=0;i<n;++i){
        const double u=x.fields[i];
        const double v=x.fields[(i+1)%n];
        r.signed_raw[i]=u*u*u+0.5*u+0.02*v;
    }
    // A multiple-branch limiter surrogate: min on both equal base states.
    // A common nearby differentiable state has a single B-subdifferential.
    r.signed_raw[0]+=std::min(x.fields[0],x.fields[1]);
    r.normalized=r.signed_raw;
    r.merit=0.;
    for(double value:r.normalized)r.merit=std::max(r.merit,std::abs(value));
    r.flux_closed=true;r.flux_fixed_point_defect=0.;
    return r;
}
// Controlled cross-equation surrogate. Derivatives are chosen so that a
// momentum-SST coupling dropped from the 6N Jacobian is detected, including
// a pressure/gauge row; this is NOT an SST physical validation case.
static PhysicalResidualV51 coupled_sst_oracle(const PhysicalCandidateV51& x){
    auto r=oracle(x);
    const auto n=x.cells;
    const double ux=x.fields[0],p=x.fields[3*n],k=x.fields[4*n];
    const double omega=x.fields[5*n];
    r.signed_raw[0]+=0.3*ux*k/omega+0.07*p;
    r.signed_raw[4*n]+=0.1*ux*omega;
    r.signed_raw[5*n]+=0.2*ux*k;
    r.normalized=r.signed_raw;
    r.merit=0.;
    for(double v:r.normalized)r.merit=std::max(r.merit,std::abs(v));
    r.face_flux[0]=ux+0.01*p+0.001*k;
    return r;
}
int main(){
    PhysicalCandidateV51 x;
    x.cells=2;
    x.fields={.2,.2,.3,.5,0.,0.,.1,.2,.7,.9,1.4,2.2};
    const std::array<double,6> scales{{1.,1.,1.,1.,1.,1.}};
    const auto base=oracle(x);
    const auto selector=generalized_anchor_candidate_v62(x,scales,true,0,1e-3);
    must(selector.fields[2*2]==0. && selector.fields[2*2+1]==0.,
         "V62 introduced fictitious Uz into the extruded 2D BFS");
    must(selector.fields[3*2]==x.fields[3*2],"V62 shifted a pinned pressure gauge");
    must(selector.fields[0]!=selector.fields[1],"V62 did not split active-set tie");
    must(selector.fields[4*2]>=x.fields[4*2] &&
         selector.fields[5*2]>x.fields[5*2],"V62 left SST admissible region");
    FrozenPhysicalJacobianV54 center(x,base,scales,scales,true,0,0.1,1e-5,
                                     oracle,0.,0.,true,0.05);
    FrozenPhysicalJacobianV54 anchored(x,base,scales,scales,true,0,0.1,1e-5,
                                       oracle,0.,0.,true,0.05,true);
    FrozenPhysicalJacobianV54 replay(x,base,scales,scales,true,0,0.1,1e-5,
                                     oracle,0.,0.,true,0.05,true);
    must(!center.statistics().v61_columns_trusted,
         "V61 unexpectedly classified exact multiple-bound tie as smooth");
    must(anchored.statistics().v62_anchor_requested &&
         anchored.statistics().v62_anchor_within_guard &&
         anchored.statistics().v61_columns_trusted &&
         anchored.statistics().v61_rejected_columns==0,
         "V62 did not obtain a consistent common smooth active-set branch");
    must(anchored.anchor_candidate_v62().fields==replay.anchor_candidate_v62().fields,
         "V62 common anchor is not reproducible");
    const auto anchor=anchored.anchor_candidate_v62();
    std::vector<double> e0(x.fields.size(),0.),e1(x.fields.size(),0.);
    e0[0]=1.;e1[1]=1.;
    const auto d0=anchored.apply(e0);
    const auto d1=anchored.apply(e1);
    const double known0=3.*anchor.fields[0]*anchor.fields[0]+.5;
    const double tie_gradient=(d0[0]-known0)+(d1[0]-.02);
    must(std::abs(tie_gradient-1.)<1e-7,
         "V62 selected inconsistent column-wise active sets for min limiter");
    if(anchor.fields[0]<anchor.fields[1]){
        must(std::abs(d0[0]-known0-1.)<1e-7 &&
             std::abs(d1[0]-.02)<1e-7,"V62 selected wrong common min branch");
    }else{
        must(std::abs(d0[0]-known0)<1e-7 &&
             std::abs(d1[0]-.02-1.)<1e-7,"V62 selected wrong common min branch");
    }
    std::vector<double> u(x.fields.size(),0.),v(x.fields.size(),0.);
    for(std::size_t i=0;i<u.size();++i){
        u[i]=std::sin(static_cast<double>(i)+.6);
        v[i]=std::cos(static_cast<double>(i)*.6+.3);
    }
    std::vector<double> sum=u;
    for(std::size_t i=0;i<sum.size();++i)sum[i]+=v[i];
    const auto au=anchored.apply(u),av=anchored.apply(v),as=anchored.apply(sum);
    for(std::size_t i=0;i<sum.size();++i)
        must(std::abs(as[i]-au[i]-av[i])<1e-11,
             "V62 anchored Amat is not additive under Arnoldi");
    const auto physical=coupled_sst_oracle(x);
    FrozenPhysicalJacobianV54 cross(x,physical,scales,scales,true,0,0.1,1e-5,
                                   coupled_sst_oracle,0.,0.,true,0.05,true);
    must(cross.statistics().v61_columns_trusted,"V62 lost a smooth momentum-SST cross derivative");
    const auto at=cross.anchor_candidate_v62();
    e0.assign(x.fields.size(),0.);e1.assign(x.fields.size(),0.);
    e0[4*x.cells]=1.;e1[5*x.cells]=1.;
    const auto d_k=cross.apply(e0),d_omega=cross.apply(e1);
    const double u0=at.fields[0],k0=at.fields[4*x.cells],w0=at.fields[5*x.cells];
    must(std::abs(d_k[0]-.3*u0/w0)<1e-7,
         "V62 silently dropped d(momentum)/dk from the full 6N reassembly");
    must(std::abs(d_omega[0]+.3*u0*k0/(w0*w0))<1e-7,
         "V62 silently dropped d(momentum)/domega from the full 6N reassembly");
    must(std::abs(d_k[5*x.cells]-.2*u0)<1e-7,
         "V62 silently dropped d(omega)/dk");
    must(std::abs(d_omega[4*x.cells]-.1*u0)<1e-7,
         "V62 silently dropped d(k)/domega");
    e0.assign(x.fields.size(),0.);e0[3*x.cells]=1.;
    const auto dp=cross.apply(e0);
    must(std::abs(dp[0]-.07)<1e-7 &&
         std::abs(dp[3*x.cells]-1.)<1e-7,
         "V62 lost physical pressure coupling or changed pressure gauge");
    GlobalPhysicalNewtonControlsV52 controls;
    controls.v54_frozen_jacobian=true;
    controls.v61_richardson_frozen=true;
    controls.v62_common_active_set_anchor=true;
    controls.v57_true_krylov_each_step=true;
    controls.pin_pressure=true;
    controls.reference_pressure=0.1;
    controls.reference_cell=0;
    controls.difference_step=1e-5;
    controls.krylov_dimension=3;
    const auto result=global_physical_newton_v52(x,oracle,controls);
    must(!result.production_committed && !result.physically_qualified,
         "V62 diagnostic prematurely published a nonlinear production state");
    controls.v61_richardson_frozen=false;
    bool rejected=false;
    try{(void)global_physical_newton_v52(x,oracle,controls);}
    catch(const std::invalid_argument&){rejected=true;}
    must(rejected,"V62 bypassed frozen Richardson/physical truth prerequisites");
    std::cout<<"V62_COMMON_ANCHOR_PASS tie_gradient="<<tie_gradient
             <<" anchor_distance="<<anchored.statistics().v62_anchor_normalized_distance
             <<" anchor_physical_shift="<<anchored.statistics().v62_anchor_physical_displacement
             <<" newton_acceptances="<<result.accepted_steps<<" qualified=0\n";
}
