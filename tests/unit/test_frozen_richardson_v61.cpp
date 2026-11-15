#include "cfdx/physics/global_physical_newton_v52.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
static void require(bool ok,const char* why){if(!ok) throw std::runtime_error(why);}

static PhysicalResidualV51 oracle(const PhysicalCandidateV51& x, bool nonsmooth=false){
    PhysicalResidualV51 r;
    const auto n=x.fields.size();
    r.signed_raw.assign(n,0.);r.normalized.assign(n,0.);
    r.eddy_viscosity.assign(x.cells,1.);
    r.face_flux.assign(x.cells+1,0.);
    for(std::size_t i=0;i<n;++i){
        const auto u=x.fields[i];
        r.signed_raw[i]=u*u*u+.5*u+.02*x.fields[(i+1)%n];
        r.normalized[i]=r.signed_raw[i];
    }
    if(nonsmooth){
        // A multiple branch at the base has no unique classical derivative.
        r.signed_raw[0]+=std::abs(x.fields[0]-.2);
        r.normalized[0]=r.signed_raw[0];
    }
    r.merit=0.;for(const auto v:r.normalized)r.merit=std::max(r.merit,std::abs(v));
    r.flux_closed=true;r.flux_fixed_point_defect=0.;
    return r;
}
int main(){
    PhysicalCandidateV51 x;
    x.cells=2;
    x.fields={.2,.4,.3,.5,.5,.6,.1,.2,.7,.9,1.4,2.2};
    const std::array<double,6> units{{1.,1.,1.,1.,1.,1.}};
    auto smooth=[&](const PhysicalCandidateV51& y){return oracle(y);};
    const auto base=smooth(x);
    FrozenPhysicalJacobianV54 old(x,base,units,units,false,0,0.,1e-3,smooth);
    FrozenPhysicalJacobianV54 fixed(x,base,units,units,false,0,0.,1e-3,
                                    smooth,0.,0.,true,0.05);
    require(old.statistics().evaluations==2*x.fields.size(),
            "V54 initial evaluation count changed by default");
    require(fixed.statistics().evaluations==4*x.fields.size() &&
            fixed.statistics().v61_columns_trusted &&
            fixed.statistics().v61_richardson_columns==x.fields.size() &&
            fixed.statistics().v61_rejected_columns==0,
            "V61 did not construct a certified smooth Richardson matrix");
    std::vector<double> e(x.fields.size(),0.);e[0]=1.;
    const double exact=3.*x.fields[0]*x.fields[0]+.5;
    const double e54=std::abs(old.apply(e)[0]-exact);
    const double e61=std::abs(fixed.apply(e)[0]-exact);
    require(e54>1e-7 && e61<e54/100. && e61<1e-8,
            "V61 did not remove the O(h^2) error of a smooth cubic");
    std::vector<double> u(e.size(),0.),v(e.size(),0.),w(e.size(),0.);
    for(std::size_t i=0;i<u.size();++i){
        const auto index=static_cast<double>(i);
        u[i]=std::sin(index+1.);v[i]=std::cos(.7*(index+1.));w[i]=u[i]-1.2*v[i];
    }
    const auto au=fixed.apply(u),av=fixed.apply(v),aw=fixed.apply(w);
    for(std::size_t i=0;i<w.size();++i)
        require(std::abs(aw[i]-au[i]+1.2*av[i])<1e-11,
                "V61 frozen Richardson action not additive");
    const auto directional=fixed.audit_direction(u,smooth);
    require(directional.fidelity<1e-4 && directional.halving<1e-3 &&
            directional.sided_nonsmoothness<0.05,
            "V61 smooth independent physical fidelity failed");
    auto kink=[&](const PhysicalCandidateV51& y){return oracle(y,true);};
    FrozenPhysicalJacobianV54 ambiguous(x,kink(x),units,units,false,0,0.,1e-4,
                                        kink,0.,0.,true,0.05);
    require(!ambiguous.statistics().v61_columns_trusted &&
            ambiguous.statistics().v61_rejected_columns>=1 &&
            ambiguous.statistics().v61_maximum_sided_disagreement>0.05,
            "V61 accepted a multiple active-set kink as a smooth Jacobian");
    auto settings=GlobalPhysicalNewtonControlsV52{};
    settings.v54_frozen_jacobian=true;
    settings.v61_richardson_frozen=true;
    settings.v57_true_krylov_each_step=true;
    settings.krylov_dimension=6;
    const auto rejected=global_physical_newton_v52(x,kink,settings);
    require(rejected.steps.size()==1 && !rejected.steps[0].accepted &&
            !rejected.steps[0].v61_column_gate_ok &&
            rejected.steps[0].v61_rejected_columns>=1 &&
            rejected.accepted_steps==0 && !rejected.production_committed,
            "V61 kink was promoted into an accepted physical Newton step");
    settings.v54_frozen_jacobian=false;
    bool threw=false;
    try{(void)global_physical_newton_v52(x,smooth,settings);}
    catch(const std::invalid_argument&){threw=true;}
    require(threw,"V61 optional mode must require frozen Jacobian");
    std::cout<<"V61_RICHARDSON_PASS old_error="<<e54
             <<" new_error="<<e61
             <<" ambiguous_columns="<<ambiguous.statistics().v61_rejected_columns
             <<" newton_committed=0\n";
}
