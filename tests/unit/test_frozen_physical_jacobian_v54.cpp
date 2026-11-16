#include "cfdx/physics/frozen_physical_jacobian_v54.h"
#include "cfdx/physics/global_physical_newton_v52.h"
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}

// A genuinely coupled, dense, candidate-owned six-field function on 2 cells.
// Its exact Jacobian has nonzero momentum <-> SST off-diagonal derivatives.
static PhysicalResidualV51 linear_physical(const PhysicalCandidateV51& x){
    PhysicalResidualV51 r;
    const auto n=x.fields.size();
    r.signed_raw.assign(n,0.);
    r.merit=0.;
    r.normalized.assign(n,0.);
    r.face_flux={x.fields[0]+x.fields[3],x.fields[1]-x.fields[3]};
    r.eddy_viscosity.assign(x.cells,1.);
    for(std::size_t i=0;i<n;++i){
        r.signed_raw[i]=2.*x.fields[i]-0.3;
        for(std::size_t j=0;j<n;++j)
            if(i!=j)r.signed_raw[i]+=0.005*(1.0+static_cast<double>((i+j)%4))*x.fields[j];
        r.normalized[i]=r.signed_raw[i];
        r.merit=std::max(r.merit,std::abs(r.normalized[i]));
    }
    r.flux_closed=true;r.flux_fixed_point_defect=0.;
    return r;
}
int main(){
    PhysicalCandidateV51 x; x.cells=2;
    x.fields={0.2,-0.1,0.1,0.2,0.,0.,0.1,0.2,1.2,1.4,2.2,2.5};
    const auto base=linear_physical(x);
    const std::array<double,6> state{{1.,1.,1.,1.,2.,3.}};
    const std::array<double,6> rows{{1.,1.,1.,1.,1.,1.}};
    FrozenPhysicalJacobianV54 a(x,base,state,rows,true,0,0.,1e-5,
                                 linear_physical);
    require(a.size()==12 && a.statistics().rows==12 &&
            a.statistics().stored_bytes==12*12*sizeof(double),
            "V54 explicit matrix dimension or storage wrong");
    require(a.statistics().evaluations>=12 &&
            a.statistics().centered_columns+a.statistics().one_sided_columns==12,
            "V54 FD columns not checked");
    std::vector<double> u(12),v(12),w(12);
    for(std::size_t i=0;i<12;++i){
        u[i]=std::sin(static_cast<double>(i+1));
        v[i]=std::cos(0.3*static_cast<double>(i+1));
        w[i]=0.7*u[i]-1.3*v[i];
    }
    const auto au=a.apply(u),av=a.apply(v),aw=a.apply(w);
    FrozenPhysicalJacobianV54 repeat(x,base,state,rows,true,0,0.,1e-5,
                                    linear_physical);
    require(au==repeat.apply(u),"V54 same-state Jacobian replay was not deterministic");
    for(std::size_t i=0;i<12;++i)
        require(std::abs(aw[i]-0.7*au[i]+1.3*av[i])<1e-11,
                "V54 frozen Amat is not additive on Krylov directions");
    const auto evidence=a.audit_direction(u,linear_physical);
    require(evidence.fidelity<1e-8 && evidence.halving<1e-8 &&
            evidence.sided_nonsmoothness<1e-8,
            "V54 exact affine Jacobian disagrees with physical F");
    // A general dense Jacobian should work with the SAME frozen Arnoldi action
    // for the true residual (not a newly chosen finite difference step).
    const auto identity=[](const std::vector<double>& r){return r;};
    std::vector<double> rhs(12,1.),solution;
    double true_residual=0.,arnoldi=0.,nonadditive=0.;
    const bool accepted=global_physical_fgmres_v52(
        [&](const std::vector<double>& d){return a.apply(d);},
        identity,rhs,12,1e-10,solution,true_residual,&arnoldi,&nonadditive);
    require(accepted && true_residual<1e-10 &&
            nonadditive<1e-10 && arnoldi<1e-10,
            "V54 Krylov true residual/Arnoldi disagreement");
    auto invalid=u;invalid[3]=std::numeric_limits<double>::quiet_NaN();
    bool threw=false;
    try{(void)a.apply(invalid);}catch(const std::invalid_argument&){threw=true;}
    require(threw,"V54 nonfinite input silently accepted");
    threw=false;
    try{(void)a.apply({1.,2.});}catch(const std::invalid_argument&){threw=true;}
    require(threw,"V54 incomplete Krylov input silently accepted");
    // Test the full V52 loop with the V54 operator AND independent audits.
    auto controls=GlobalPhysicalNewtonControlsV52{};
    controls.v54_frozen_jacobian=true;
    controls.maximum_steps=1;controls.krylov_dimension=12;
    controls.linear_relative_tolerance=1e-5;
    auto result=global_physical_newton_v52(x,linear_physical,controls);
    require(result.steps.size()==1 && result.steps[0].v54_frozen_jacobian_ready &&
            result.steps[0].v54_basis_audits>=1 &&
            result.steps[0].v54_krylov_fidelity_ok &&
            result.steps[0].linear_true_relative<controls.linear_relative_tolerance &&
            !result.production_committed,"V54 frozen six-field Newton contract failed");
    // A deliberately non-smooth function must never be promoted on the
    // strength of a pretty Arnoldi contraction alone.
    auto nonsmooth=[](const PhysicalCandidateV51& candidate){
        auto r=linear_physical(candidate);
        r.signed_raw[0]+=std::abs(candidate.fields[0]-0.2);
        r.normalized[0]=r.signed_raw[0];
        r.merit=0.;for(double row:r.normalized)r.merit=std::max(r.merit,std::abs(row));
        return r;
    };
    FrozenPhysicalJacobianV54 kink(x,nonsmooth(x),state,rows,false,0,0.,1e-5,nonsmooth);
    const auto directional=kink.audit_direction(u,nonsmooth);
    require(std::isfinite(directional.fidelity) &&
            std::isfinite(directional.halving) &&
            directional.sided_nonsmoothness>0.05,
            "V54 kink probe should remain a finite, nonqualified diagnostic");
    auto zero=u;for(double& c:zero)c=0.;
    threw=false;
    try{(void)kink.audit_direction(zero,nonsmooth);}
    catch(const std::invalid_argument&){threw=true;}
    require(threw,"V54 zero-direction fidelity must fail closed");
}
