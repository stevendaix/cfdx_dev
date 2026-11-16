#include "cfdx/physics/global_sst_block_v53.h"
#include "cfdx/physics/jv_fidelity_v53.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool value,const char* text){if(!value)throw std::runtime_error(text);}

static PhysicalResidualV51 evaluate(const PhysicalCandidateV51& x) {
    const auto n=x.cells;
    PhysicalResidualV51 r;
    r.signed_raw.assign(6*n,0.0);
    r.normalized.assign(6*n,0.0);
    r.eddy_viscosity.assign(n,0.0);
    r.face_flux.assign(n,0.0);
    for(std::size_t i=0;i<n;++i) {
        const auto at=[&](std::size_t b,std::size_t j){return x.fields[b*n+j];};
        r.signed_raw[i]=at(0,i)+0.04*at(4,i)-1.;
        r.signed_raw[n+i]=at(1,i)+0.06*at(5,i)+0.3;
        r.signed_raw[2*n+i]=at(2,i)-0.1;
        r.signed_raw[3*n+i]=at(3,i)+0.03*at(4,i)-0.2;
        const double left_k=i?at(4,i-1):0.0,right_k=i+1<n?at(4,i+1):0.0;
        const double left_w=i?at(5,i-1):0.0,right_w=i+1<n?at(5,i+1):0.0;
        r.signed_raw[4*n+i]=3.*at(4,i)-0.2*left_k-0.2*right_k+
             0.1*at(5,i)+0.08*at(0,i)+0.04*at(3,i)-2.;
        r.signed_raw[5*n+i]=2.*at(5,i)-0.1*left_w-0.1*right_w+
             0.15*at(4,i)+0.05*at(1,i)+0.06*at(3,i)-1.;
        r.eddy_viscosity[i]=at(4,i)/at(5,i);
        r.face_flux[i]=at(0,i)+at(1,i);
    }
    r.normalized=r.signed_raw;r.merit=0.;
    for(double a:r.normalized)r.merit=std::max(r.merit,std::abs(a));
    r.flux_closed=true;r.flux_fixed_point_defect=0.;
    return r;
}
int main() {
    constexpr std::size_t n=9;
    PhysicalCandidateV51 base;
    base.cells=n;base.fields.resize(6*n);
    for(std::size_t b=0;b<6;++b)
        for(std::size_t i=0;i<n;++i)
            base.fields[b*n+i]=(b==5?1.8:b==4?0.7:0.1);
    cfdx::core::SparseMatrix k(n,n),w(n,n);
    std::vector<std::vector<std::size_t>> graph(n);
    for(std::size_t i=0;i<n;++i){
        k.push_back(i,i,3.);w.push_back(i,i,2.);
        if(i){k.push_back(i,i-1,-0.2);w.push_back(i,i-1,-0.1);graph[i].push_back(i-1);}
        if(i+1<n){k.push_back(i,i+1,-0.2);w.push_back(i,i+1,-0.1);graph[i].push_back(i+1);}
    }
    k.finalize();w.finalize();
    const GlobalSSTBlockV53::Scales scales{{1,1,1,1,1,1}};
    const auto native=[](const std::vector<double>& v,
                         const GlobalSSTBlockV53::Scales&,
                         const GlobalSSTBlockV53::Scales&,std::size_t){return v;};
    GlobalSSTBlockV53 m(base,evaluate(base),scales,scales,k,w,graph,
                       evaluate,native,false,0);
    require(m.statistics().colors>1 && m.statistics().candidate_evaluations>6,
            "V53 did not probe a colored local Jacobian");
    require(m.statistics().transport_nnz>2*n,"V53 lost FV transport neighbors");
    require(m.statistics().maximum_locality_relative_error<1e-7,
            "V53 grouping is not equivalent to independent columns");
    std::vector<double> rhs(6*n,0.);
    for(std::size_t i=0;i<n;++i){
        rhs[i]=0.15;rhs[n+i]=-0.07;rhs[2*n+i]=0.03;rhs[3*n+i]=0.12;
        rhs[4*n+i]=1.2;rhs[5*n+i]=0.7;
    }
    const auto z=m.apply(rhs);
    // The test function is linear: F(base+z)-F(base) is exactly J*z.
    auto shifted=base;
    for(std::size_t i=0;i<6*n;++i)shifted.fields[i]+=z[i];
    const auto fp=evaluate(shifted),f0=evaluate(base);
    double error=0.,rnorm=0.;
    for(std::size_t i=0;i<6*n;++i) {
        const double result=fp.signed_raw[i]-f0.signed_raw[i];
        error+=std::pow(result-rhs[i],2);
        rnorm+=std::pow(rhs[i],2);
    }
    require(std::sqrt(error/rnorm)<0.002,"V53 coupled block inverse too inaccurate");
    const auto jv=audit_physical_jv_v53(base,rhs,evaluate,1e-5);
    require(jv.both_sides_valid && jv.central_halving_relative<1e-8 &&
            jv.sided_nonsmooth_relative<1e-8,
            "V53 smooth directional derivative is inconsistent");
    // Failure must never cause implicit identity fallback.
    bool rejected=false;
    try{
        auto invalid=graph;invalid[0].push_back(n+1);
        GlobalSSTBlockV53 bad(base,evaluate(base),scales,scales,k,w,invalid,
                              evaluate,native,false,0);
        (void)bad;
    }catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"V53 invalid mesh graph admitted");
    rejected=false;
    try{
        GlobalSSTBlockV53 bad(base,evaluate(base),scales,scales,k,w,graph,
                              evaluate,{},false,0);
        (void)bad;
    }catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"V53 missing 4N Schur admitted");
    // Non-smooth |x| branch is deliberately detected rather than called
    // an accurate classical two-sided Jacobian.
    auto nonsmooth=[](const PhysicalCandidateV51& x){
        auto r=evaluate(x);
        r.signed_raw[0]=std::abs(x.fields[0]-0.1);
        r.normalized[0]=r.signed_raw[0];
        r.merit=0.;for(double a:r.normalized)r.merit=std::max(r.merit,std::abs(a));
        return r;
    };
    std::vector<double> dir(6*n,0.);dir[0]=1.;
    const auto corner=audit_physical_jv_v53(base,dir,nonsmooth,1e-5);
    require(corner.sided_nonsmooth_relative>0.5,
            "V53 failed to identify a nonsmooth branch");
}
