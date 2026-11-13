#pragma once
// V55 experimental *sparse/matrix-free* correction, not a dense reference:
// use V49's AMG as a preconditioner, but apply the full V54 WLS stencil in
// portable two-double arithmetic. Accepted updates decrease the reassembled
// two-double defect. The original V48 DOUBLE defect is checked separately.
// No state is ever published to CFDX's production solver or to HO-LO.
#include "cfdx/physics/dom_wls_hilo_gate_v55.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct DomV55RefinementOptions {
    std::size_t inner_budget=240;
    int restart=80;
    unsigned max_updates=8;
    double inner_tolerance=1e-5;
    double strict=1e-13;
};
struct DomV55Refinement {
    bool amg_setup=false,split_strict=false,original_double_strict=false;
    bool publish=false; // independent angular/CHT gates NOT supplied
    double initial_split_l2=std::numeric_limits<double>::infinity();
    double final_split_l2=std::numeric_limits<double>::infinity();
    double final_original_double_l2=std::numeric_limits<double>::infinity();
    std::size_t krylov_iterations=0,accepted_updates=0,rejected_updates=0;
    std::vector<DomHiLoV55> diagnostic_state;
};
inline DomV55Refinement solve_wls_hilo_refinement_v55(
        const DomStableWlsV54& stable,
        const DomSmmGrayWallV49& original,
        const std::vector<double>& rhs,
        const std::vector<double>& initial,
        const DomV55RefinementOptions& cfg={}) {
    using namespace cfdx::core;
    const std::size_t n=rhs.size();
    if(n<8 || original.dimension()!=n || initial.size()!=n ||
       !(cfg.strict>0.) || !(cfg.inner_tolerance>0.) ||
       cfg.inner_tolerance>=1. || cfg.restart<2 ||
       cfg.inner_budget<static_cast<std::size_t>(cfg.restart) ||
       cfg.max_updates==0 || cfg.max_updates>16)
        throw std::invalid_argument("V55 WLS refine: invalid controls");
    long double den=0.;
    for(std::size_t i=0;i<n;++i){
        if(!std::isfinite(rhs[i]) || !std::isfinite(initial[i]) ||
           initial[i]<0.)
            throw std::invalid_argument("V55 WLS refine: nonphysical inputs");
        den+=static_cast<long double>(rhs[i])*rhs[i];
    }
    if(!(den>0.))throw std::invalid_argument("V55 WLS refine: zero RHS");
    LinearOperator op;
    op.size=n;
    op.apply=[&stable,n](const Vector& v,Vector& y){
        std::vector<DomHiLoV55> direction(n);
        for(std::size_t i=0;i<n;++i)direction[i]=DomHiLoV55(v(i));
        const auto result=stable.apply_hilo_v55(direction);
        if(y.size()!=n)y.resize(n);
        for(std::size_t i=0;i<n;++i){
            y(i)=result[i].hi+result[i].lo;
            if(!std::isfinite(y(i)))
                throw std::runtime_error("V55 WLS refine: nonfinite matvec");
        }
    };
    FunctionalLinearOperator fine(n,op.apply);
    MatrixFreeVcyclePreconditioner amg(fine,.55,2,2,.25,12,
        AMGInterpolationPolicy::DirectionalSemiCoarsening,false,false,
        AMGSmootherPolicy::AlgebraicLineGaussSeidel);
    DomV55Refinement out;
    if(!amg.setup(original.auxiliary_matrix()))return out;
    out.amg_setup=true;
    std::vector<DomHiLoV55> state(n);
    for(std::size_t i=0;i<n;++i)state[i]=DomHiLoV55(initial[i]);
    const auto residual=[&](const std::vector<DomHiLoV55>& x,
                            std::vector<DomHiLoV55>* r){
        const auto y=stable.apply_hilo_v55(x);
        long double acc=0.;
        if(r)r->resize(n);
        for(std::size_t i=0;i<n;++i){
            const DomHiLoV55 defect=DomHiLoV55(rhs[i])-y[i];
            if(!dom_finite_v55(defect))
                throw std::runtime_error("V55 WLS refine: nonfinite true defect");
            const long double d=dom_to_long_double_v55(defect);
            acc+=d*d;
            if(r)(*r)[i]=defect;
        }
        return static_cast<double>(std::sqrt(acc/den));
    };
    out.initial_split_l2=residual(state,nullptr);
    out.final_split_l2=out.initial_split_l2;
    for(unsigned update=0;update<cfg.max_updates &&
                          out.final_split_l2>cfg.strict;++update){
        std::vector<DomHiLoV55> true_r;
        const double before=residual(state,&true_r);
        Vector b(n,0.),dx(n,0.);
        for(std::size_t i=0;i<n;++i)b(i)=true_r[i].hi+true_r[i].lo;
        const auto info=solve_gmres(op,b,dx,cfg.restart,cfg.inner_budget,
                                    cfg.inner_tolerance,&amg);
        out.krylov_iterations+=info.iterations;
        if(info.status!=SolverStatus::CONVERGED){
            ++out.rejected_updates;break;
        }
        bool improved=false;
        for(double alpha=1.;alpha>=1./64.;alpha*=.5){
            auto candidate=state;
            bool physical=true;
            for(std::size_t i=0;i<n;++i){
                candidate[i]+=DomHiLoV55(alpha)*DomHiLoV55(dx(i));
                if(!dom_finite_v55(candidate[i]) ||
                   candidate[i].hi<0. ||
                   (candidate[i].hi==0. && candidate[i].lo<0.)){
                    physical=false;break;
                }
            }
            if(!physical)continue;
            const double after=residual(candidate,nullptr);
            if(std::isfinite(after) && after<before){
                out.final_split_l2=after;
                state.swap(candidate);
                ++out.accepted_updates;
                improved=true;
                break;
            }
        }
        if(!improved){++out.rejected_updates;break;}
    }
    out.final_split_l2=residual(state,nullptr);
    out.split_strict=out.final_split_l2<=cfg.strict;
    std::vector<double> rounded(n);
    for(std::size_t i=0;i<n;++i)rounded[i]=state[i].hi;
    const auto orig=original.apply(rounded);
    long double num=0.;
    for(std::size_t i=0;i<n;++i){
        const long double d=static_cast<long double>(rhs[i])-orig[i];
        num+=d*d;
    }
    out.final_original_double_l2=static_cast<double>(std::sqrt(num/den));
    out.original_double_strict=out.final_original_double_l2<=cfg.strict;
    out.diagnostic_state=std::move(state);
    out.publish=false;
    return out;
}
} // namespace cfdx::physics
