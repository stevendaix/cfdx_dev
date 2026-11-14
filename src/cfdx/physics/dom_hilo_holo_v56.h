#pragma once
// CFDX V56 (research only). Carry a two-double cell state through the
// refreshed, globally accepted HO-SMM/AMG iteration. Angular HO is still
// double-based: independent V54 angular certification is mandatory and the
// invisible low component is bounded before HO projection. No production API
// or historic double convergence contract is silently replaced.
#include "cfdx/physics/dom_holo_certificate_v54.h"
#include "cfdx/physics/dom_wls_iterative_refinement_v55.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct DomV56FaceClosure {
    std::vector<DomHiLoV55> corrected_rhs, ho_minus_lo_face;
    double identity=std::numeric_limits<double>::infinity();
    double physical_ho=std::numeric_limits<double>::infinity();
};

// Independent HO currents enter *once* per face. The same two-double face
// correction is added to the owner and subtracted from the neighbour. Unlike
// a double LO correction, its cancellation cannot discard the low part of J.
// The resulting identity is an assembly check, NOT the angular certificate.
inline DomV56FaceClosure build_dom_hilo_face_closure_v56(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const DomSmmGrayWallV49& lo,const DomStableWlsV54& stable,
    const std::vector<DomHiLoV55>& projected,
    const std::vector<double>& gas_ib,
    const std::vector<double>& absorption,
    const std::vector<double>& ho_face_current) {
    using H=DomHiLoV55;
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(projected.size()!=nc || gas_ib.size()!=nc || absorption.size()!=nc ||
       ho_face_current.size()!=nf || lo.dimension()!=nc)
        throw std::invalid_argument("V56 closure: incompatible dimensions");
    const auto L=stable.evaluate_hilo_v55(projected);
    DomV56FaceClosure report;
    report.ho_minus_lo_face.resize(nf);
    report.corrected_rhs.resize(nc);
    std::vector<H> low_wall(nf),corrected_div(nc),ho_div(nc),base(nc);
    for(const auto& w:lo.walls()) {
        low_wall[w.face]=H(w.conductance)*(projected[w.owner]-H(w.ib));
        base[w.owner]+=H(w.conductance)*H(w.ib);
    }
    const H fourpi=dom_four_pi_v55();
    for(std::size_t c=0;c<nc;++c) {
        if(!dom_finite_v55(projected[c]) || !std::isfinite(gas_ib[c]) ||
           !std::isfinite(absorption[c]) || absorption[c]<=0. ||
           !(geometry.cell_volumes[c]>0.))
            throw std::invalid_argument("V56 closure: invalid optical source");
        base[c]+=fourpi*H(absorption[c])*H(geometry.cell_volumes[c])*H(gas_ib[c]);
    }
    for(std::size_t f=0;f<nf;++f) {
        if(!std::isfinite(ho_face_current[f]))
            throw std::invalid_argument("V56 closure: nonfinite angular flux");
        const auto owner=mesh.ownership().owner(f);
        const auto nb=mesh.ownership().neighbour(f);
        const H low=nb>=0 ? L.face_watts[f] : low_wall[f];
        const H high(ho_face_current[f]);
        const H delta=high-low;
        report.ho_minus_lo_face[f]=delta;
        corrected_div[owner]+=delta;
        ho_div[owner]+=high;
        if(nb>=0) {
            const auto j=static_cast<std::size_t>(nb);
            corrected_div[j]-=delta;
            ho_div[j]-=high;
        }
    }
    report.identity=0.;
    report.physical_ho=0.;
    // Reassemble HO energy directly from its own signed face currents and
    // true absorption. Never use HO-LO identity as a substitute for V54.
    for(std::size_t c=0;c<nc;++c) {
        report.corrected_rhs[c]=base[c]-corrected_div[c];
        const H algebraic=L.divergence_watts[c]-report.corrected_rhs[c];
        const H physical=ho_div[c]+fourpi*H(absorption[c])*
            H(geometry.cell_volumes[c])*(projected[c]-H(gas_ib[c]));
        const long double ar=dom_to_long_double_v55(algebraic);
        const long double pr=dom_to_long_double_v55(physical);
        const long double scale=std::max({1.L,std::abs(ar),std::abs(pr)});
        report.identity=std::max(report.identity,
            static_cast<double>(std::abs(ar-pr)/scale));
        report.physical_ho=std::max(report.physical_ho,
            static_cast<double>(std::abs(pr)/scale));
    }
    return report;
}

struct DomV56Linear {
    bool amg_setup=false,strict=false;
    double initial=std::numeric_limits<double>::infinity();
    double final=std::numeric_limits<double>::infinity();
    std::size_t krylov=0,accepted_updates=0,rejected_updates=0;
    std::vector<DomHiLoV55> split;
};

// V55 algebraic refinement generalized to a true hi+lo RHS and hi+lo start.
// Arnoldi is only an inner direction mechanism: every accepted step is gated
// on a reassembled, non-preconditioned hi+lo residual.
inline DomV56Linear solve_dom_hilo_smm_amg_v56(
    const DomSmmGrayWallV49& lo,const DomStableWlsV54& stable,
    const std::vector<DomHiLoV55>& rhs,
    const std::vector<DomHiLoV55>& initial,
    const DomV55RefinementOptions& cfg={}) {
    using namespace cfdx::core;
    const std::size_t n=lo.dimension();
    if(n<8 || rhs.size()!=n || initial.size()!=n ||
       !(cfg.strict>0.) || cfg.strict>1e-13 ||
       !(cfg.inner_tolerance>0.) || cfg.inner_tolerance>=1. ||
       cfg.restart<2 || cfg.inner_budget<static_cast<std::size_t>(cfg.restart) ||
       cfg.max_updates==0 || cfg.max_updates>16)
        throw std::invalid_argument("V56: invalid / loosened linear controls");
    long double den=0.;
    for(std::size_t c=0;c<n;++c) {
        if(!dom_finite_v55(rhs[c]) || !dom_finite_v55(initial[c]) ||
           initial[c].hi<0. || (initial[c].hi==0. && initial[c].lo<0.))
            throw std::invalid_argument("V56: nonphysical split inputs");
        const long double v=dom_to_long_double_v55(rhs[c]);
        den+=v*v;
    }
    if(!(den>0.))throw std::invalid_argument("V56: zero RHS norm");
    LinearOperator op;
    op.size=n;
    op.apply=[&stable,n](const Vector& v,Vector& y){
        std::vector<DomHiLoV55> direction(n);
        for(std::size_t i=0;i<n;++i)direction[i]=DomHiLoV55(v(i));
        const auto Av=stable.apply_hilo_v55(direction);
        if(y.size()!=n)y.resize(n);
        for(std::size_t i=0;i<n;++i){
            y(i)=Av[i].hi+Av[i].lo;
            if(!std::isfinite(y(i)))throw std::runtime_error("V56: bad Krylov matvec");
        }
    };
    FunctionalLinearOperator fine(n,op.apply);
    MatrixFreeVcyclePreconditioner amg(fine,.55,2,2,.25,12,
        AMGInterpolationPolicy::DirectionalSemiCoarsening,false,false,
        AMGSmootherPolicy::AlgebraicLineGaussSeidel);
    DomV56Linear result;
    if(!amg.setup(lo.auxiliary_matrix()))return result;
    result.amg_setup=true;
    auto x=initial;
    const auto true_residual=[&](const std::vector<DomHiLoV55>& state,
                                  std::vector<DomHiLoV55>* defects) {
        const auto applied=stable.apply_hilo_v55(state);
        long double sum=0.;
        if(defects)defects->resize(n);
        for(std::size_t i=0;i<n;++i){
            const auto r=rhs[i]-applied[i];
            if(!dom_finite_v55(r))throw std::runtime_error("V56: nonfinite true residual");
            const auto v=dom_to_long_double_v55(r);
            sum+=v*v;
            if(defects)(*defects)[i]=r;
        }
        return static_cast<double>(std::sqrt(sum/den));
    };
    result.initial=true_residual(x,nullptr);
    result.final=result.initial;
    for(unsigned update=0;update<cfg.max_updates && result.final>cfg.strict;++update){
        std::vector<DomHiLoV55> defects;
        const double before=true_residual(x,&defects);
        Vector b(n,0.),dx(n,0.);
        for(std::size_t i=0;i<n;++i)b(i)=defects[i].hi+defects[i].lo;
        const auto info=solve_gmres(op,b,dx,cfg.restart,cfg.inner_budget,
                                     cfg.inner_tolerance,&amg);
        result.krylov+=info.iterations;
        if(info.status!=SolverStatus::CONVERGED){++result.rejected_updates;break;}
        bool committed=false;
        for(double alpha=1.;alpha>=1./64.;alpha*=.5){
            auto candidate=x;
            bool physical=true;
            for(std::size_t i=0;i<n;++i){
                candidate[i]+=DomHiLoV55(alpha)*DomHiLoV55(dx(i));
                if(!dom_finite_v55(candidate[i]) || candidate[i].hi<0. ||
                   (candidate[i].hi==0. && candidate[i].lo<0.)) {
                    physical=false;break;
                }
            }
            if(!physical)continue;
            const auto after=true_residual(candidate,nullptr);
            if(std::isfinite(after) && after<before){
                x.swap(candidate);result.final=after;
                ++result.accepted_updates;committed=true;break;
            }
        }
        if(!committed){++result.rejected_updates;break;}
    }
    result.final=true_residual(x,nullptr);
    result.strict=result.final<=cfg.strict;
    if(result.strict)result.split=std::move(x);
    return result;
}

struct DomV56Options {
    DomV53Options v53;
    DomV55RefinementOptions linear;
    std::size_t simulated_partitions=1;
};
struct DomV56Result {
    bool experimental_certified=false,production_published=false;
    std::size_t refreshed=0,accepted=0,rejected=0;
    std::size_t lo_attempts=0,lo_strict=0,lo_accepted=0,picard_accepted=0;
    std::size_t krylov=0;
    double final_ho=std::numeric_limits<double>::infinity();
    double final_lo=std::numeric_limits<double>::infinity();
    double final_angular=std::numeric_limits<double>::infinity();
    double final_original_double_corrected=std::numeric_limits<double>::infinity();
    double last_linear=std::numeric_limits<double>::infinity();
    double largest_hidden_component=0.;
    double max_tensor_refresh=0.;
    double max_face_correction_refresh=0.;
    double last_identity=std::numeric_limits<double>::infinity();
    std::string diagnostic;
    std::vector<DomHiLoV55> experimental_state;
    std::vector<double> published_double_state; // deliberately never populated
    std::vector<double> ho_history,lo_history;
};

inline DomV56Result solve_dom_holo_hilo_v56(
    DomMomentSchurSystemV44& ho,
    const std::vector<DomHiLoV55>& initial,const DomV56Options& options={}){
    using H=DomHiLoV55;
    const auto& ctrl=options.v53.base.base;
    const auto nc=ho.n_cells();
    if(initial.size()!=ho.dimension() || ctrl.max_outer==0 ||
       ctrl.max_backtracks==0 || !(ctrl.physical_tolerance>0.) ||
       ctrl.physical_tolerance>=1. ||
       !(ctrl.required_contraction>0.) || ctrl.required_contraction>=1. ||
       options.simulated_partitions==0 || options.simulated_partitions>128 ||
       options.linear.strict>1e-13 || !(options.linear.strict>0.))
        throw std::invalid_argument("V56: invalid exterior / precision controls");
    DomV56Result out;
    auto accepted=initial;
    const auto make_lo=[&](const DomV50Snapshot& snap){
        std::vector<double> extinction(nc);
        for(std::size_t c=0;c<nc;++c)
            extinction[c]=ho.cell_absorption()[c]+ho.cell_scattering()[c];
        return DomSmmGrayWallV49(ho.source_mesh(),ho.source_geometry(),
            ho.directions(),ho.source_walls(),snap.eddington,
            extinction,ho.cell_absorption());
    };
    struct State {
        DomV50Snapshot snap;
        DomV54AngularAudit angular;
        DomV56FaceClosure closure;
        double ho=std::numeric_limits<double>::infinity();
        double lo=std::numeric_limits<double>::infinity();
        double original_double_corrected=std::numeric_limits<double>::infinity();
        double hidden=0.;
    };
    const auto refresh=[&](const std::vector<H>& x) -> State {
        std::vector<double> rounded(x.size());
        double hidden=0.;
        for(std::size_t i=0;i<x.size();++i){
            if(!dom_finite_v55(x[i]) || x[i].hi<0. ||
               (x[i].hi==0. && x[i].lo<0.))
                throw std::invalid_argument("V56: invalid split candidate");
            rounded[i]=x[i].hi;
            hidden=std::max(hidden,std::abs(x[i].lo)/
                std::max(1.,std::abs(x[i].hi)));
        }
        // HO transport still consumes rounded double J/R. A mathematically
        // meaningful low component is sub-ulp; never permit an O(1%) hidden
        // moment which the angular equations cannot see.
        if(hidden>ctrl.physical_tolerance)
            throw std::invalid_argument("V56: low component invisible to HO exceeds physical gate");
        State result;
        result.hidden=hidden;
        result.snap=evaluate_dom_variable_v50(ho,rounded);
        result.angular=audit_dom_angular_v54(ho,rounded,
            ctrl.physical_tolerance,options.simulated_partitions);
        const auto lo=make_lo(result.snap);
        std::vector<double> kt(nc);
        for(std::size_t c=0;c<nc;++c)
            kt[c]=ho.cell_absorption()[c]+ho.cell_scattering()[c];
        const DomStableWlsV54 stable(ho.source_mesh(),ho.source_geometry(),
            lo,result.snap.eddington,kt,ho.cell_absorption());
        std::vector<H> projected(nc),current(nc);
        for(std::size_t c=0;c<nc;++c){
            projected[c]=H(result.snap.projected[c]);current[c]=x[c];
        }
        result.closure=build_dom_hilo_face_closure_v56(ho.source_mesh(),
            ho.source_geometry(),lo,stable,projected,result.snap.gas_ib,
            ho.cell_absorption(),result.snap.face_current);
        const auto evaluated=stable.apply_hilo_v55(current);
        long double norm=0.,den=0.;
        for(std::size_t i=0;i<nc;++i){
            const auto r=evaluated[i]-result.closure.corrected_rhs[i];
            const auto v=dom_to_long_double_v55(r);
            const auto b=dom_to_long_double_v55(result.closure.corrected_rhs[i]);
            norm+=v*v;den+=b*b;
        }
        result.lo=den>0.?static_cast<double>(std::sqrt(norm/den)):
            std::numeric_limits<double>::infinity();
        // Do not quietly replace the historical V48 double residual with the
        // new compensated one. Rebuild the original correction from the same
        // current angular snapshot and evaluate original.apply(hi) separately.
        std::vector<double> projected_cells(result.snap.projected.begin(),
            result.snap.projected.begin()+static_cast<std::ptrdiff_t>(nc));
        const auto original_source=lo.consistent_source(projected_cells,
            result.snap.gas_ib,result.snap.face_current);
        std::vector<double> original_J(nc);
        for(std::size_t c=0;c<nc;++c)original_J[c]=rounded[c];
        const auto original_action=lo.apply(original_J);
        long double dnum=0.,dden=0.;
        for(std::size_t c=0;c<nc;++c){
            const long double diff=static_cast<long double>(original_action[c])-
                original_source.corrected_rhs_watts[c];
            const long double b=original_source.corrected_rhs_watts[c];
            dnum+=diff*diff;dden+=b*b;
        }
        result.original_double_corrected=dden>0.?
            static_cast<double>(std::sqrt(dnum/dden)):
            std::numeric_limits<double>::infinity();
        result.ho=std::max({result.snap.true_relative,
            result.angular.max_ordinate_equation,
            result.angular.max_wall_reflection,
            result.angular.max_cell_energy,
            result.angular.global_energy,
            result.angular.max_fixed_moment,
            result.angular.face_partition_error,hidden});
        ++out.refreshed;
        out.largest_hidden_component=std::max(out.largest_hidden_component,hidden);
        return result;
    };
    auto state=refresh(accepted);
    out.ho_history.push_back(state.ho);
    out.lo_history.push_back(state.lo);
    const auto certify=[&]() {
        if(out.lo_accepted==0 || out.last_linear>options.linear.strict ||
           !state.angular.accepted || state.ho>ctrl.physical_tolerance ||
           state.lo>ctrl.physical_tolerance ||
           state.closure.identity>ctrl.physical_tolerance ||
           state.hidden>ctrl.physical_tolerance)return false;
        // One more independent angular and freshly assembled WLS check is
        // performed by refresh() before setting this certificate.
        out.experimental_certified=true;
        out.experimental_state=accepted;
        out.diagnostic="V56: refreshed split/HO angular/gray-wall gates satisfied; CPU research only";
        return true;
    };
    for(std::size_t outer=0;outer<ctrl.max_outer &&
                           !out.experimental_certified;++outer){
        std::vector<H> picard(accepted.size());
        for(std::size_t i=0;i<accepted.size();++i)
            picard[i]=H(state.snap.projected[i])-accepted[i];
        const auto trial=[&](const std::vector<H>& step,bool is_lo) {
            for(unsigned t=0;t<ctrl.max_backtracks;++t){
                const H alpha(std::ldexp(1.,-static_cast<int>(t)));
                auto candidate=accepted;
                bool valid=true;
                for(std::size_t i=0;i<candidate.size();++i){
                    candidate[i]+=alpha*step[i];
                    if(!dom_finite_v55(candidate[i]) || candidate[i].hi<0. ||
                       (candidate[i].hi==0. && candidate[i].lo<0.)){
                        valid=false;break;
                    }
                }
                if(!valid){++out.rejected;continue;}
                try{
                    auto next=refresh(candidate);
                    const bool pass=next.angular.accepted &&
                        next.ho<=ctrl.physical_tolerance &&
                        next.lo<=ctrl.physical_tolerance;
                    const bool contraction=next.ho<ctrl.required_contraction*state.ho &&
                        (next.lo<ctrl.required_contraction*state.lo ||
                         next.lo<=ctrl.physical_tolerance) &&
                        next.closure.identity<=ctrl.physical_tolerance &&
                        next.hidden<=ctrl.physical_tolerance;
                    if(pass || contraction){
                        // Evidence that HO E tensors and HO-LO face corrections
                        // are rebuilt from every *accepted* angular state.
                        for(std::size_t i=0;i<state.snap.eddington.size();++i)
                            for(std::size_t k=0;k<6;++k)
                                out.max_tensor_refresh=std::max(out.max_tensor_refresh,
                                    std::abs(next.snap.eddington[i][k]-
                                             state.snap.eddington[i][k]));
                        for(std::size_t f=0;f<state.closure.ho_minus_lo_face.size();++f)
                            out.max_face_correction_refresh=std::max(
                                out.max_face_correction_refresh,
                                std::abs(static_cast<double>(dom_to_long_double_v55(
                                    next.closure.ho_minus_lo_face[f]-
                                    state.closure.ho_minus_lo_face[f]))));
                        accepted.swap(candidate);state=std::move(next);
                        ++out.accepted;
                        if(is_lo)++out.lo_accepted;
                        else ++out.picard_accepted;
                        out.ho_history.push_back(state.ho);
                        out.lo_history.push_back(state.lo);
                        if(pass && out.lo_accepted>0 &&
                           out.last_linear<=options.linear.strict){
                            state=refresh(accepted); // independent final recheck
                            (void)certify();
                        }
                        return true;
                    }
                }catch(const std::exception&) { /* invalid trial; accepted state intact */ }
                ++out.rejected;
            }
            return false;
        };
        bool committed=false;
        {
            const auto lo=make_lo(state.snap);
            std::vector<double> kt(nc);
            for(std::size_t c=0;c<nc;++c)
                kt[c]=ho.cell_absorption()[c]+ho.cell_scattering()[c];
            const DomStableWlsV54 stable(ho.source_mesh(),ho.source_geometry(),
                lo,state.snap.eddington,kt,ho.cell_absorption());
            ++out.lo_attempts;
            std::vector<H> before(accepted.begin(),accepted.begin()+
                                     static_cast<std::ptrdiff_t>(nc));
            const auto solve=solve_dom_hilo_smm_amg_v56(lo,stable,
                state.closure.corrected_rhs,before,options.linear);
            out.krylov+=solve.krylov;
            out.last_linear=solve.final;
            if(solve.strict){
                ++out.lo_strict;
                std::vector<H> step(accepted.size());
                for(std::size_t i=0;i<nc;++i)
                    step[i]=solve.split[i]-accepted[i];
                for(std::size_t i=nc;i<step.size();++i)
                    step[i]=picard[i];
                if(options.v53.base.use_wall_response && ho.n_reflecting_faces()){
                    auto probe=accepted;
                    for(std::size_t i=0;i<nc;++i)probe[i]=solve.split[i];
                    std::vector<double> rounded(probe.size());
                    for(std::size_t i=0;i<probe.size();++i)rounded[i]=probe[i].hi;
                    const auto response=ho.project(rounded,true);
                    for(std::size_t i=nc;i<step.size();++i)
                        step[i]=H(response[i])-accepted[i];
                }
                committed=trial(step,true);
            }
        }
        if(!committed && ctrl.allow_picard_fallback)
            committed=trial(picard,false);
        if(!committed){out.diagnostic="V56: true HO/LO noncontracting; rolled back";break;}
    }
    out.final_ho=state.ho;
    out.final_lo=state.lo;
    out.final_original_double_corrected=state.original_double_corrected;
    out.last_identity=state.closure.identity;
    out.final_angular=std::max({state.angular.max_ordinate_equation,
        state.angular.max_wall_reflection,state.angular.max_cell_energy,
        state.angular.global_energy,state.angular.max_fixed_moment,
        state.angular.face_partition_error});
    if(!out.experimental_certified && out.diagnostic.empty())
        out.diagnostic="V56: budget exhausted; no state published";
    out.production_published=false;
    out.published_double_state.clear();
    return out;
}
} // namespace cfdx::physics
