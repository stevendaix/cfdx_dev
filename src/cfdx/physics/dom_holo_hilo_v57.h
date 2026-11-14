#pragma once
// CFDX V57 research-only: angular intensities, HO face powers, moments and
// refreshed SMM correction are all computed in two-double arithmetic.
// The accepted state is the split vector, never a hidden double projection.
#include "cfdx/physics/dom_angular_hilo_v57.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <chrono>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct DomV57FaceClosure {
    std::vector<DomHiLoV55> corrected_rhs,ho_minus_lo_face;
    double identity=std::numeric_limits<double>::infinity();
    double physical_ho=std::numeric_limits<double>::infinity();
};

inline DomV57FaceClosure build_dom_hilo_face_closure_v57(
    const cfdx::core::Mesh& mesh,const FvGeometry& geo,
    const DomSmmGrayWallV49& lo,const DomStableWlsV54& stable,
    const std::vector<DomHiLoV55>& moments,const std::vector<double>& gas_ib,
    const std::vector<double>& ka,const std::vector<DomHiLoV55>& ho_faces) {
    using H=DomHiLoV55;
    const std::size_t nc=mesh.n_cells(),nf=mesh.n_faces();
    if(moments.size()!=nc || gas_ib.size()!=nc || ka.size()!=nc ||
       ho_faces.size()!=nf || lo.dimension()!=nc)
        throw std::invalid_argument("V57 closure: dimensions");
    const auto L=stable.evaluate_hilo_v55(moments);
    DomV57FaceClosure r;
    r.corrected_rhs.assign(nc,H(0.));
    r.ho_minus_lo_face.resize(nf);
    std::vector<H> wall_flux(nf),corrected_div(nc),ho_div(nc),source(nc),abs_ho(nc);
    for(const auto& w:lo.walls()){
        wall_flux[w.face]=H(w.conductance)*(moments[w.owner]-H(w.ib));
        source[w.owner]+=H(w.conductance)*H(w.ib);
    }
    for(std::size_t c=0;c<nc;++c){
        if(!(ka[c]>0.) || !dom_finite_v55(moments[c]) ||
           !(geo.cell_volumes[c]>0.))
            throw std::invalid_argument("V57 closure: invalid optical inputs");
        source[c]+=dom_four_pi_v55()*H(ka[c])*H(geo.cell_volumes[c])*H(gas_ib[c]);
    }
    for(std::size_t f=0;f<nf;++f){
        if(!dom_finite_v55(ho_faces[f]))
            throw std::invalid_argument("V57 closure: nonfinite HO flux");
        const auto o=mesh.ownership().owner(f);
        const auto nb=mesh.ownership().neighbour(f);
        const H qlo=nb>=0?L.face_watts[f]:wall_flux[f];
        const H delta=ho_faces[f]-qlo;
        r.ho_minus_lo_face[f]=delta;
        corrected_div[o]+=delta;
        ho_div[o]+=ho_faces[f];
        abs_ho[o]+=dom_abs_v57(ho_faces[f]);
        if(nb>=0){
            const auto j=static_cast<std::size_t>(nb);
            corrected_div[j]-=delta;
            ho_div[j]-=ho_faces[f];
            abs_ho[j]+=dom_abs_v57(ho_faces[f]);
        }
    }
    r.identity=0.;r.physical_ho=0.;
    for(std::size_t c=0;c<nc;++c){
        r.corrected_rhs[c]=source[c]-corrected_div[c];
        const H algebra=L.divergence_watts[c]-r.corrected_rhs[c];
        const H physical=ho_div[c]+dom_four_pi_v55()*H(ka[c])*
            H(geo.cell_volumes[c])*(moments[c]-H(gas_ib[c]));
        const H identity=algebra-physical;
        const long double scale=std::max({1.L,
            std::abs(dom_to_long_double_v55(abs_ho[c])),
            std::abs(dom_to_long_double_v55(physical))});
        r.identity=std::max(r.identity,
            static_cast<double>(std::abs(dom_to_long_double_v55(identity))/scale));
        r.physical_ho=std::max(r.physical_ho,
            static_cast<double>(std::abs(dom_to_long_double_v55(physical))/scale));
    }
    return r;
}

struct DomV57Options {
    DomV56Options v56;
    // A non-monotone component may decrease the joint max-norm; the accepted
    // merit must still strictly contract, and ALL physical gates apply at exit.
    bool allow_merit_globalization=true;
    std::function<void(std::size_t,double,double)> progress; // optional audit trace
};
struct DomV57Result {
    bool experimental_certified=false,production_published=false;
    std::size_t accepted=0,rejected=0,refreshed=0;
    std::size_t lo_attempts=0,lo_strict=0,lo_accepted=0,picard_accepted=0;
    std::size_t krylov=0,angular_sweeps=0;
    double final_ho=std::numeric_limits<double>::infinity();
    double final_lo=std::numeric_limits<double>::infinity();
    double final_angular=std::numeric_limits<double>::infinity();
    double final_energy=std::numeric_limits<double>::infinity();
    double last_linear=std::numeric_limits<double>::infinity();
    double original_double_corrected=std::numeric_limits<double>::infinity();
    double final_identity=std::numeric_limits<double>::infinity();
    double max_tensor_refresh=0.,max_face_refresh=0.;
    double largest_low=0.,max_angular_csr=0.;
    double angular_seconds=0.,refresh_seconds=0.,amg_seconds=0.;
    std::string diagnostic;
    std::vector<DomHiLoV55> certified_split_state;
    std::vector<double> production_double_state; // intentionally always empty
    std::vector<double> ho_history,lo_history,energy_history,merit_history;
};

inline DomV57Result solve_dom_holo_hilo_v57(
    DomMomentSchurSystemV44& ho,
    const std::vector<DomHiLoV55>& initial,
    const DomV57Options& options={}) {
    using H=DomHiLoV55;
    const auto& opt=options.v56;
    const auto& ctrl=opt.v53.base.base;
    const std::size_t nc=ho.n_cells();
    if(initial.size()!=ho.dimension() || ctrl.max_outer==0 ||
       ctrl.max_backtracks==0 || !(ctrl.physical_tolerance>0.) ||
       ctrl.physical_tolerance>1e-9 ||
       !(ctrl.required_contraction>0.) || ctrl.required_contraction>=1. ||
       opt.linear.strict>1e-13 || !(opt.linear.strict>0.) ||
       opt.simulated_partitions==0 || opt.simulated_partitions>128)
        throw std::invalid_argument("V57: invalid/loosened physical or linear controls");
    DomAngularHiloTransportV57 angular_transport(ho);
    DomV57Result out;
    auto accepted=initial;
    std::vector<double> kt(nc);
    for(std::size_t c=0;c<nc;++c)kt[c]=ho.cell_absorption()[c]+ho.cell_scattering()[c];
    struct State {
        DomV57AngularSolution ho;
        DomV57IndependentAudit audit;
        DomV57FaceClosure closure;
        double physical=std::numeric_limits<double>::infinity();
        double low=std::numeric_limits<double>::infinity();
        double merit=std::numeric_limits<double>::infinity();
    };
    const auto make_lo=[&](const DomV57AngularSolution& s){
        return DomSmmGrayWallV49(ho.source_mesh(),ho.source_geometry(),
            ho.directions(),ho.source_walls(),s.eddington,kt,ho.cell_absorption());
    };
    const auto refresh=[&](const std::vector<H>& x)->State {
        const auto refresh_start=std::chrono::steady_clock::now();
        State s;
        const auto angular_start=std::chrono::steady_clock::now();
        s.ho=angular_transport.solve(x);
        out.angular_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-angular_start).count();
        ++out.angular_sweeps;
        s.audit=audit_dom_angular_hilo_v57(ho,x,s.ho,ctrl.physical_tolerance);
        const auto lo=make_lo(s.ho);
        const DomStableWlsV54 stable(ho.source_mesh(),ho.source_geometry(),
            lo,s.ho.eddington,kt,ho.cell_absorption());
        std::vector<H> projected(s.ho.moment.begin(),
            s.ho.moment.begin()+static_cast<std::ptrdiff_t>(nc));
        s.closure=build_dom_hilo_face_closure_v57(ho.source_mesh(),
            ho.source_geometry(),lo,stable,projected,s.ho.gas_ib,
            ho.cell_absorption(),s.ho.face_current);
        std::vector<H> current(x.begin(),x.begin()+static_cast<std::ptrdiff_t>(nc));
        const auto applied=stable.apply_hilo_v55(current);
        long double num=0.,den=0.;
        for(std::size_t i=0;i<nc;++i){
            const auto r=applied[i]-s.closure.corrected_rhs[i];
            const long double a=dom_to_long_double_v55(r);
            const long double b=dom_to_long_double_v55(s.closure.corrected_rhs[i]);
            num+=a*a;den+=b*b;
        }
        s.low=den>0.?static_cast<double>(std::sqrt(num/den)):
            std::numeric_limits<double>::infinity();
        s.physical=std::max({s.audit.max_ordinate_equation,
            s.audit.max_wall_reflection,s.audit.max_cell_energy,
            s.audit.global_energy,s.audit.max_fixed_moment,
            s.audit.max_face_recompute,s.ho.max_true_csr_relative,
            s.closure.physical_ho});
        s.merit=std::max(s.physical,s.low);
        if(!std::isfinite(s.merit) || s.closure.identity>1e-13 ||
           s.ho.max_true_csr_relative>1e-13)
            throw std::runtime_error("V57: invalid refreshed angular/face ledger");
        out.max_angular_csr=std::max(out.max_angular_csr,s.ho.max_true_csr_relative);
        out.largest_low=std::max(out.largest_low,s.ho.max_noncanonical_low);
        ++out.refreshed;
        out.refresh_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-refresh_start).count();
        return s;
    };
    auto state=refresh(accepted);
    const auto record=[&](){
        out.ho_history.push_back(state.physical);
        out.lo_history.push_back(state.low);
        out.energy_history.push_back(state.audit.global_energy);
        out.merit_history.push_back(state.merit);
    };
    record();
    const auto certify=[&]()->bool {
        if(out.lo_accepted==0 || out.last_linear>opt.linear.strict ||
           !state.audit.accepted || state.physical>ctrl.physical_tolerance ||
           state.low>ctrl.physical_tolerance ||
           state.closure.identity>1e-13 ||
           state.ho.max_true_csr_relative>1e-13)return false;
        // Explicitly retain the original double angular oracle as a separate
        // regression; it cannot certify the split flow by itself.
        std::vector<double> rounded(accepted.size());
        for(std::size_t i=0;i<rounded.size();++i)rounded[i]=accepted[i].hi;
        const auto double_audit=audit_dom_angular_v54(ho,rounded,
            ctrl.physical_tolerance,opt.simulated_partitions);
        if(!double_audit.accepted){out.diagnostic="V57: double V54 cross-audit rejected";return false;}
        // A final independent angular sweep plus LO assembly on THE stored
        // vector; never rely on a cached accepted-trial score.
        state=refresh(accepted);
        if(!state.audit.accepted || state.low>ctrl.physical_tolerance ||
           state.physical>ctrl.physical_tolerance ||
           state.closure.identity>1e-13)return false;
        out.experimental_certified=true;
        out.certified_split_state=accepted;
        out.diagnostic="V57: split HO and LO independently certified; research CPU only";
        return true;
    };
    for(std::size_t outer=0;outer<ctrl.max_outer &&
         !out.experimental_certified;++outer){
        std::vector<H> picard(accepted.size());
        for(std::size_t i=0;i<accepted.size();++i)
            picard[i]=state.ho.moment[i]-accepted[i];
        const auto trial=[&](const std::vector<H>& step,bool is_lo) {
            for(unsigned backtrack=0;backtrack<ctrl.max_backtracks;++backtrack){
                const double alpha=std::ldexp(1.,-static_cast<int>(backtrack));
                auto candidate=accepted;
                bool positive=true;
                for(std::size_t i=0;i<candidate.size();++i){
                    candidate[i]+=H(alpha)*step[i];
                    if(!dom_finite_v55(candidate[i]) || candidate[i].hi<0. ||
                       (candidate[i].hi==0. && candidate[i].lo<0.)){
                        positive=false;break;
                    }
                }
                if(!positive){++out.rejected;continue;}
                try{
                    auto next=refresh(candidate);
                    const bool all=next.audit.accepted &&
                        next.physical<=ctrl.physical_tolerance &&
                        next.low<=ctrl.physical_tolerance;
                    const bool individual=next.physical<ctrl.required_contraction*state.physical &&
                        (next.low<ctrl.required_contraction*state.low ||
                         next.low<=ctrl.physical_tolerance);
                    const bool merit=options.allow_merit_globalization &&
                        next.merit<(1.-1e-4*alpha)*state.merit;
                    if(all || individual || merit){
                        for(std::size_t i=0;i<nc;++i)
                            for(unsigned k=0;k<6;++k)
                                out.max_tensor_refresh=std::max(out.max_tensor_refresh,
                                    std::abs(next.ho.eddington[i][k]-
                                             state.ho.eddington[i][k]));
                        for(std::size_t f=0;f<state.closure.ho_minus_lo_face.size();++f)
                            out.max_face_refresh=std::max(out.max_face_refresh,
                                static_cast<double>(std::abs(dom_to_long_double_v55(
                                    next.closure.ho_minus_lo_face[f]-
                                    state.closure.ho_minus_lo_face[f]))));
                        accepted.swap(candidate);
                        state=std::move(next);
                        ++out.accepted;
                        if(is_lo)++out.lo_accepted;
                        else ++out.picard_accepted;
                        record();
                        if(options.progress)options.progress(out.accepted,state.physical,state.low);
                        if(all)(void)certify();
                        return true;
                    }
                }catch(const std::exception&){/* transaction rejects malformed trial */}
                ++out.rejected;
            }
            return false;
        };
        bool committed=false;
        {
            const auto lo=make_lo(state.ho);
            const DomStableWlsV54 stable(ho.source_mesh(),ho.source_geometry(),
                lo,state.ho.eddington,kt,ho.cell_absorption());
            ++out.lo_attempts;
            std::vector<H> before(accepted.begin(),
                accepted.begin()+static_cast<std::ptrdiff_t>(nc));
            const auto linear_start=std::chrono::steady_clock::now();
            const auto solved=solve_dom_hilo_smm_amg_v56(lo,stable,
                state.closure.corrected_rhs,before,opt.linear);
            out.amg_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-linear_start).count();
            out.last_linear=solved.final;
            out.krylov+=solved.krylov;
            if(solved.strict){
                ++out.lo_strict;
                std::vector<H> step(accepted.size());
                for(std::size_t i=0;i<nc;++i)
                    step[i]=solved.split[i]-accepted[i];
                for(std::size_t i=nc;i<step.size();++i)
                    step[i]=picard[i];
                // The wall response is also transported in hi+lo, rather
                // than V56's rounded-double ho.project() call.
                if(opt.v53.base.use_wall_response && ho.n_reflecting_faces()){
                    auto probed=accepted;
                    for(std::size_t i=0;i<nc;++i)probed[i]=solved.split[i];
                    const auto response=angular_transport.solve(probed);
                    ++out.angular_sweeps;
                    for(std::size_t i=nc;i<step.size();++i)
                        step[i]=response.moment[i]-accepted[i];
                }
                committed=trial(step,true);
            }
        }
        if(!committed && ctrl.allow_picard_fallback)
            committed=trial(picard,false);
        if(!committed){out.diagnostic="V57: no accepted true HO/LO merit step";break;}
    }
    out.final_ho=state.physical;
    out.final_lo=state.low;
    out.final_angular=std::max({state.audit.max_ordinate_equation,
        state.audit.max_wall_reflection,state.audit.max_fixed_moment,
        state.audit.max_face_recompute,state.audit.csr_residual});
    out.final_energy=std::max(state.audit.max_cell_energy,state.audit.global_energy);
    out.final_identity=state.closure.identity;
    // Final ORIGINAL V48 double residual is diagnostic only. Its source is
    // explicitly the current rounded split HO flux, not an 80-bit result.
    const auto lo=make_lo(state.ho);
    std::vector<double> original_mom(nc),original_x(nc),original_face(state.ho.face_current.size());
    for(std::size_t c=0;c<nc;++c){
        original_mom[c]=state.ho.moment[c].hi;original_x[c]=accepted[c].hi;
    }
    for(std::size_t f=0;f<original_face.size();++f)
        original_face[f]=state.ho.face_current[f].hi;
    const auto source=lo.consistent_source(original_mom,state.ho.gas_ib,original_face);
    const auto action=lo.apply(original_x);
    long double nr=0.,dr=0.;
    for(std::size_t i=0;i<nc;++i){
        const long double e=static_cast<long double>(action[i])-source.corrected_rhs_watts[i];
        const long double b=source.corrected_rhs_watts[i];
        nr+=e*e;dr+=b*b;
    }
    out.original_double_corrected=dr>0.?static_cast<double>(std::sqrt(nr/dr)):
        std::numeric_limits<double>::infinity();
    if(out.diagnostic.empty())out.diagnostic="V57: budget reached without all physical gates";
    out.production_published=false;
    out.production_double_state.clear();
    return out;
}

} // namespace cfdx::physics
