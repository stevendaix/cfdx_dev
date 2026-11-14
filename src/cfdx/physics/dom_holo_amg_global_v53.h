#pragma once
// V53 research-only, no production DOM changes. Coupled angular HO +
// continually rebuilt LO second-moment predictor with directional AMG.
#include "cfdx/physics/dom_holo_global_v51.h"
#include "cfdx/physics/dom_smm_anisotropic_gate_v52.h"
#include <deque>
#include <limits>
#include <string>
#include <vector>
#include <cmath>

namespace cfdx::physics {
struct DomV53Options {
    DomV51Options base;
    double strict_linear=1e-13;
    std::size_t max_krylov=300, wide_budget=450;
    int restart=30,wide_restart=120;
    bool wide_retry=true,require_amg_step=true;
    unsigned iterative_refinements=1;
    cfdx::core::AMGInterpolationPolicy interpolation=
        cfdx::core::AMGInterpolationPolicy::DirectionalSemiCoarsening;
    cfdx::core::AMGSmootherPolicy smoother=
        cfdx::core::AMGSmootherPolicy::AlgebraicLineGaussSeidel;
};
struct DomV53Linear {
    bool strict=false, wide_used=false;
    double true_l2=std::numeric_limits<double>::infinity();
    std::size_t iterations=0,levels=0;
    std::vector<double> J;
};
struct DomV53Audit {
    bool ho_pass=false,lo_pass=false;
    double ho=std::numeric_limits<double>::infinity();
    double lo=std::numeric_limits<double>::infinity();
    double global=std::numeric_limits<double>::infinity();
    double identity=std::numeric_limits<double>::infinity();
    double reflector=std::numeric_limits<double>::infinity();
};
struct DomV53Result {
    bool converged=false,coupled_amg=false;
    std::size_t accepted=0,lo_attempts=0,lo_accepted=0,amg_strict=0;
    std::size_t amg_failed=0,anderson_accepted=0,picard_accepted=0,rejected=0;
    std::size_t refreshes=0;
    double final_ho=std::numeric_limits<double>::infinity();
    double final_lo=std::numeric_limits<double>::infinity();
    double final_global=std::numeric_limits<double>::infinity();
    double last_linear=std::numeric_limits<double>::infinity();
    std::string reason;
    std::vector<double> state,ho_history,lo_history;
};

inline DomV53Linear solve_smm_amg_corrected_v53(
    const DomSmmGrayWallV49& lo,const std::vector<double>& rhs,
    const std::vector<double>& initial,const DomV53Options& options) {
    using namespace cfdx::core;
    const auto n=lo.dimension();
    if(n==0 || rhs.size()!=n || initial.size()!=n ||
       options.max_krylov==0 || options.restart<2 ||
       options.wide_budget==0 || options.wide_restart<2 ||
       options.iterative_refinements>3 ||
       !(options.strict_linear>0.) || !std::isfinite(options.strict_linear))
        throw std::invalid_argument("V53: malformed LO system or controls");
    LinearOperator op;
    op.size=n;
    op.apply=[&lo,n](const Vector& v,Vector& y){
        std::vector<double> a(n);
        for(std::size_t i=0;i<n;++i)a[i]=v(i);
        const auto Av=lo.apply(a);
        if(y.size()!=n)y.resize(n);
        for(std::size_t i=0;i<n;++i)y(i)=Av[i];
    };
    FunctionalLinearOperator fine(n,op.apply);
    MatrixFreeVcyclePreconditioner amg(fine,.55,2,2,.25,12,
        options.interpolation,false,false,options.smoother);
    DomV53Linear result;
    if(!amg.setup(lo.auxiliary_matrix()))return result;
    result.levels=amg.hierarchy_level_sizes().size();
    Vector b(n,0.),x(n,0.);
    for(std::size_t i=0;i<n;++i){
        if(!std::isfinite(rhs[i]) || !std::isfinite(initial[i]))
            throw std::invalid_argument("V53: nonfinite RHS or initial moment");
        b(i)=rhs[i];x(i)=initial[i];
    }
    const auto audit=[&](){
        std::vector<double> v(n);
        for(std::size_t i=0;i<n;++i)v[i]=x(i);
        const auto Ax=lo.apply(v);
        long double num=0.,den=0.;
        for(std::size_t i=0;i<n;++i){
            const long double r=static_cast<long double>(rhs[i])-Ax[i];
            num+=r*r;den+=static_cast<long double>(rhs[i])*rhs[i];
        }
        return den>0.?std::sqrt(static_cast<double>(num/den)):
            std::numeric_limits<double>::infinity();
    };
    const auto first=solve_gmres(op,b,x,options.restart,options.max_krylov,
                                 options.strict_linear,&amg);
    result.iterations+=first.iterations;
    result.true_l2=audit();
    result.strict=first.status==SolverStatus::CONVERGED &&
                  result.true_l2<=options.strict_linear;
    for(unsigned j=0;j<options.iterative_refinements && !result.strict;++j){
        std::vector<double> v(n);
        for(std::size_t i=0;i<n;++i)v[i]=x(i);
        const auto Av=lo.apply(v);
        Vector r(n,0.),dx(n,0.);
        for(std::size_t i=0;i<n;++i)r(i)=rhs[i]-Av[i];
        // Internal correction tolerance does not change the strict final gate.
        const auto inner=solve_gmres(op,r,dx,options.restart,
                                     options.max_krylov,1e-3,&amg);
        result.iterations+=inner.iterations;
        if(inner.status!=SolverStatus::CONVERGED)break;
        for(std::size_t i=0;i<n;++i)x(i)+=dx(i);
        result.true_l2=audit();
        result.strict=result.true_l2<=options.strict_linear;
    }
    if(!result.strict && options.wide_retry &&
       (options.wide_restart>options.restart ||
        options.wide_budget>options.max_krylov)){
        result.wide_used=true;
        for(std::size_t i=0;i<n;++i)x(i)=initial[i];
        const auto widened=solve_gmres(op,b,x,options.wide_restart,
            options.wide_budget,options.strict_linear,&amg);
        result.iterations+=widened.iterations;
        result.true_l2=audit();
        result.strict=widened.status==SolverStatus::CONVERGED &&
            result.true_l2<=options.strict_linear;
    }
    if(result.strict){
        result.J.resize(n);
        for(std::size_t i=0;i<n;++i){
            result.J[i]=x(i);
            if(!std::isfinite(result.J[i]))result.strict=false;
        }
        if(!result.strict)result.J.clear();
    }
    return result;
}

// The LO cell balance is algebraically linked to HO face currents, NOT an
// independent physical oracle; HO angular transport + gray-wall Schur +
// global boundary energy are independently evaluated in snapshot.
inline DomV53Audit audit_holo_amg_v53(const DomMomentSchurSystemV44& ho,
    const DomV50Snapshot& snap,const DomSmmGrayWallV49& lo,double tol){
    const auto nc=ho.n_cells();
    std::vector<double> J(nc);
    for(std::size_t i=0;i<nc;++i)J[i]=snap.projected[i];
    const auto src=lo.consistent_source(J,snap.gas_ib,snap.face_current);
    const auto Ax=lo.apply(J);
    DomV53Audit a;
    a.ho=snap.true_relative;
    a.identity=src.max_discrete_identity_relative;
    a.lo=0.;
    constexpr double fourpi=12.56637061435917295385;
    long double global=0.,scale=1.;
    for(std::size_t i=0;i<nc;++i){
        const double emission=fourpi*ho.cell_absorption()[i]*
            ho.source_geometry().cell_volumes[i]*(J[i]-snap.gas_ib[i]);
        const double cell=Ax[i]-src.corrected_rhs_watts[i];
        a.lo=std::max(a.lo,std::abs(cell)/
            std::max({1.,std::abs(snap.cell_divergence[i]),std::abs(emission)}));
        global+=emission;scale+=std::abs(emission);
    }
    const auto& mesh=ho.source_mesh();
    for(std::size_t f=0;f<mesh.n_faces();++f)
        if(mesh.ownership().neighbour(f)<0){
            global+=snap.face_current[f];scale+=std::abs(snap.face_current[f]);
        }
    a.global=std::abs(static_cast<double>(global/scale));
    a.reflector=0.;
    for(const auto& w:ho.wall_faces())if(w.emissivity==0.){
        const double area=ho.source_geometry().face_area_vectors[w.face].mag();
        const double flux_scale=std::max(1.,fourpi*area*
            std::max(J[w.owner],w.ib));
        a.reflector=std::max(a.reflector,
                             std::abs(snap.face_current[w.face])/flux_scale);
    }
    a.ho_pass=snap.fixed_relative<=tol && snap.transport_relative<=tol &&
        snap.physical_cell_relative<=tol && snap.global_energy_relative<=tol &&
        snap.maximum_trace_error<=tol;
    a.lo_pass=a.lo<=tol && a.global<=tol && a.identity<=tol && a.reflector<=tol;
    return a;
}

inline DomV53Result solve_dom_holo_amg_v53(
    DomMomentSchurSystemV44& ho,const std::vector<double>& initial,
    const DomV53Options& options={}){
    const auto& controls=options.base.base;
    if(initial.size()!=ho.dimension() || controls.max_outer==0 ||
       controls.max_backtracks==0 || !(controls.physical_tolerance>0.) ||
       !(controls.required_contraction>0.) || controls.required_contraction>=1. ||
       options.base.anderson_depth>12)
        throw std::invalid_argument("V53: invalid HO global controls");
    DomV53Result result;
    auto accepted=initial;
    const auto nc=ho.n_cells();
    auto create_lo=[&](const DomV50Snapshot& snap){
        std::vector<double> kt(nc);
        for(std::size_t i=0;i<nc;++i)
            kt[i]=ho.cell_absorption()[i]+ho.cell_scattering()[i];
        return DomSmmGrayWallV49(ho.source_mesh(),ho.source_geometry(),
            ho.directions(),ho.source_walls(),snap.eddington,kt,
            ho.cell_absorption());
    };
    struct State {DomV50Snapshot ho;DomV53Audit dual;};
    const auto refresh=[&](const std::vector<double>& x){
        auto snap=evaluate_dom_variable_v50(ho,x);
        const auto lo=create_lo(snap);
        auto audit=audit_holo_amg_v53(ho,snap,lo,controls.physical_tolerance);
        ++result.refreshes;
        return State{std::move(snap),std::move(audit)};
    };
    auto state=refresh(accepted);
    result.ho_history.push_back(state.dual.ho);
    result.lo_history.push_back(state.dual.lo);
    std::deque<dom_v51_detail::Secant> history;
    const auto publish=[&](){
        if(!state.dual.ho_pass || !state.dual.lo_pass ||
           (options.require_amg_step && !result.lo_accepted))return false;
        result.converged=true;
        result.coupled_amg=result.lo_accepted>0 && result.amg_strict>0;
        result.state=accepted;
        result.reason="both true HO physical and corrected LO balances accepted";
        return true;
    };
    (void)publish();
    for(std::size_t outer=0;outer<controls.max_outer && !result.converged;++outer){
        std::vector<double> picard(accepted.size());
        for(std::size_t i=0;i<accepted.size();++i)
            picard[i]=state.ho.projected[i]-accepted[i];
        history.push_back({accepted,picard});
        while(history.size()>options.base.anderson_depth+1)history.pop_front();
        const auto trial=[&](const std::vector<double>& direction,int kind){
            if(direction.size()!=accepted.size())return false;
            for(unsigned k=0;k<controls.max_backtracks;++k){
                const double alpha=std::ldexp(1.,-static_cast<int>(k));
                auto candidate=accepted;
                bool finite=true;
                for(std::size_t i=0;i<candidate.size();++i){
                    candidate[i]+=alpha*direction[i];
                    if(!std::isfinite(candidate[i]) || candidate[i]<0.)finite=false;
                }
                if(!finite){++result.rejected;continue;}
                try{
                    auto next=refresh(candidate);
                    const auto contraction=controls.required_contraction;
                    const bool improving=next.dual.ho<contraction*state.dual.ho &&
                        (next.dual.lo<contraction*state.dual.lo ||
                         next.dual.lo<=controls.physical_tolerance) &&
                        (next.dual.reflector<contraction*state.dual.reflector ||
                         next.dual.reflector<=controls.physical_tolerance) &&
                         next.dual.identity<=controls.physical_tolerance;
                    if(improving || (next.dual.ho_pass && next.dual.lo_pass)){
                        accepted=std::move(candidate);state=std::move(next);
                        ++result.accepted;
                        if(kind==0)++result.lo_accepted;
                        else if(kind==1)++result.anderson_accepted;
                        else ++result.picard_accepted;
                        result.ho_history.push_back(state.dual.ho);
                        result.lo_history.push_back(state.dual.lo);
                        (void)publish();
                        return true;
                    }
                }catch(const std::runtime_error&){/* candidate was invalid */}
                ++result.rejected;
            }
            return false;
        };
        bool committed=false;
        {
            const auto lo=create_lo(state.ho);
            std::vector<double> J(nc),begin(nc);
            for(std::size_t i=0;i<nc;++i){
                J[i]=state.ho.projected[i];begin[i]=accepted[i];
            }
            const auto corrected=lo.consistent_source(J,state.ho.gas_ib,
                                                       state.ho.face_current);
            if(corrected.max_discrete_identity_relative<=controls.physical_tolerance){
                ++result.lo_attempts;
                const auto low=solve_smm_amg_corrected_v53(lo,
                    corrected.corrected_rhs_watts,begin,options);
                result.last_linear=low.true_l2;
                if(low.strict){
                    ++result.amg_strict;
                    std::vector<double> step(accepted.size());
                    for(std::size_t i=0;i<nc;++i)step[i]=low.J[i]-accepted[i];
                    for(std::size_t i=nc;i<step.size();++i)step[i]=picard[i];
                    if(options.base.use_wall_response && ho.n_reflecting_faces()){
                        auto probe=accepted;
                        for(std::size_t i=0;i<nc;++i)probe[i]=low.J[i];
                        const auto response=ho.project(probe,true);
                        for(std::size_t i=nc;i<step.size();++i)
                            step[i]=response[i]-accepted[i];
                    }
                    committed=trial(step,0);
                }else ++result.amg_failed;
            }
        }
        if(!committed && options.base.use_anderson &&
           options.base.anderson_depth){
            const auto step=dom_v51_detail::direction(history,
                options.base.anderson_ridge,options.base.max_anderson_step_ratio);
            if(!step.empty())committed=trial(step,1);
        }
        if(!committed && controls.allow_picard_fallback)
            committed=trial(picard,2);
        if(!committed){result.reason="non-contracting HO/LO trial; rollback";break;}
    }
    result.final_ho=state.dual.ho;
    result.final_lo=state.dual.lo;
    result.final_global=state.dual.global;
    if(!result.converged){
        result.state.clear();
        if(result.reason.empty())result.reason="budget exceeded; no state published";
    }
    return result;
}
} // namespace cfdx::physics
