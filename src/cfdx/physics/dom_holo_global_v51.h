#pragma once
// V51 research-only: safeguarded joint gas / gray-wall HOLO acceleration.
// The historical V50 driver and production DOM remain completely unchanged.
// A candidate is accepted exclusively using the independent V50 full angular
// FV residual, including the same reflection, transport and global energy gates.
#include "cfdx/physics/dom_holo_variable_v50.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct DomV51Options {
    DomV50Options base;
    unsigned anderson_depth=5;
    double anderson_ridge=1.e-10;
    double max_anderson_step_ratio=30.;
    bool use_anderson=true;
    bool use_wall_response=true;
    bool use_low_order=true;
};

struct DomV51Result {
    bool converged=false;
    std::size_t accepted_steps=0, rejected_trials=0, lo_solves=0, lo_failures=0;
    std::size_t lo_steps=0, anderson_steps=0, picard_steps=0;
    std::size_t wall_response_maps=0, transport_maps=0, acceleration_restarts=0;
    double initial_true=std::numeric_limits<double>::infinity();
    double final_true=std::numeric_limits<double>::infinity();
    double seconds=0.;
    std::string diagnostic;
    std::vector<double> final_moments;
    std::vector<double> accepted_defect_history;
};

namespace dom_v51_detail {
struct Secant {
    std::vector<double> x,f;
};

inline double dot(const std::vector<double>& x,const std::vector<double>& y) {
    if(x.size()!=y.size()) throw std::invalid_argument("V51 Anderson: dot dimension");
    long double s=0.;
    for(std::size_t i=0;i<x.size();++i)s+=static_cast<long double>(x[i])*y[i];
    return static_cast<double>(s);
}

// Anderson II, with a relative Tikhonov ridge and pivot rejection. The full
// x=(J,R) and f=T(x)-x participate. No angle or wall equation is discarded.
// Returning an empty vector means that Picard must be considered instead.
inline std::vector<double> direction(const std::deque<Secant>& history,
                                     double ridge,double max_factor) {
    if(history.size()<2)return {};
    const auto n=history.back().x.size();
    const std::size_t m=history.size()-1;
    const auto& f=history.back().f;
    if(n==0 || f.size()!=n || m>12)return {};
    std::vector<std::vector<double>> dX(m,std::vector<double>(n));
    std::vector<std::vector<double>> dF(m,std::vector<double>(n));
    for(std::size_t j=0;j<m;++j) {
        if(history[j].x.size()!=n || history[j].f.size()!=n ||
           history[j+1].x.size()!=n || history[j+1].f.size()!=n)return {};
        for(std::size_t i=0;i<n;++i) {
            dX[j][i]=history[j+1].x[i]-history[j].x[i];
            dF[j][i]=history[j+1].f[i]-history[j].f[i];
        }
    }
    std::vector<std::vector<double>> a(m,std::vector<double>(m+1,0.));
    double diagonal_max=0.;
    for(std::size_t j=0;j<m;++j) {
        a[j][m]=dot(dF[j],f);
        for(std::size_t k=0;k<m;++k)a[j][k]=dot(dF[j],dF[k]);
        diagonal_max=std::max(diagonal_max,a[j][j]);
    }
    if(!(diagonal_max>0.) || !std::isfinite(diagonal_max))return {};
    const double regularization=ridge*diagonal_max;
    for(std::size_t j=0;j<m;++j)a[j][j]+=regularization;
    // Small regularized normal equations: partial pivot with finite gates.
    for(std::size_t k=0;k<m;++k) {
        std::size_t pivot=k;
        for(std::size_t row=k+1;row<m;++row)
            if(std::abs(a[row][k])>std::abs(a[pivot][k]))pivot=row;
        if(!(std::abs(a[pivot][k])>1.e-13*diagonal_max) ||
           !std::isfinite(a[pivot][k]))return {};
        if(pivot!=k)std::swap(a[pivot],a[k]);
        for(std::size_t row=k+1;row<m;++row) {
            const double factor=a[row][k]/a[k][k];
            for(std::size_t col=k+1;col<=m;++col)
                a[row][col]-=factor*a[k][col];
            a[row][k]=0.;
        }
    }
    std::vector<double> gamma(m);
    for(std::size_t row=m;row-->0;) {
        double rhs=a[row][m];
        for(std::size_t k=row+1;k<m;++k)rhs-=a[row][k]*gamma[k];
        gamma[row]=rhs/a[row][row];
        if(!std::isfinite(gamma[row]))return {};
    }
    std::vector<double> step=f;
    for(std::size_t j=0;j<m;++j)
        for(std::size_t i=0;i<n;++i)
            step[i]-=gamma[j]*(dX[j][i]+dF[j][i]);
    const double step_norm=std::sqrt(std::max(0.,dot(step,step)));
    const double picard_norm=std::sqrt(std::max(0.,dot(f,f)));
    if(!std::isfinite(step_norm) || !(picard_norm>0.) ||
       step_norm>max_factor*picard_norm)return {};
    return step;
}
} // namespace dom_v51_detail

// Transactional: final_moments is populated ONLY when every physical HO gate
// is satisfied. All trial state, E, source correction and wall R are local.
inline DomV51Result solve_dom_holo_global_v51(
        DomMomentSchurSystemV44& ho,const std::vector<double>& initial,
        const DomV51Options& o={}) {
    const auto& opt=o.base;
    if(initial.size()!=ho.dimension() || opt.max_outer==0 ||
       opt.max_lo_krylov==0 || opt.max_backtracks==0 ||
       !(opt.physical_tolerance>0.) || !std::isfinite(opt.physical_tolerance) ||
       !(opt.linear_tolerance>0.) || !std::isfinite(opt.linear_tolerance) ||
       !(opt.required_contraction>0.) || !(opt.required_contraction<1.) ||
       o.anderson_depth>12 || !(o.anderson_ridge>=0.) ||
       !std::isfinite(o.anderson_ridge) || !(o.max_anderson_step_ratio>=1.) ||
       !std::isfinite(o.max_anderson_step_ratio))
        throw std::invalid_argument("V51: invalid globalized HOLO controls");
    const auto start=std::chrono::steady_clock::now();
    DomV51Result result;
    std::vector<double> accepted=initial;
    auto s=evaluate_dom_variable_v50(ho,accepted);
    result.initial_true=s.true_relative;
    result.final_true=s.true_relative;
    result.accepted_defect_history.push_back(s.true_relative);
    const auto qualifies=[&](const DomV50Snapshot& a) {
        return a.fixed_relative<=opt.physical_tolerance &&
               a.transport_relative<=opt.physical_tolerance &&
               a.physical_cell_relative<=opt.physical_tolerance &&
               a.global_energy_relative<=opt.physical_tolerance &&
               a.maximum_trace_error<=opt.physical_tolerance;
    };
    auto publish=[&]() {
        result.converged=true;
        result.final_moments=accepted;
        result.diagnostic="all independent gas, gray-wall and transport gates accepted";
    };
    if(qualifies(s)) publish();
    std::deque<dom_v51_detail::Secant> history;
    for(std::size_t outer=0;outer<opt.max_outer && !result.converged;++outer) {
        std::vector<double> f(accepted.size());
        for(std::size_t k=0;k<f.size();++k)f[k]=s.projected[k]-accepted[k];
        history.push_back({accepted,f});
        while(history.size()>static_cast<std::size_t>(o.anderson_depth)+1)
            history.pop_front();
        bool committed=false;
        const auto attempt=[&](const std::vector<double>& step,int kind)->bool {
            if(step.size()!=accepted.size())return false;
            for(unsigned t=0;t<opt.max_backtracks;++t) {
                const double alpha=std::ldexp(1.,-static_cast<int>(t));
                auto candidate=accepted;
                bool positive=true;
                for(std::size_t k=0;k<candidate.size();++k) {
                    candidate[k]+=alpha*step[k];
                    if(!(candidate[k]>=0.) || !std::isfinite(candidate[k]))positive=false;
                }
                if(!positive){++result.rejected_trials;continue;}
                try {
                    auto next=evaluate_dom_variable_v50(ho,candidate);
                    if(next.true_relative<opt.required_contraction*s.true_relative ||
                       qualifies(next)) {
                        accepted=std::move(candidate);
                        s=std::move(next);
                        ++result.accepted_steps;
                        if(kind==0)++result.lo_steps;
                        else if(kind==1)++result.anderson_steps;
                        else ++result.picard_steps;
                        result.final_true=s.true_relative;
                        result.accepted_defect_history.push_back(s.true_relative);
                        if(qualifies(s))publish();
                        return true;
                    }
                } catch(const std::runtime_error&) {
                    // The HO angular solve or candidate produced an invalid state.
                }
                ++result.rejected_trials;
            }
            return false;
        };
        // Try an extrapolation of the full HO fixed-point including reflected
        // face moments BEFORE the frozen low-order update. The correction can
        // be rejected without modifying any accepted moment or energy field.
        if(o.use_anderson && o.anderson_depth>0) {
            const auto step=dom_v51_detail::direction(
                history,o.anderson_ridge,o.max_anderson_step_ratio);
            if(!step.empty())committed=attempt(step,1);
            else if(history.size()>1)++result.acceleration_restarts;
        }
        if(committed)continue;
        const auto nc=ho.n_cells();
        if(o.use_low_order) {
            std::vector<double> kt(nc),guess(nc);
            for(std::size_t c=0;c<nc;++c) {
                guess[c]=accepted[c];
                kt[c]=ho.cell_absorption()[c]+ho.cell_scattering()[c];
            }
            DomSmmGrayWallV49 lo(ho.source_mesh(),ho.source_geometry(),
                ho.directions(),ho.source_walls(),s.eddington,kt,ho.cell_absorption());
            const auto consistent=lo.consistent_source(guess,s.gas_ib,s.face_current);
            if(consistent.max_discrete_identity_relative>opt.physical_tolerance) {
                result.diagnostic="V51: refreshed HO/LO face identity rejected";
                break;
            }
            ++result.lo_solves;
            const auto low=solve_corrected_smm_v50(lo,consistent.corrected_rhs_watts,
                guess,opt.preconditioner,opt.linear_tolerance,opt.max_lo_krylov);
            if(low.converged) {
                std::vector<double> step(accepted.size());
                for(std::size_t c=0;c<nc;++c)step[c]=low.J[c]-accepted[c];
                for(std::size_t k=nc;k<step.size();++k)step[k]=f[k];
                if(o.use_wall_response && ho.n_reflecting_faces()>0) {
                    // The variation of wall R due to a changed volume J is
                    // part of the SAME angular transport operator. Recompute
                    // R with the LO candidate J and the current accepted R.
                    auto probe=accepted;
                    for(std::size_t c=0;c<nc;++c)probe[c]=low.J[c];
                    const auto response=ho.project(probe,true);
                    ++result.wall_response_maps;
                    for(std::size_t k=nc;k<step.size();++k)
                        step[k]=response[k]-accepted[k];
                }
                committed=attempt(step,0);
            } else ++result.lo_failures;
        }
        if(committed)continue;
        if(!opt.allow_picard_fallback) {
            result.diagnostic="V51: both accelerated directions rejected; rolled back";
            break;
        }
        committed=attempt(f,2);
        if(!committed) {
            result.diagnostic="V51: no contracting full-HO step; rolled back";
            break;
        }
    }
    if(!result.converged && result.diagnostic.empty())
        result.diagnostic="V51: iteration budget exceeded; no state published";
    result.transport_maps=ho.transport_maps();
    result.seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    return result;
}
} // namespace cfdx::physics
