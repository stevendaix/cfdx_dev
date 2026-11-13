#pragma once
// CFDX V50: experimental cellwise HO transport and refreshed HO/LO outer loop.
// No production entry point is modified. The HO angular FV matrix is the
// authority; the auxiliary second-moment system is never a physical oracle.
// An iterate is published ONLY if the independently reconstructed HO currents,
// gray-wall moments, true Schur fixed point, and global energy all satisfy
// the unmodified caller-supplied physical tolerance.
#include "cfdx/physics/dom_smm_gray_wall_v49.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct DomV50Snapshot {
    std::vector<double> projected, face_current, cell_divergence, gas_ib;
    std::vector<DomSmmSecondMomentV48::Tensor> eddington;
    double fixed_relative=0., transport_relative=0., physical_cell_relative=0.;
    double global_energy_relative=0., maximum_trace_error=0.;
    double true_relative=std::numeric_limits<double>::infinity();
};

// Reconstructs face fluxes without using the angular assembly residual.
// This is an independent physical check: every interior face contributes a
// single oriented number to owner and its negative to neighbour. Inhomogeneous
// absorption/scattering is evaluated cell-by-cell in the SAME HO stencil.
inline DomV50Snapshot evaluate_dom_variable_v50(
        DomMomentSchurSystemV44& transport, const std::vector<double>& x) {
    const auto& mesh=transport.source_mesh();
    const auto& geo=transport.source_geometry();
    const auto& angles=transport.directions();
    const auto& walls=transport.wall_faces();
    const auto& ka=transport.cell_absorption();
    const auto& ks=transport.cell_scattering();
    const auto nc=mesh.n_cells(), nf=mesh.n_faces();
    if(x.size()!=transport.dimension() || ka.size()!=nc || ks.size()!=nc)
        throw std::invalid_argument("V50: HO state/optical field dimensions");
    for(double v:x) if(!std::isfinite(v) || v<0.)
        throw std::invalid_argument("V50: negative/nonfinite candidate moment");
    std::vector<DomMomentSchurSystemV44::CellField> intensities;
    DomV50Snapshot s;
    s.projected=transport.project(x,true,&intensities);
    s.face_current.assign(nf,0.);
    s.cell_divergence.assign(nc,0.);
    s.gas_ib.resize(nc);
    s.eddington.assign(nc,{});
    std::vector<const DomMomentSchurSystemV44::WallFace*> boundary(nf,nullptr);
    for(const auto& w:walls) {
        if(w.face>=nf || boundary[w.face])
            throw std::invalid_argument("V50: duplicate HO wall geometry");
        boundary[w.face]=&w;
    }
    const double fourpi=4.*std::acos(-1.);
    for(std::size_t m=0;m<angles.size();++m) {
        const auto& d=angles[m];
        const double dir[6]={d.dx*d.dx,d.dy*d.dy,d.dz*d.dz,
                             d.dx*d.dy,d.dx*d.dz,d.dy*d.dz};
        for(std::size_t c=0;c<nc;++c) {
            const double value=intensities[m](c);
            if(!(value>=0.) || !std::isfinite(value))
                throw std::runtime_error("V50: nonphysical angular intensity");
            for(std::size_t k=0;k<6;++k)
                s.eddington[c][k]+=d.weight*value*dir[k];
        }
        for(std::size_t f=0;f<nf;++f) {
            const auto& sf=geo.face_area_vectors[f];
            const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
            const auto owner=mesh.ownership().owner(f);
            const auto nb=mesh.ownership().neighbour(f);
            double upwind=0.;
            if(mu>=0.) upwind=intensities[m](owner);
            else if(nb>=0) upwind=intensities[m](static_cast<std::size_t>(nb));
            else {
                const auto* w=boundary[f];
                if(!w || !(w->hminus>0.) || !(w->hplus>0.))
                    throw std::invalid_argument("V50: incomplete gray-wall geometry");
                const double reflection=w->reflect_index<transport.n_reflecting_faces()
                    ? (1.-w->emissivity)*(w->hplus/w->hminus)*
                      x[nc+w->reflect_index] : 0.;
                upwind=w->emissivity*w->ib+reflection;
            }
            s.face_current[f]+=d.weight*mu*upwind;
        }
    }
    for(std::size_t f=0;f<nf;++f) {
        const auto o=mesh.ownership().owner(f);
        const auto n=mesh.ownership().neighbour(f);
        s.cell_divergence[o]+=s.face_current[f];
        if(n>=0) s.cell_divergence[static_cast<std::size_t>(n)]-=s.face_current[f];
    }
    long double boundary_total=0.,vol_total=0.,scale_total=1.;
    for(std::size_t f=0;f<nf;++f) if(mesh.ownership().neighbour(f)<0) {
        boundary_total+=s.face_current[f];scale_total+=std::abs(s.face_current[f]);
    }
    for(std::size_t c=0;c<nc;++c) {
        const double J=s.projected[c];
        if(!(J>0.) || !std::isfinite(J))
            throw std::runtime_error("V50: nonpositive projected HO scalar moment");
        for(auto& tensor:s.eddington[c]) tensor/=(fourpi*J);
        const auto& e=s.eddington[c];
        s.maximum_trace_error=std::max(s.maximum_trace_error,
            std::abs(e[0]+e[1]+e[2]-1.));
        s.gas_ib[c]=blackbody_intensity(transport.source_temperature()(c));
        const double volume=geo.cell_volumes[c];
        const double absorption=fourpi*volume*ka[c]*(J-s.gas_ib[c]);
        const double lag=fourpi*volume*ks[c]*(J-x[c]);
        const double physical=s.cell_divergence[c]+absorption;
        const double angular=physical+lag;
        const double scale=std::max({1.,std::abs(s.cell_divergence[c]),
            std::abs(absorption),std::abs(lag),fourpi*volume*(ka[c]+ks[c])*J});
        s.transport_relative=std::max(s.transport_relative,std::abs(angular)/scale);
        s.physical_cell_relative=std::max(s.physical_cell_relative,std::abs(physical)/scale);
        vol_total+=absorption;
        scale_total+=std::abs(absorption);
    }
    s.global_energy_relative=std::abs(static_cast<double>(boundary_total+vol_total)) /
                             static_cast<double>(scale_total);
    for(std::size_t i=0;i<x.size();++i)
        s.fixed_relative=std::max(s.fixed_relative,
            std::abs(x[i]-s.projected[i])/
            std::max({1.,std::abs(x[i]),std::abs(s.projected[i])}));
    s.true_relative=std::max({s.fixed_relative,s.transport_relative,
         s.physical_cell_relative,s.global_energy_relative,s.maximum_trace_error});
    return s;
}

struct DomV50LinearCorrection {
    bool converged=false;
    std::size_t iterations=0,amg_levels=0;
    double relative=std::numeric_limits<double>::infinity();
    std::string diagnostic;
    std::vector<double> J;
};

// Linear solve for a FROZEN HO-consistent LO right-hand side; its residual is
// algebraic only. Publication / physical acceptance belongs to the outer HO.
inline DomV50LinearCorrection solve_corrected_smm_v50(
    const DomSmmGrayWallV49& system,const std::vector<double>& rhs,
    const std::vector<double>& initial,DomSmmV49Preconditioner method,
    double linear_tol,std::size_t budget=180) {
    const auto n=system.dimension();
    if(rhs.size()!=n || initial.size()!=n || !(linear_tol>0.) ||
       !std::isfinite(linear_tol) || budget==0)
        throw std::invalid_argument("V50 LO: invalid corrected system");
    DomV50LinearCorrection result;
    cfdx::core::Vector b(n,0.),x(n,0.);
    for(std::size_t c=0;c<n;++c) {
        if(!std::isfinite(rhs[c]) || !std::isfinite(initial[c]))
            throw std::invalid_argument("V50 LO: nonfinite corrected source");
        b(c)=rhs[c];x(c)=initial[c];
    }
    cfdx::core::LinearOperator op;
    op.size=n;
    op.apply=[&system,n](const cfdx::core::Vector& in,cfdx::core::Vector& out) {
        std::vector<double> v(n);
        for(std::size_t i=0;i<n;++i)v[i]=in(i);
        const auto y=system.apply(v);
        if(out.size()!=n)out.resize(n);
        for(std::size_t i=0;i<n;++i)out(i)=y[i];
    };
    std::unique_ptr<cfdx::core::FunctionalLinearOperator> fine;
    std::unique_ptr<cfdx::core::MatrixFreeVcyclePreconditioner> amg;
    cfdx::core::JacobiPreconditioner jacobi;
    cfdx::core::Preconditioner* pre=nullptr;
    if(method!=DomSmmV49Preconditioner::None) {
        const auto aux=system.auxiliary_matrix();
        if(method==DomSmmV49Preconditioner::JacobiAuxiliary) {
            if(!jacobi.setup(aux)){result.diagnostic="Jacobi setup rejected";return result;}
            pre=&jacobi;
        } else if(method==DomSmmV49Preconditioner::AmgAuxiliary) {
            fine=std::make_unique<cfdx::core::FunctionalLinearOperator>(n,op.apply);
            amg=std::make_unique<cfdx::core::MatrixFreeVcyclePreconditioner>(
                *fine,.55,2,2,.25,12,cfdx::core::AMGInterpolationPolicy::SmoothedAggregation,
                false,false,cfdx::core::AMGSmootherPolicy::DampedJacobi);
            if(!amg->setup(aux)) {
                result.diagnostic="AMG setup rejected: "+amg->last_setup_error();return result;
            }
            result.amg_levels=amg->hierarchy_level_sizes().size();
            pre=amg.get();
        } else throw std::invalid_argument("V50 LO: unsupported preconditioner");
    }
    const auto solved=cfdx::core::solve_gmres(op,b,x,30,budget,linear_tol,pre);
    result.iterations=solved.iterations;
    if(solved.status!=cfdx::core::SolverStatus::CONVERGED) {
        result.diagnostic="corrected low-order GMRES did not converge";return result;
    }
    result.J.resize(n);
    for(std::size_t c=0;c<n;++c) result.J[c]=x(c);
    const auto actual=system.apply(result.J);
    double num=0.,den=1.;
    for(std::size_t c=0;c<n;++c) {
        num=std::max(num,std::abs(actual[c]-rhs[c]));
        den=std::max(den,std::abs(rhs[c]));
    }
    result.relative=num/den;
    result.converged=result.relative<=linear_tol;
    if(!result.converged) {
        result.J.clear();result.diagnostic="true corrected LO residual rejected";
    }
    return result;
}

struct DomV50Options {
    std::size_t max_outer=60, max_lo_krylov=180;
    unsigned max_backtracks=10;
    double physical_tolerance=1e-9, linear_tolerance=1e-12;
    double required_contraction=.999;
    DomSmmV49Preconditioner preconditioner=DomSmmV49Preconditioner::JacobiAuxiliary;
    bool allow_picard_fallback=true;
};
struct DomV50Result {
    bool converged=false;
    std::size_t accepted_steps=0,rejected_trials=0,lo_solves=0,lo_failures=0;
    std::size_t lo_accepted_steps=0,picard_fallbacks=0,transport_maps=0;
    double initial_true=std::numeric_limits<double>::infinity();
    double final_true=std::numeric_limits<double>::infinity();
    double seconds=0.;
    std::string diagnostic;
    std::vector<double> final_moments;
    std::vector<double> accepted_defect_history; // diagnostics, not physical states
};

// Transactional HO/LO outer loop. E, HO wall R and HO currents are refreshed
// on every accepted state and EVERY line-search trial. A rejected trial can
// never replace the last accepted state. A converged state is published only
// after a fresh full angular transport and all physical checks.
inline DomV50Result solve_dom_holo_v50(DomMomentSchurSystemV44& ho,
                                      const std::vector<double>& initial,
                                      const DomV50Options& options={}) {
    if(!(options.physical_tolerance>0.) || !std::isfinite(options.physical_tolerance) ||
       !(options.linear_tolerance>0.) || !std::isfinite(options.linear_tolerance) ||
       !(options.required_contraction>0.) || !(options.required_contraction<1.) ||
       options.max_outer==0 || options.max_lo_krylov==0 || options.max_backtracks==0 ||
       initial.size()!=ho.dimension())
        throw std::invalid_argument("V50 HOLO: invalid controls or state dimensions");
    const auto t0=std::chrono::steady_clock::now();
    DomV50Result result;
    std::vector<double> accepted=initial;
    auto s=evaluate_dom_variable_v50(ho,accepted);
    result.initial_true=s.true_relative;
    result.final_true=s.true_relative;
    result.accepted_defect_history.push_back(s.true_relative);
    const auto qualifies=[&](const DomV50Snapshot& a) {
        return a.fixed_relative<=options.physical_tolerance &&
               a.transport_relative<=options.physical_tolerance &&
               a.physical_cell_relative<=options.physical_tolerance &&
               a.global_energy_relative<=options.physical_tolerance &&
               a.maximum_trace_error<=options.physical_tolerance;
    };
    auto publish=[&]() {
        result.converged=true;
        result.final_moments=accepted;
        result.diagnostic="all independently reconstructed HO balances accepted";
    };
    if(qualifies(s))publish();
    for(std::size_t outer=0;outer<options.max_outer && !result.converged;++outer) {
        const auto nc=ho.n_cells();
        std::vector<double> guess(nc),extinction(nc);
        for(std::size_t c=0;c<nc;++c) {
            guess[c]=accepted[c];
            extinction[c]=ho.cell_absorption()[c]+ho.cell_scattering()[c];
        }
        DomSmmGrayWallV49 lo(ho.source_mesh(),ho.source_geometry(),ho.directions(),
                             ho.source_walls(),s.eddington,extinction,
                             ho.cell_absorption());
        const auto source=lo.consistent_source(guess,s.gas_ib,s.face_current);
        if(source.max_discrete_identity_relative>options.physical_tolerance) {
            result.diagnostic="refreshed HO/LO face identity rejected";break;
        }
        ++result.lo_solves;
        const auto correction=solve_corrected_smm_v50(lo,source.corrected_rhs_watts,
                       guess,options.preconditioner,options.linear_tolerance,
                       options.max_lo_krylov);
        std::vector<double> step(accepted.size(),0.);
        bool use_lo=correction.converged;
        if(use_lo) {
            for(std::size_t c=0;c<nc;++c)step[c]=correction.J[c]-accepted[c];
            for(std::size_t k=nc;k<step.size();++k)
                step[k]=s.projected[k]-accepted[k];
        } else {
            ++result.lo_failures;
            if(!options.allow_picard_fallback) {
                result.diagnostic="corrected HO/LO Krylov failed; no state published";break;
            }
        }
        bool step_accepted=false;
        // Trial states always reassemble the HO transport, including gray R.
        for(unsigned branch=0;branch<(use_lo && options.allow_picard_fallback?2u:1u) &&
                            !step_accepted;++branch) {
            if(!use_lo || branch==1) {
                ++result.picard_fallbacks;
                for(std::size_t k=0;k<step.size();++k)
                    step[k]=s.projected[k]-accepted[k];
            }
            for(unsigned trial=0;trial<options.max_backtracks;++trial) {
                const double alpha=std::ldexp(1.,-static_cast<int>(trial));
                std::vector<double> candidate=accepted;
                bool valid=true;
                for(std::size_t k=0;k<candidate.size();++k) {
                    candidate[k]+=alpha*step[k];
                    if(!(candidate[k]>=0.) || !std::isfinite(candidate[k]))valid=false;
                }
                if(!valid){++result.rejected_trials;continue;}
                try {
                    auto next=evaluate_dom_variable_v50(ho,candidate);
                    if(next.true_relative<options.required_contraction*s.true_relative ||
                       qualifies(next)) {
                        accepted=std::move(candidate);
                        s=std::move(next);
                        ++result.accepted_steps;
                        if(use_lo && branch==0)++result.lo_accepted_steps;
                        result.final_true=s.true_relative;
                        result.accepted_defect_history.push_back(s.true_relative);
                        if(qualifies(s))publish();
                        step_accepted=true;
                        break;
                    }
                } catch(const std::runtime_error&) {
                    // Signed/nonphysical angular intensities: reject candidate.
                }
                ++result.rejected_trials;
            }
        }
        if(!step_accepted) {
            result.diagnostic="no HO-consistent contracting step; state rolled back";
            break;
        }
    }
    if(!result.converged && result.diagnostic.empty())
        result.diagnostic="outer HO/LO budget exceeded; no state published";
    result.transport_maps=ho.transport_maps();
    result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    return result;
}
} // namespace cfdx::physics
