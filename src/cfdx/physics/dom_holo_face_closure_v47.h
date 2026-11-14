#pragma once
// V47 research-only, source-consistent HOLO closure for the discrete UPWIND
// DOM transport. This is not a new production radiation solver. At each state
// x=(cell J, gray-wall R), recompute angular intensities and oriented face
// currents; a frozen-drift low-order diffusion step is accepted ONLY after
// the full high-order fixed-point defect has contracted. In particular a
// small low-order algebraic residual cannot certify the radiative energy.
#include "cfdx/physics/dom_face_moment_v45.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct DomHoloSnapshotV47 {
    std::vector<double> moments;
    std::vector<double> projected;
    std::vector<double> face_current;
    std::vector<double> face_diffusion;
    std::vector<double> face_closure;
    std::vector<double> low_order_energy_watts; // one signed residual per cell
    std::vector<double> schur_defect; // (J-J_HO, R-R_HO)
    // Symmetric Eddington tensor E = sum_m(w I_m s_m s_m)/(4*pi J_HO).
    // 6 per cell: xx, yy, zz, xy, xz, yz. It is a diagnostic HO
    // second moment, not a separately discretized elliptic operator.
    std::vector<double> eddington_symmetric;
    double max_eddington_trace_error=0.;
    DomFaceMomentAuditV45 independent;
    double max_discrete_identity_relative=0.;
    double max_fixed_point_relative=0.;
};
struct DomHoloRunV47 {
    bool converged=false;
    std::size_t accepted_steps=0;
    std::size_t rejected_steps=0;
    double initial_defect=std::numeric_limits<double>::infinity();
    double final_defect=std::numeric_limits<double>::infinity();
    double physical_energy=std::numeric_limits<double>::infinity();
    double seconds=0.;
    // Fields are exported only after physical acceptance.
    std::vector<double> accepted_moments;
};
struct DomHoloStepV47 {
    bool contracted=false;
    bool qualified=false;
    unsigned line_search_trials=0;
    double before=std::numeric_limits<double>::infinity();
    double after=std::numeric_limits<double>::infinity();
    // No candidate is published on failure.
    std::vector<double> accepted_moments;
};

class DomHoloFaceClosureV47 {
public:
    using CellField=DomMomentSchurSystemV44::CellField;
    DomHoloFaceClosureV47(DomMomentSchurSystemV44& transport,
                         const cfdx::core::Mesh& mesh,
                         const FvGeometry& geo,
                         const CellField& temperature,
                         const std::vector<DiscreteDirection>& directions,
                         const DomWallBoundaryConditions& walls,
                         const RadiationTransportControls& controls)
        : transport_(transport),mesh_(mesh),geo_(geo),temperature_(temperature),
          directions_(directions),walls_(walls),controls_(controls) {
        if(&mesh!=&transport.source_mesh() ||
           &geo!=&transport.source_geometry() ||
           &temperature!=&transport.source_temperature() ||
           &walls!=&transport.source_walls() ||
           mesh.n_cells()!=transport.n_cells() || geo_.cell_volumes.size()!=mesh.n_cells() ||
           temperature_.size()!=mesh.n_cells() ||
           !(controls_.absorption>0.) || !(controls_.scattering>=0.) ||
           !std::isfinite(controls_.absorption) || !std::isfinite(controls_.scattering))
            throw std::invalid_argument("V47 HOLO: invalid geometry / coefficients");
        const auto& tc=transport_.transport_controls();
        if(tc.absorption!=controls.absorption || tc.scattering!=controls.scattering ||
           transport_.directions().size()!=directions.size())
            throw std::invalid_argument("V47 HOLO: optical operator mismatch");
        for(std::size_t m=0;m<directions.size();++m) {
            const auto& a=transport_.directions()[m];const auto& b=directions[m];
            if(a.dx!=b.dx || a.dy!=b.dy || a.dz!=b.dz || a.weight!=b.weight)
                throw std::invalid_argument("V47 HOLO: angular operator mismatch");
        }
        // Validate the same wall geometry as the authoritative Schur system.
        const auto& faces=transport_.wall_faces();
        wall_at_face_.assign(mesh.n_faces(),nullptr);
        for(const auto& w:faces) {
            if(w.face>=wall_at_face_.size() || wall_at_face_[w.face])
                throw std::invalid_argument("V47 HOLO: wall partition mismatch");
            wall_at_face_[w.face]=&w;
        }
        for(std::size_t f=0;f<mesh.n_faces();++f)
            if(mesh.ownership().neighbour(f)<0 && !wall_at_face_[f])
                throw std::invalid_argument("V47 HOLO: missing boundary face");
    }

    // One fresh HO transport solve per snapshot. It is deliberately NOT a
    // static correction carried over from a previous Picard state.
    DomHoloSnapshotV47 evaluate(const std::vector<double>& x) {
        validate_moments(x);
        std::vector<CellField> intensity;
        DomHoloSnapshotV47 s;
        s.moments=x;
        s.projected=transport_.project(x,true,&intensity);
        s.independent=audit_dom_face_moments_v45(
            mesh_,geo_,temperature_,directions_,walls_,controls_,transport_,x,intensity);
        s.face_current=s.independent.high_order_face_current;
        s.face_diffusion=s.independent.diffusive_face_current;
        s.face_closure=s.independent.closure_face_current;
        s.eddington_symmetric.assign(6*mesh_.n_cells(),0.);
        for(std::size_t m=0;m<directions_.size();++m) {
            const auto& d=directions_[m];
            const double terms[6]={d.dx*d.dx,d.dy*d.dy,d.dz*d.dz,
                                   d.dx*d.dy,d.dx*d.dz,d.dy*d.dz};
            for(std::size_t c=0;c<mesh_.n_cells();++c)
                for(std::size_t k=0;k<6;++k)
                    s.eddington_symmetric[6*c+k]+=d.weight*intensity[m](c)*terms[k];
        }
        const double fourpi=4.*std::acos(-1.);
        for(std::size_t c=0;c<mesh_.n_cells();++c) {
            const double G=fourpi*s.projected[c];
            if(!(G>0.) || !std::isfinite(G))
                throw std::runtime_error("V47 HOLO: nonpositive high-order intensity moment");
            for(std::size_t k=0;k<6;++k)
                s.eddington_symmetric[6*c+k]/=G;
            const double trace=s.eddington_symmetric[6*c]+
                s.eddington_symmetric[6*c+1]+s.eddington_symmetric[6*c+2];
            s.max_eddington_trace_error=std::max(
                s.max_eddington_trace_error,std::abs(trace-1.));
        }
        s.low_order_energy_watts.assign(mesh_.n_cells(),0.);
        accumulate_divergence(s.face_current,s.low_order_energy_watts);
        s.schur_defect.resize(x.size());
        const double pi=std::acos(-1.);
        for(std::size_t c=0;c<mesh_.n_cells();++c) {
            const double V=geo_.cell_volumes[c];
            s.low_order_energy_watts[c]+=4.*pi*controls_.absorption*V*
                (x[c]-blackbody_intensity(temperature_(c)));
            s.schur_defect[c]=x[c]-s.projected[c];
            const double as_cell_res=s.low_order_energy_watts[c]/
                (4.*pi*(controls_.absorption+controls_.scattering)*V);
            const double scale=std::max({1.,std::abs(x[c]),std::abs(s.projected[c])});
            s.max_discrete_identity_relative=std::max(
                s.max_discrete_identity_relative,
                std::abs(as_cell_res-s.schur_defect[c])/scale);
        }
        for(std::size_t i=0;i<x.size();++i) {
            s.schur_defect[i]=x[i]-s.projected[i];
            s.max_fixed_point_relative=std::max(s.max_fixed_point_relative,
                std::abs(s.schur_defect[i])/
                std::max({1.,std::abs(x[i]),std::abs(s.projected[i])}));
        }
        return s;
    }

    // Linear, fully refreshed face-current derivative. The wall emission is
    // ZERO in this action; gray wall reflection depends on the varied R.
    // Therefore this is the actual Jacobian action of the affine discrete
    // transport, not a finite-difference approximation to frozen delta-Q.
    std::vector<double> linearized_low_order_defect(const std::vector<double>& v,
                                                     double* identity_error=nullptr) {
        validate_moments(v,false);
        std::vector<CellField> angular;
        const auto projected=transport_.project(v,false,&angular);
        const auto current=linear_face_current(v,angular);
        std::vector<double> div(mesh_.n_cells(),0.);
        accumulate_divergence(current,div);
        const double pi=std::acos(-1.);
        std::vector<double> out(v.size(),0.);
        double error=0.;
        for(std::size_t c=0;c<mesh_.n_cells();++c) {
            const double V=geo_.cell_volumes[c];
            out[c]=(div[c]+4.*pi*controls_.absorption*V*v[c])/
                (4.*pi*(controls_.absorption+controls_.scattering)*V);
            const double exact=v[c]-projected[c];
            error=std::max(error,std::abs(out[c]-exact)/
                std::max({1.,std::abs(v[c]),std::abs(projected[c])}));
        }
        for(std::size_t i=mesh_.n_cells();i<v.size();++i)
            out[i]=v[i]-projected[i];
        if(identity_error) *identity_error=error;
        return out;
    }

    // Solve one frozen-closure LO correction with the same metrics and
    // emissivities used by V44's conservative FV diffusion matrix. Gray R
    // variables receive their own Schur defect (not folded into gas sources).
    // A true HO reevaluation and backtracking are mandatory for acceptance.
    DomHoloStepV47 attempt_corrected_low_order_step(const std::vector<double>& x,
                        double required_contraction=0.995,unsigned max_trials=12,
                        const DomLocalAngularJacobiV45* local=nullptr,
                        double diffusion_fraction=0.) {
        return attempt_from_snapshot(x,evaluate(x),required_contraction,
                                     max_trials,local,diffusion_fraction);
    }

private:
    DomHoloStepV47 attempt_from_snapshot(const std::vector<double>& x,
                        const DomHoloSnapshotV47& initial,
                        double required_contraction,unsigned max_trials,
                        const DomLocalAngularJacobiV45* local,
                        double diffusion_fraction) {
        if(!(required_contraction>0.) || !(required_contraction<1.) ||
           max_trials==0 || max_trials>30 ||
           !std::isfinite(diffusion_fraction) || diffusion_fraction<0. ||
           diffusion_fraction>1. || (!local && diffusion_fraction!=0.))
            throw std::invalid_argument("V47 HOLO: invalid line-search parameters");
        DomHoloStepV47 result;
        result.before=initial.max_fixed_point_relative;
        if(!std::isfinite(result.before) || !(result.before>0.)) return result;
        std::vector<double> step(x.size(),0.);
        if(local) {
            const auto local_response=local->apply(initial.schur_defect);
            if(local_response.size()!=x.size())
                throw std::invalid_argument("V47 HOLO: local preconditioner size");
            for(std::size_t i=0;i<x.size();++i) step[i]=-local_response[i];
        }
        if(!local || diffusion_fraction>0.) {
            ScalarEquation diffusion;
            diffusion.matrix=transport_.diffusion_matrix();
            diffusion.rhs=cfdx::core::Vector(mesh_.n_cells(),0.);
            const double fourpi=4.*std::acos(-1.);
            for(std::size_t c=0;c<mesh_.n_cells();++c)
                diffusion.rhs(c)=-initial.low_order_energy_watts[c]/fourpi;
            cfdx::core::Vector delta(mesh_.n_cells(),0.);
            ScalarSolveControls linear;
            linear.tolerance=controls_.linear_tolerance;
            linear.max_iterations=controls_.linear_max_iterations;
            linear.relaxation=1.;
            const auto solved=solve_scalar_equation(diffusion,delta,linear);
            if(solved.status!=cfdx::core::SolverStatus::CONVERGED) return result;
            if(local) {
                for(std::size_t c=0;c<mesh_.n_cells();++c)
                    step[c]=(1.-diffusion_fraction)*step[c]+
                        diffusion_fraction*delta(c);
            } else {
                for(std::size_t c=0;c<mesh_.n_cells();++c) step[c]=delta(c);
                for(std::size_t i=mesh_.n_cells();i<x.size();++i)
                    step[i]=-initial.schur_defect[i];
            }
        }
        for(unsigned trial=0;trial<max_trials;++trial) {
            result.line_search_trials=trial+1;
            const double alpha=std::ldexp(1.,-static_cast<int>(trial));
            std::vector<double> candidate=x;
            bool positive=true;
            for(std::size_t i=0;i<x.size();++i) {
                candidate[i]+=alpha*step[i];
                if(!std::isfinite(candidate[i]) || candidate[i]<0.) positive=false;
            }
            if(!positive) continue;
            const auto after=evaluate(candidate);
            if(after.max_fixed_point_relative<required_contraction*result.before &&
               after.max_discrete_identity_relative<=controls_.tolerance) {
                result.contracted=true;
                result.after=after.max_fixed_point_relative;
                // Distinguish contraction from final physical qualification.
                result.qualified=dom_face_moment_qualified_v45(
                    after.independent,controls_.tolerance) &&
                    after.max_discrete_identity_relative<=controls_.tolerance &&
                    after.max_fixed_point_relative<=controls_.tolerance;
                result.accepted_moments=std::move(candidate);
                return result;
            }
        }
        return result;
    }

public:
    // Outer HOLO iteration: every low-order trial is followed by a new HO
    // angular solve and an independently recomputed face-energy balance.
    // Failure has no externally published field. This is a research method;
    // no claim of unconditional convergence or production performance.
    DomHoloRunV47 iterate_refreshed_closure(std::vector<double> x,
                                             std::size_t max_steps=100,
                                             const DomLocalAngularJacobiV45* local=nullptr,
                                             double diffusion_fraction=0.) {
        validate_moments(x);
        const auto start=std::chrono::steady_clock::now();
        if(max_steps==0 || max_steps>1000)
            throw std::invalid_argument("V47 HOLO: invalid outer budget");
        DomHoloRunV47 run;
        for(std::size_t k=0;k<=max_steps;++k) {
            const auto current=evaluate(x);
            if(k==0) run.initial_defect=current.max_fixed_point_relative;
            run.final_defect=current.max_fixed_point_relative;
            run.physical_energy=current.independent.global_energy_relative;
            if(current.max_fixed_point_relative<=controls_.tolerance &&
               current.max_discrete_identity_relative<=controls_.tolerance &&
               dom_face_moment_qualified_v45(current.independent,controls_.tolerance)) {
                run.converged=true;
                run.accepted_moments=std::move(x);
                run.seconds=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-start).count();
                return run;
            }
            if(k==max_steps) break;
            const auto next=attempt_from_snapshot(
                x,current,.995,12,local,diffusion_fraction);
            if(!next.contracted) {
                ++run.rejected_steps;
                break;
            }
            x=next.accepted_moments;
            ++run.accepted_steps;
        }
        run.seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
        return run;
    }

private:
    void validate_moments(const std::vector<double>& x,bool positive=true) const {
        if(x.size()!=transport_.dimension())
            throw std::invalid_argument("V47 HOLO: invalid moment dimension");
        for(const double value:x) if(!std::isfinite(value) || (positive && value<0.))
            throw std::invalid_argument("V47 HOLO: invalid moment data");
    }
    void accumulate_divergence(const std::vector<double>& q,
                               std::vector<double>& div) const {
        if(q.size()!=mesh_.n_faces() || div.size()!=mesh_.n_cells())
            throw std::invalid_argument("V47 HOLO: face divergence size");
        for(std::size_t f=0;f<q.size();++f) {
            const auto owner=mesh_.ownership().owner(f);
            const auto nei=mesh_.ownership().neighbour(f);
            if(owner>=div.size() || !std::isfinite(q[f]))
                throw std::invalid_argument("V47 HOLO: invalid current");
            div[owner]+=q[f];
            if(nei>=0) {
                if(static_cast<std::size_t>(nei)>=div.size())
                    throw std::invalid_argument("V47 HOLO: bad neighbor");
                div[static_cast<std::size_t>(nei)]-=q[f];
            }
        }
    }
    std::vector<double> linear_face_current(const std::vector<double>& moments,
                                           const std::vector<CellField>& angular) const {
        if(angular.size()!=directions_.size())
            throw std::invalid_argument("V47 HOLO: angular dimension mismatch");
        std::vector<double> current(mesh_.n_faces(),0.);
        for(std::size_t m=0;m<directions_.size();++m) {
            const auto& d=directions_[m];
            if(angular[m].size()!=mesh_.n_cells())
                throw std::invalid_argument("V47 HOLO: invalid angular cell count");
            for(std::size_t f=0;f<mesh_.n_faces();++f) {
                const auto& sf=geo_.face_area_vectors[f];
                const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                const auto own=mesh_.ownership().owner(f);
                const auto nei=mesh_.ownership().neighbour(f);
                double upwind=0.;
                if(mu>=0.) upwind=angular[m](own);
                else if(nei>=0) upwind=angular[m](static_cast<std::size_t>(nei));
                else {
                    const auto* w=wall_at_face_[f];
                    if(!w || !(w->hminus>0.))
                        throw std::invalid_argument("V47 HOLO: invalid gray wall");
                    upwind=w->reflect_index<transport_.n_reflecting_faces()
                        ? (1.-w->emissivity)*w->hplus/w->hminus*
                          moments[mesh_.n_cells()+w->reflect_index] : 0.;
                }
                if(!std::isfinite(upwind))
                    throw std::invalid_argument("V47 HOLO: nonfinite face intensity");
                current[f]+=d.weight*mu*upwind;
            }
        }
        return current;
    }
    DomMomentSchurSystemV44& transport_;
    const cfdx::core::Mesh& mesh_;
    const FvGeometry& geo_;
    const CellField& temperature_;
    const std::vector<DiscreteDirection>& directions_;
    const DomWallBoundaryConditions& walls_;
    RadiationTransportControls controls_;
    std::vector<const DomMomentSchurSystemV44::WallFace*> wall_at_face_;
};
} // namespace cfdx::physics
