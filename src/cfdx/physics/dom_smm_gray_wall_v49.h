#pragma once
// CFDX V49 - opt-in, research-only second-moment gray-wall low-order system.
// The gray Robin coefficient is derived from the ACTUAL discrete half-space
// angular moments and the owner-to-face resistance. This is an autonomous
// P1 half-range wall APPROXIMATION, not the exact high-order gray wall Schur.
// E and optical properties are frozen for each matrix-free Krylov solve.
// Rebuild after EVERY high-order update: no silent stale HO/LO closure.
// Production radiation and historical tolerances are unchanged.
#include "cfdx/physics/dom_smm_second_moment_v48.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct DomSmmGrayWallV49Result {
    bool accepted=false;
    bool amg_setup_ok=false;
    std::string diagnostic;
    std::size_t iterations=0;
    std::size_t preconditioner_levels=0;
    double true_linear_relative=std::numeric_limits<double>::infinity();
    double max_cell_balance_relative=std::numeric_limits<double>::infinity();
    double global_energy_relative=std::numeric_limits<double>::infinity();
    double maximum_quadrature_pi_error=0.;
    double seconds=0.;
    std::vector<double> accepted_J;
    std::vector<double> accepted_wall_flux_watts;
    std::vector<double> accepted_wall_R;
};

enum class DomSmmV49Preconditioner {None, JacobiAuxiliary, AmgAuxiliary};

struct DomSmmHoloSourceV49 {
    std::vector<double> corrected_rhs_watts;
    std::vector<double> face_ho_minus_lo_watts;
    double max_discrete_identity_relative=0.;
    double max_ho_energy_residual_relative=0.;
    // Diagnostics only. Computing the corrected source is NOT a successful
    // HOLO solve and MUST NOT be used as a convergence certificate.
};

class DomSmmGrayWallV49 {
public:
    using Tensor=DomSmmSecondMomentV48::Tensor;
    struct Wall {
        std::size_t face=0,owner=0;
        double conductance=0.,ib=0.,hemi=0.,moment2=0.,
               distance=0.,D=0.,emissivity=0.,area=0.;
    };

    DomSmmGrayWallV49(const cfdx::core::Mesh& mesh,const FvGeometry& geo,
                       const std::vector<DiscreteDirection>& directions,
                       const DomWallBoundaryConditions& bcs,
                       std::vector<Tensor> E,std::vector<double> extinction,
                       std::vector<double> absorption,
                       double maximum_hemisphere_error=0.03,
                       const std::vector<double>& wall_face_blackbody={})
       :mesh_(mesh),geo_(geo),smm_(mesh,geo),E_(std::move(E)),
        extinction_(std::move(extinction)),absorption_(std::move(absorption)) {
        const auto nc=mesh.n_cells(),nf=mesh.n_faces();
        if(E_.size()!=nc || extinction_.size()!=nc || absorption_.size()!=nc ||
           !std::isfinite(maximum_hemisphere_error) ||
           maximum_hemisphere_error<0. || maximum_hemisphere_error>=1. ||
           (!wall_face_blackbody.empty() && wall_face_blackbody.size()!=nf))
            throw std::invalid_argument("V49 SMM: dimensions / quadrature gate");
        validate_discrete_directions(directions);
        validate_dom_wall_conditions(mesh,bcs);
        for(std::size_t c=0;c<nc;++c) {
            const auto& e=E_[c];
            if(!(extinction_[c]>0.) || !std::isfinite(extinction_[c]) ||
               !(absorption_[c]>0.) || !std::isfinite(absorption_[c]) ||
               absorption_[c]>extinction_[c] ||
               !v37::tensor_is_spd({e[0],e[1],e[2],e[3],e[4],e[5]}) ||
               std::abs(e[0]+e[1]+e[2]-1.)>1.e-9)
                throw std::invalid_argument("V49 SMM: nonphysical optical coefficients or E");
        }
        std::vector<unsigned char> seen(nf,0);
        for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
            const auto& patch=mesh.boundary().patch(p);
            const auto& bc=bcs.at(patch.name);
            for(const auto f:patch.face_ids) {
                if(f>=nf || seen[f]++ || mesh.ownership().neighbour(f)>=0)
                    throw std::invalid_argument("V49 SMM: duplicate/inverted wall face");
                const auto o=mesh.ownership().owner(f);
                const auto sf=geo.face_area_vectors[f];
                const double area=sf.mag();
                const auto n=sf*(1./area);
                const double distance=(geo.face_centres[f]-geo.cell_centres[o]).dot(n);
                if(!(distance>0.) || !std::isfinite(distance))
                    throw std::invalid_argument("V49 SMM: invalid wall distance");
                const auto hm=dom_hemisphere_moments(directions,{n.x,n.y,n.z});
                max_pi_error_=std::max(max_pi_error_,hm.relative_pi_error);
                if(hm.relative_pi_error>maximum_hemisphere_error)
                    throw std::invalid_argument("V49 SMM: wall hemisphere quadrature not qualified");
                long double m2=0.;
                for(const auto& d:directions) {
                    const double mu=d.dx*n.x+d.dy*n.y+d.dz*n.z;
                    if(mu>0.) m2+=static_cast<long double>(d.weight)*mu*mu;
                }
                const double B=3.*static_cast<double>(m2)/(4.*std::acos(-1.));
                if(!(B>=0.) || !std::isfinite(B) || !(B<1.))
                    throw std::invalid_argument("V49 SMM: invalid half-range second moment");
                const auto& e=E_[o];
                const double Enn=e[0]*n.x*n.x+e[1]*n.y*n.y+e[2]*n.z*n.z+
                    2.*(e[3]*n.x*n.y+e[4]*n.x*n.z+e[5]*n.y*n.z);
                const double D=4.*std::acos(-1.)*Enn/extinction_[o];
                if(!(D>0.) || !std::isfinite(D))
                    throw std::invalid_argument("V49 SMM: degenerate wall second moment");
                // q_n = eps*H*(R-Ib), R=J_w+(B/H)*q_n and
                // J_w=J_owner-q_n*d/D; exactly eliminating q_n and R gives:
                // q_n = [eps*H/(1+eps*H*d/D-eps*B)]*(J_owner-Ib).
                const double denom=1.+bc.emissivity*hm.outgoing*distance/D-
                                   bc.emissivity*B;
                if(!(denom>0.) || !std::isfinite(denom))
                    throw std::invalid_argument("V49 SMM: nonpositive gray-wall Robin denominator");
                const double C=area*bc.emissivity*hm.outgoing/denom;
                if(!std::isfinite(C) || C<0.)
                    throw std::invalid_argument("V49 SMM: invalid gray-wall conductance");
                const double wall_ib=wall_face_blackbody.empty() ?
                    blackbody_intensity(bc.temperature) : wall_face_blackbody[f];
                if(!std::isfinite(wall_ib) || wall_ib<0.)
                    throw std::invalid_argument("V49 SMM: invalid face-resolved wall emission");
                walls_.push_back({f,o,C,wall_ib,
                                  hm.outgoing,static_cast<double>(m2),distance,D,
                                  bc.emissivity,area});
            }
        }
        for(std::size_t f=0;f<nf;++f)
            if(mesh.ownership().neighbour(f)<0 && !seen[f])
                throw std::invalid_argument("V49 SMM: incomplete radiative boundary partition");
        if(walls_.empty())throw std::invalid_argument("V49 SMM: empty walls");
    }

    const std::vector<Wall>& walls() const noexcept {return walls_;}
    double max_pi_error() const noexcept {return max_pi_error_;}
    std::size_t dimension() const noexcept {return mesh_.n_cells();}

    // V52 audit guard: a reconstructed physical ledger must use precisely the
    // geometry and frozen E/kappa arrays of the operator solved by GMRES.
    // Bitwise-equivalent numeric vectors and reference identity are required;
    // otherwise two different operators could yield a false certificate.
    bool frozen_operator_matches_v52(const cfdx::core::Mesh& mesh,
            const FvGeometry& geometry,const std::vector<Tensor>& E,
            const std::vector<double>& extinction,
            const std::vector<double>& absorption) const noexcept {
        return &mesh==&mesh_ && &geometry==&geo_ &&
               E==E_ && extinction==extinction_ && absorption==absorption_;
    }

    // True frozen-E matrix-free action, including the Robin wall derivative.
    // The V48 WLS tensor operator is not guaranteed SPD on every mesh, so use
    // GMRES rather than assuming CG symmetry/positivity.
    std::vector<double> apply(const std::vector<double>& direction) const {
        auto y=smm_.elliptic_action(direction,E_,extinction_,absorption_);
        for(const auto& w:walls_)y[w.owner]+=w.conductance*direction[w.owner];
        for(double v:y) if(!std::isfinite(v))
            throw std::runtime_error("V49 SMM: nonfinite matrix-free action");
        return y;
    }

    std::vector<double> rhs(const std::vector<double>& gas_ib,
                            const std::vector<double>& volumetric_source_watts={}) const {
        const auto nc=mesh_.n_cells();
        if(gas_ib.size()!=nc ||
           (!volumetric_source_watts.empty() && volumetric_source_watts.size()!=nc))
            throw std::invalid_argument("V49 SMM: source dimension mismatch");
        const double fourpi=4.*std::acos(-1.);
        std::vector<double> b(nc,0.);
        for(std::size_t c=0;c<nc;++c) {
            if(!std::isfinite(gas_ib[c]) || gas_ib[c]<0. ||
               (!volumetric_source_watts.empty() &&
                !std::isfinite(volumetric_source_watts[c])))
                throw std::invalid_argument("V49 SMM: invalid gas emission/source");
            b[c]=fourpi*absorption_[c]*geo_.cell_volumes[c]*gas_ib[c]+
                 (volumetric_source_watts.empty()?0.:volumetric_source_watts[c]);
        }
        for(const auto& w:walls_)b[w.owner]+=w.conductance*w.ib;
        return b;
    }

    // HO-to-LO discrete residual correction for a FRESH angular transport
    // snapshot. This is a face-conservative EXACT identity, not a frozen
    // closure that may be reused across nonlinear/HO iterations. It includes
    // gray walls (even purely reflecting ones) through their actual HO flux.
    // A physically unconverged HO snapshot remains unconverged after this call.
    DomSmmHoloSourceV49 consistent_source(const std::vector<double>& ho_J,
            const std::vector<double>& gas_ib,
            const std::vector<double>& ho_face_current_watts) const {
        if(ho_J.size()!=dimension() || ho_face_current_watts.size()!=mesh_.n_faces())
            throw std::invalid_argument("V49 HO/LO: incorrect high-order dimensions");
        std::vector<double> low_wall(mesh_.n_faces(),0.);
        for(const auto& w:walls_)
            low_wall[w.face]=w.conductance*(ho_J[w.owner]-w.ib);
        const auto lo=smm_.evaluate(ho_J,E_,extinction_,low_wall);
        const auto correction=smm_.reconcile(lo,ho_face_current_watts);
        DomSmmHoloSourceV49 result;
        result.face_ho_minus_lo_watts=correction.face_ho_correction;
        result.corrected_rhs_watts=rhs(gas_ib);
        std::vector<double> correction_divergence(dimension(),0.);
        for(std::size_t f=0;f<mesh_.n_faces();++f) {
            const double dq=result.face_ho_minus_lo_watts[f];
            const auto o=mesh_.ownership().owner(f);
            correction_divergence[o]+=dq;
            const auto n=mesh_.ownership().neighbour(f);
            if(n>=0)correction_divergence[static_cast<std::size_t>(n)]-=dq;
        }
        const auto a=apply(ho_J);
        const double fourpi=4.*std::acos(-1.);
        for(std::size_t c=0;c<dimension();++c) {
            result.corrected_rhs_watts[c]-=correction_divergence[c];
            const double independent=correction.reconstructed_divergence[c]+
                fourpi*absorption_[c]*geo_.cell_volumes[c]*(ho_J[c]-gas_ib[c]);
            const double reconstructed=a[c]-result.corrected_rhs_watts[c];
            const double scale=std::max({1.,std::abs(independent),
                std::abs(correction.reconstructed_divergence[c]),
                std::abs(a[c]),std::abs(result.corrected_rhs_watts[c])});
            result.max_discrete_identity_relative=std::max(
                result.max_discrete_identity_relative,
                std::abs(reconstructed-independent)/scale);
            result.max_ho_energy_residual_relative=std::max(
                result.max_ho_energy_residual_relative,std::abs(independent)/
                std::max({1.,std::abs(correction.reconstructed_divergence[c]),
                    std::abs(fourpi*absorption_[c]*geo_.cell_volumes[c]*
                             (ho_J[c]-gas_ib[c]))}));
        }
        return result;
    }

    // Sparse SPD AUXILIARY graph for preconditioning only. This 2-point graph
    // is NOT the V48 tensor WLS operator, and never replaces its true matvec.
    cfdx::core::SparseMatrix auxiliary_matrix() const {
        const auto nc=mesh_.n_cells();
        cfdx::core::SparseMatrix A(nc,nc);
        const double fourpi=4.*std::acos(-1.);
        for(std::size_t c=0;c<nc;++c)
            A.push_back(c,c,fourpi*absorption_[c]*geo_.cell_volumes[c]);
        for(std::size_t f=0;f<mesh_.n_faces();++f) {
            const auto owner=mesh_.ownership().owner(f);
            const auto nei=mesh_.ownership().neighbour(f);
            if(nei<0)continue;
            const auto j=static_cast<std::size_t>(nei);
            const auto sf=geo_.face_area_vectors[f];
            const double area=sf.mag();
            const auto n=sf*(1./area);
            const double d=(geo_.cell_centres[j]-geo_.cell_centres[owner]).dot(n);
            if(!(d>0.))throw std::invalid_argument("V49 SMM: inverted auxiliary face");
            const auto& a=E_[owner];const auto& b=E_[j];
            const Tensor e={(a[0]+b[0])*.5,(a[1]+b[1])*.5,(a[2]+b[2])*.5,
                            (a[3]+b[3])*.5,(a[4]+b[4])*.5,(a[5]+b[5])*.5};
            const double Enn=e[0]*n.x*n.x+e[1]*n.y*n.y+e[2]*n.z*n.z+
                              2.*(e[3]*n.x*n.y+e[4]*n.x*n.z+e[5]*n.y*n.z);
            const double extinction=2./(1./extinction_[owner]+1./extinction_[j]);
            const double C=fourpi*area*Enn/(extinction*d);
            if(!(C>0.) || !std::isfinite(C))
                throw std::invalid_argument("V49 SMM: invalid auxiliary face strength");
            A.push_back(owner,owner,C);A.push_back(j,j,C);
            A.push_back(owner,j,-C);A.push_back(j,owner,-C);
        }
        for(const auto& w:walls_)A.push_back(w.owner,w.owner,w.conductance);
        A.finalize();
        return A;
    }

    DomSmmGrayWallV49Result solve(const std::vector<double>& gas_ib,
            DomSmmV49Preconditioner policy=DomSmmV49Preconditioner::None,
            double tolerance=1.e-10,std::size_t max_iterations=250,
            int restart=30) const {
        if(!(tolerance>0.) || !std::isfinite(tolerance) ||
           max_iterations==0 || restart<2)
            throw std::invalid_argument("V49 SMM: invalid Krylov controls");
        const auto clock_start=std::chrono::steady_clock::now();
        DomSmmGrayWallV49Result result;
        result.maximum_quadrature_pi_error=max_pi_error_;
        const auto bvalues=rhs(gas_ib);
        const auto nc=dimension();
        cfdx::core::Vector b(nc,0.),x(nc,0.);
        for(std::size_t c=0;c<nc;++c){b(c)=bvalues[c];x(c)=gas_ib[c];}
        cfdx::core::LinearOperator op;
        op.size=nc;
        op.apply=[this,nc](const cfdx::core::Vector& v,cfdx::core::Vector& y){
            std::vector<double> direction(nc);
            for(std::size_t c=0;c<nc;++c)direction[c]=v(c);
            const auto out=apply(direction);
            if(y.size()!=nc)y.resize(nc);
            for(std::size_t c=0;c<nc;++c)y(c)=out[c];
        };
        // AMG setup on an explicitly documented SPD auxiliary graph; fine
        // apply uses the REAL tensor WLS operator. Reject setup failure.
        std::unique_ptr<cfdx::core::FunctionalLinearOperator> actual_operator;
        std::unique_ptr<cfdx::core::MatrixFreeVcyclePreconditioner> amg;
        cfdx::core::JacobiPreconditioner jacobi;
        cfdx::core::Preconditioner* pre=nullptr;
        if(policy!=DomSmmV49Preconditioner::None) {
            const auto aux=auxiliary_matrix();
            if(policy==DomSmmV49Preconditioner::JacobiAuxiliary) {
                if(!jacobi.setup(aux)) {
                    result.diagnostic="V49 SMM: auxiliary Jacobi setup failed";
                    return result;
                }
                pre=&jacobi;
            } else if(policy==DomSmmV49Preconditioner::AmgAuxiliary) {
                actual_operator=std::make_unique<cfdx::core::FunctionalLinearOperator>(nc,op.apply);
                amg=std::make_unique<cfdx::core::MatrixFreeVcyclePreconditioner>(
                    *actual_operator,0.55,2,2,0.25,12,
                    cfdx::core::AMGInterpolationPolicy::SmoothedAggregation,
                    false,false,cfdx::core::AMGSmootherPolicy::DampedJacobi);
                if(!amg->setup(aux)) {
                    result.diagnostic="V49 SMM: AMG failed: "+amg->last_setup_error();
                    return result;
                }
                result.amg_setup_ok=true;
                result.preconditioner_levels=amg->hierarchy_level_sizes().size();
                pre=amg.get();
            } else {
                throw std::invalid_argument("V49 SMM: unknown preconditioner");
            }
        }
        const auto linear=cfdx::core::solve_gmres(op,b,x,restart,max_iterations,
                    std::min(tolerance*0.001,1.e-13),
                                                  pre);
        result.iterations=linear.iterations;
        if(linear.status!=cfdx::core::SolverStatus::CONVERGED) {
            result.diagnostic="V49 SMM: GMRES did not close true algebraic residual";
            return result;
        }
        std::vector<double> candidate(nc);
        for(std::size_t c=0;c<nc;++c)candidate[c]=x(c);
        for(const auto value:candidate)
            if(!std::isfinite(value)||value<0.) {
                result.diagnostic="V49 SMM: negative/nonfinite J";return result;
            }
        const auto Ax=apply(candidate);
        double lin_num=0.,lin_den=1.e-300;
        for(std::size_t c=0;c<nc;++c) {
            lin_num=std::max(lin_num,std::abs(Ax[c]-bvalues[c]));
            lin_den=std::max(lin_den,std::abs(bvalues[c]));
        }
        result.true_linear_relative=lin_num/lin_den;
        std::vector<double> wall_flux(mesh_.n_faces(),0.),wall_R(mesh_.n_faces(),0.);
        for(const auto& w:walls_) {
            const double q=w.conductance*(candidate[w.owner]-w.ib);
            wall_flux[w.face]=q;
            const double q_area=q/w.area;
            const double Jface=candidate[w.owner]-q_area*w.distance/w.D;
            const double B=3.*w.moment2/(4.*std::acos(-1.));
            wall_R[w.face]=Jface+B*q_area/w.hemi;
            const double gray_power=w.area*w.emissivity*w.hemi*(wall_R[w.face]-w.ib);
            if(std::abs(gray_power-q)>1.e-10*std::max({1.,std::abs(q),std::abs(gray_power)})) {
                result.diagnostic="V49 SMM: gray wall moment power mismatch";
                return result;
            }
        }
        const auto fv=smm_.evaluate(candidate,E_,extinction_,wall_flux);
        const double fourpi=4.*std::acos(-1.);

        double max_cell=0.;
        for(std::size_t c=0;c<nc;++c) {
            const double volume=fourpi*absorption_[c]*geo_.cell_volumes[c]*
                               (candidate[c]-gas_ib[c]);
            const double defect=fv.divergence_watts[c]+volume;
            const double scale=std::max({1.,std::abs(fv.divergence_watts[c]),
                                            std::abs(volume)});
            max_cell=std::max(max_cell,std::abs(defect)/scale);
        }
        result.max_cell_balance_relative=max_cell;
        // diagnostic: independent global boundary and absorption ledger.
        long double ext_power=0.,volume_power=0.,scale_power=1.;
        for(const auto& w:walls_) {ext_power+=wall_flux[w.face];scale_power+=std::abs(wall_flux[w.face]);}
        for(std::size_t c=0;c<nc;++c){
          const long double v=static_cast<long double>(fourpi)*absorption_[c]*geo_.cell_volumes[c]*(candidate[c]-gas_ib[c]);
          volume_power+=v;scale_power+=std::abs(v);
        }
        result.global_energy_relative=std::abs(static_cast<double>(ext_power+volume_power))/static_cast<double>(scale_power);
        // Global ledger is evaluated separately from cell-divergence summation.
        // Never conceal a failing cell balance by correcting the global total.
        if(!(result.true_linear_relative<=tolerance) ||
           !(result.max_cell_balance_relative<=tolerance) ||
           !(result.global_energy_relative<=tolerance)) {
            result.diagnostic="V49 SMM: independent physical balance rejected";
            return result;
        }
        result.accepted=true;
        result.diagnostic="accepted (research-only, not HO-consistent wall closure)";
        result.accepted_J=std::move(candidate);
        result.accepted_wall_flux_watts=std::move(wall_flux);
        result.accepted_wall_R=std::move(wall_R);
        result.seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-clock_start).count();
        return result;
    }

private:
    const cfdx::core::Mesh& mesh_;
    const FvGeometry& geo_;
    DomSmmSecondMomentV48 smm_;
    std::vector<Tensor> E_;
    std::vector<double> extinction_,absorption_;
    std::vector<Wall> walls_;
    double max_pi_error_=0.;
};
} // namespace cfdx::physics
