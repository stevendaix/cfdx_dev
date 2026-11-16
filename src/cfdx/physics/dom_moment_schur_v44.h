#pragma once
// V44 opt-in research operator. The reduced linear system is the EXACT
// cell-mean / diffuse-wall-outgoing moment Schur complement of the existing
// upwind finite-volume angular matrices (at fixed T and optical properties).
// This is NOT a claim that the optional diffusion preconditioner itself is a
// discretely consistent DSA closure, nor a new default radiation solver.
//
// State: x = (J_c, R_f), J=G/(4*pi), R=G_out/Hplus. The R block is present
// only on reflecting boundary faces (emissivity < 1). Eliminating the angular
// intensities gives (I-P A^-1 S)x = P A^-1 b, with the same face areas,
// incoming ordinate values, and bounded upwind FV stencil as production DOM.
#include "cfdx/physics/dom_dsa_study_v43.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct DomMomentSchurV44Result {
    bool converged = false;
    std::size_t krylov_iterations = 0;
    std::size_t transport_maps = 0;
    std::size_t diffusion_solves = 0;
    std::size_t reflecting_faces = 0;
    // V46 research-only: additional true-residual restarts requested by a
    // strict independent physical gate, never by the default V44 path.
    std::size_t physical_refinements_v46 = 0;
    std::size_t physical_rejections_v46 = 0;
    double true_schur_relative = std::numeric_limits<double>::infinity();
    double local_transport_relative = std::numeric_limits<double>::infinity();
    double global_energy_relative = std::numeric_limits<double>::infinity();
    double seconds = 0.;
    std::vector<double> scalar_mean_intensity;
    std::vector<double> reflected_face_mean_intensity;
};

inline double dom_schur_dot_v44(const std::vector<double>& x,
                               const std::vector<double>& y) {
    if (x.size()!=y.size()) throw std::invalid_argument("DOM Schur V44: dot size");
    long double s=0.;
    for (std::size_t i=0;i<x.size();++i) s+=static_cast<long double>(x[i])*y[i];
    return static_cast<double>(s);
}
inline double dom_schur_norm_v44(const std::vector<double>& x) {
    return std::sqrt(std::max(0.,dom_schur_dot_v44(x,x)));
}
inline double dom_schur_max_relative_v44(const std::vector<double>& x,
                                        const std::vector<double>& reference) {
    if(x.size()!=reference.size())
        throw std::invalid_argument("DOM Schur V44: residual size");
    double numerator=0.,denominator=1.;
    for(std::size_t i=0;i<x.size();++i) {
        numerator=std::max(numerator,std::abs(x[i]));
        denominator=std::max(denominator,std::abs(reference[i]));
    }
    return numerator/denominator;
}

class DomMomentSchurSystemV44 {
public:
    using CellField=cfdx::core::Field<double,cfdx::core::Location::CELL>;
    using FaceField=cfdx::core::Field<double,cfdx::core::Location::FACE>;

    struct WallFace {
        std::size_t face=0,owner=0;
        std::size_t reflect_index=std::numeric_limits<std::size_t>::max();
        double area=0.,hplus=0.,hminus=0.,emissivity=1.,ib=0.;
    };

    DomMomentSchurSystemV44(const cfdx::core::Mesh& mesh,
                           const FvGeometry& geometry,
                           const CellField& temperature,
                           const std::vector<DiscreteDirection>& directions,
                           const DomWallBoundaryConditions& walls,
                           const RadiationTransportControls& controls,
                           std::vector<double> cell_absorption={},
                           std::vector<double> cell_scattering={})
        : mesh_(mesh),geo_(geometry),temperature_(temperature),
          directions_(directions),wall_bcs_(walls),controls_(controls),
          absorption_cell_(std::move(cell_absorption)),
          scattering_cell_(std::move(cell_scattering)),
          nc_(mesh.n_cells())
    {
        validate_radiation_transport_controls(controls_);
        validate_discrete_directions(directions_);
        validate_dom_wall_conditions(mesh_,wall_bcs_);
        if(nc_==0 || temperature_.size()!=nc_ ||
           !(controls_.absorption>0.) || !std::isfinite(controls_.absorption) ||
           !std::isfinite(controls_.scattering) ||
           geo_.face_area_vectors.size()!=mesh_.n_faces() ||
           geo_.face_centres.size()!=mesh_.n_faces() ||
           geo_.cell_centres.size()!=nc_ || geo_.cell_volumes.size()!=nc_)
            throw std::invalid_argument("DOM Schur V44: invalid geometry or transport data");
        // V50: optional cellwise coefficients use the exact same FV angular
        // matrices, not a postprocessed heterogeneous low-order approximation.
        // Both vectors are mandatory together. The historical uniform path
        // initializes the same numbers and preserves its original behaviour.
        const bool has_absorption=!absorption_cell_.empty();
        const bool has_scattering=!scattering_cell_.empty();
        if(has_absorption!=has_scattering ||
           (has_absorption && (absorption_cell_.size()!=nc_ ||
                               scattering_cell_.size()!=nc_)))
            throw std::invalid_argument("DOM Schur V50: mismatched cellwise optical fields");
        variable_optics_=has_absorption;
        if(!variable_optics_) {
            absorption_cell_.assign(nc_,controls_.absorption);
            scattering_cell_.assign(nc_,controls_.scattering);
        }
        for(std::size_t c=0;c<nc_;++c)
            if(!(absorption_cell_[c]>0.) || !std::isfinite(absorption_cell_[c]) ||
               !(scattering_cell_[c]>=0.) || !std::isfinite(scattering_cell_[c]) ||
               !std::isfinite(absorption_cell_[c]+scattering_cell_[c]))
                throw std::invalid_argument("DOM Schur V50: invalid cellwise optical coefficients");
        ib_.resize(nc_);
        for(std::size_t c=0;c<nc_;++c) {
            if (!(temperature_(c)>=0.) || !std::isfinite(temperature_(c)) ||
                !(geo_.cell_volumes[c]>0.) || !std::isfinite(geo_.cell_volumes[c]))
                throw std::invalid_argument("DOM Schur V44: invalid gas temperature or volume");
            ib_[c]=blackbody_intensity(temperature_(c));
        }
        std::vector<char> seen(mesh_.n_faces(),0);
        face_to_wall_.assign(mesh_.n_faces(),std::numeric_limits<std::size_t>::max());
        for(std::size_t p=0;p<mesh_.boundary().n_patches();++p) {
            const auto& patch=mesh_.boundary().patch(p);
            if(patch.type==cfdx::core::PatchType::EMPTY)
                throw std::invalid_argument("DOM Schur V44: EMPTY patches not qualified");
            const auto& bc=wall_bcs_.at(patch.name);
            for(auto f:patch.face_ids) {
                if(f>=mesh_.n_faces() || seen[f]++ ||
                   mesh_.ownership().neighbour(f)>=0)
                    throw std::invalid_argument("DOM Schur V44: invalid boundary partition");
                const auto sf=geo_.face_area_vectors[f];
                const double a=sf.mag();
                const auto owner=mesh_.ownership().owner(f);
                if(!(a>0.) || !std::isfinite(a) || owner>=nc_)
                    throw std::invalid_argument("DOM Schur V44: invalid boundary face");
                const auto hemi=dom_hemisphere_moments(directions_,{sf.x,sf.y,sf.z});
                WallFace wall;
                wall.face=f;wall.owner=owner;wall.area=a;
                wall.hplus=hemi.outgoing;wall.hminus=hemi.incoming;
                wall.emissivity=bc.emissivity;
                wall.ib=blackbody_intensity(bc.temperature);
                if(bc.emissivity<1.) wall.reflect_index=nreflect_++;
                face_to_wall_[f]=faces_.size();
                faces_.push_back(wall);
            }
        }
        if(faces_.empty()) throw std::invalid_argument("DOM Schur V44: no wall faces");
        // Assemble exactly the production bounded UPWIND angular operator.
        // Fixed incoming boundary values are set to zero for matrix assembly,
        // then emitted/reflected incident RHS terms are restored in project().
        CellField zero(nc_,"zero","W/m3/sr",1);
        CellField sp(nc_,"extinction","1/m",1);
        CellField old(nc_,"I","W/m2/sr",1);
        zero.fill(0.); old.fill(0.);
        for(std::size_t c=0;c<nc_;++c)
            sp(c)=-(absorption_cell_[c]+scattering_cell_[c]);
        ScalarBoundaryConditions bcs;
        ScalarBoundaryFaceValues zero_faces;
        for(std::size_t p=0;p<mesh_.boundary().n_patches();++p) {
            const auto& name=mesh_.boundary().patch(p).name;
            bcs[name]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
            auto& values=zero_faces.values[name];
            values.assign(mesh_.n_faces(),std::numeric_limits<double>::quiet_NaN());
            for(const auto f:mesh_.boundary().patch(p).face_ids) values[f]=0.;
        }
        angular_.reserve(directions_.size());
        for(const auto& d:directions_) {
            FaceField flux(mesh_.n_faces(),"s dot Sf","m2",1);
            for(std::size_t f=0;f<mesh_.n_faces();++f) {
                const auto sf=geo_.face_area_vectors[f];
                flux(f)=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
            }
            angular_.push_back(assemble_scalar_equation(
                mesh_,geo_,flux,0.,zero,sp,bcs,true,&zero_faces,
                nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&old));
        }
    }

    std::size_t n_cells() const { return nc_; }
    std::size_t n_reflecting_faces() const { return nreflect_; }
    std::size_t dimension() const { return nc_+nreflect_; }
    std::size_t transport_maps() const { return transport_maps_; }
    // Research-only access to the real UPWIND angular diagonal. It permits
    // a local moment-Jacobi approximation without probing full A_m^{-1}.
    double angular_cell_diagonal(std::size_t ordinate,std::size_t cell) const {
        if(ordinate>=angular_.size() || cell>=nc_)
            throw std::out_of_range("DOM Schur V44: angular diagonal index");
        return angular_[ordinate].matrix(cell,cell);
    }
    // V47: fail-closed identity checks for the HOLO geometry supplied by caller.
    const cfdx::core::Mesh& source_mesh() const { return mesh_; }
    const FvGeometry& source_geometry() const { return geo_; }
    const CellField& source_temperature() const { return temperature_; }
    const DomWallBoundaryConditions& source_walls() const { return wall_bcs_; }
    const std::vector<WallFace>& wall_faces() const { return faces_; }
    const std::vector<DiscreteDirection>& directions() const { return directions_; }
    const RadiationTransportControls& transport_controls() const { return controls_; }
    const std::vector<double>& cell_absorption() const {return absorption_cell_;}
    const std::vector<double>& cell_scattering() const {return scattering_cell_;}
    bool has_variable_optics() const {return variable_optics_;}

    // V57 research-only: read the *assembled* UPWIND angular CSR, including
    // its boundary diagonal. No mutation of production/legacy angular solves.
    const cfdx::core::SparseMatrix& angular_matrix_research_v57(
        std::size_t ordinate) const {
        if(ordinate>=angular_.size())
            throw std::out_of_range("V57: angular matrix ordinate");
        return angular_[ordinate].matrix;
    }

    // Computes M*x, or b+M*x, where M is the EXACT moment/wall projection
    // of the inverse angular FV matrices and b is the fixed emission term.
    // Signed intermediate intensities are allowed for linear Krylov actions;
    // the final physical state is screened by the independent V42 FV audit.
    std::vector<double> project(const std::vector<double>& x,bool emission,
                                std::vector<CellField>* captured=nullptr) {
        if(x.size()!=dimension())
            throw std::invalid_argument("DOM Schur V44: moment vector dimension");
        for(double v:x) if(!std::isfinite(v))
            throw std::invalid_argument("DOM Schur V44: nonfinite moment vector");
        std::vector<double> out(dimension(),0.);
        if(captured) {captured->clear();captured->reserve(directions_.size());}
        for(std::size_t m=0;m<directions_.size();++m) {
            auto& eq=angular_[m];
            for(std::size_t c=0;c<nc_;++c)
                eq.rhs(c)=geo_.cell_volumes[c]*(
                    scattering_cell_[c]*x[c]+(emission?absorption_cell_[c]*ib_[c]:0.));
            const auto& d=directions_[m];
            for(const auto& face:faces_) {
                const auto sf=geo_.face_area_vectors[face.face];
                const double normal_flux=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                if(normal_flux>=0.) continue;
                const double emission_intensity=emission?face.emissivity*face.ib:0.;
                const double reflection=face.reflect_index<nreflect_
                    ? (1.-face.emissivity)*(face.hplus/face.hminus)*
                      x[nc_+face.reflect_index] : 0.;
                eq.rhs(face.owner)+=(-normal_flux)*(emission_intensity+reflection);
            }
            cfdx::core::Vector intensity(nc_,0.);
            ScalarSolveControls lin;
            lin.tolerance=controls_.linear_tolerance;
            lin.max_iterations=controls_.linear_max_iterations;
            lin.relaxation=1.;
            const auto lr=solve_scalar_equation(eq,intensity,lin);
            if(lr.status!=cfdx::core::SolverStatus::CONVERGED)
                throw std::runtime_error("DOM Schur V44: angular linear solve failed");
            if(captured) captured->emplace_back(nc_,"I","W/m2/sr",1);
            for(std::size_t c=0;c<nc_;++c) {
                const double I=intensity(c);
                if(!std::isfinite(I))
                    throw std::runtime_error("DOM Schur V44: nonfinite angular solution");
                out[c]+=directions_[m].weight*I/(4.*M_PI);
                if(captured) captured->back()(c)=I;
            }
            for(const auto& face:faces_) if(face.reflect_index<nreflect_) {
                const auto sf=geo_.face_area_vectors[face.face];
                const double normal_flux=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                if(normal_flux>0.)
                    out[nc_+face.reflect_index]+=
                        d.weight*(normal_flux/face.area)/face.hplus*
                        intensity(face.owner);
            }
        }
        ++transport_maps_;
        for(double v:out) if(!std::isfinite(v))
            throw std::runtime_error("DOM Schur V44: nonfinite projected moment");
        return out;
    }

    std::vector<double> constant_rhs() {
        return project(std::vector<double>(dimension(),0.),true);
    }

    std::vector<double> apply_exact_schur(const std::vector<double>& x) {
        const auto response=project(x,false);
        std::vector<double> out=x;
        for(std::size_t i=0;i<x.size();++i) out[i]-=response[i];
        return out;
    }

    // Approximate diffusion preconditioner only, NOT the exact discrete Schur.
    // The wall extrapolation is intentionally recorded as a hypothesis rather
    // than presented as a proven high-order/low-order boundary equivalence.
    cfdx::core::SparseMatrix diffusion_matrix() const {
        const double ka=controls_.absorption,kt=ka+controls_.scattering;
        const double D=1./(3.*kt);
        cfdx::core::SparseMatrix B(nc_,nc_);
        std::vector<double> diag(nc_,0.);
        for(std::size_t c=0;c<nc_;++c)
            diag[c]+=ka*geo_.cell_volumes[c];
        const auto& own=mesh_.ownership();
        for(std::size_t f=0;f<mesh_.n_faces();++f) {
            const auto sf=geo_.face_area_vectors[f];
            const double area=sf.mag();
            if(!(area>0.) || !std::isfinite(area))
                throw std::invalid_argument("DOM Schur V44: invalid diffusion face area");
            const auto o=own.owner(f);
            const auto ni=own.neighbour(f);
            if(o>=nc_) throw std::invalid_argument("DOM Schur V44: invalid diffusion owner");
            if(ni>=0) {
                const auto n=static_cast<std::size_t>(ni);
                if(n>=nc_ || n==o)
                    throw std::invalid_argument("DOM Schur V44: invalid diffusion neighbour");
                const auto d=geo_.cell_centres[n]-geo_.cell_centres[o];
                const double normal=d.dot(sf)/area;
                if(!(normal>0.) || !std::isfinite(normal))
                    throw std::invalid_argument("DOM Schur V44: inverted interior face");
                const double k=D*area/normal;
                diag[o]+=k;diag[n]+=k;
                B.push_back(o,n,-k);B.push_back(n,o,-k);
            } else {
                const auto wall_id=face_to_wall_[f];
                if(wall_id>=faces_.size())
                    throw std::invalid_argument("DOM Schur V44: unassigned wall face");
                const auto d=geo_.face_centres[f]-geo_.cell_centres[o];
                const double normal=d.dot(sf)/area;
                if(!(normal>0.) || !std::isfinite(normal))
                    throw std::invalid_argument("DOM Schur V44: inverted boundary face");
                const double eps=faces_[wall_id].emissivity;
                if(eps>0.) {
                    const double extrap=2.*(2.-eps)/eps;
                    diag[o]+=area/(normal/D+extrap);
                }
            }
        }
        for(std::size_t c=0;c<nc_;++c) B.push_back(c,c,diag[c]);
        B.finalize();
        return B;
    }

    std::vector<double> precondition(const std::vector<double>& v,
                                     ScalarEquation& diffusion,double blend) {
        if(v.size()!=dimension())
            throw std::invalid_argument("DOM Schur V44: precondition dimension");
        const double scale=std::max(controls_.absorption,controls_.scattering);
        for(std::size_t c=0;c<nc_;++c)
            diffusion.rhs(c)=scale*geo_.cell_volumes[c]*v[c];
        cfdx::core::Vector x(nc_,0.);
        ScalarSolveControls lin;
        lin.tolerance=controls_.linear_tolerance;
        lin.max_iterations=controls_.linear_max_iterations;
        lin.relaxation=1.;
        const auto lr=solve_scalar_equation(diffusion,x,lin);
        if(lr.status!=cfdx::core::SolverStatus::CONVERGED)
            throw std::runtime_error("DOM Schur V44: diffusion preconditioner failed");
        std::vector<double> z=v;
        for(std::size_t c=0;c<nc_;++c)
            z[c]=(1.-blend)*v[c]+blend*x(c);
        return z;
    }

    // No mutation of external gas/radiation fields occurs in this class.
    DomMomentSchurV44Result solve(bool use_diffusion_preconditioner,
                                 std::size_t restart=12,
                                 std::size_t max_krylov_iterations=100,
                                 double diffusion_mix=1.,
        const std::function<std::vector<double>(const std::vector<double>&)>&
            custom_preconditioner={},
        // Experimental V46 callback checks the candidate *reconstructed* FV
        // angular state (including the exact face-current energy balance).
        // The default is empty, which preserves V44/V45 stopping behaviour.
        const std::function<bool(const std::vector<double>&,
            const std::vector<CellField>&)>& physical_gate_v46={}) {
        if(variable_optics_)
            throw std::invalid_argument("DOM Schur V50: legacy physical solve assumes uniform optics; use V50 variable-audit driver");
        if(restart<2 || restart>40 || max_krylov_iterations==0 ||
           !std::isfinite(diffusion_mix) || diffusion_mix<0. || diffusion_mix>1. ||
           (use_diffusion_preconditioner && custom_preconditioner))
            throw std::invalid_argument("DOM Schur V44: invalid GMRES controls");
        transport_maps_=0;
        const auto t0=std::chrono::steady_clock::now();
        const auto b=constant_rhs();
        const auto n=dimension();
        const auto rhs_l2=std::max(1.,dom_schur_norm_v44(b));
        double gate=controls_.tolerance*0.05;
        std::vector<double> x(n,0.);
        DomMomentSchurV44Result result;
        result.reflecting_faces=nreflect_;
        ScalarEquation diffusion;
        if(use_diffusion_preconditioner) {
            diffusion.matrix=diffusion_matrix();
            diffusion.rhs=cfdx::core::Vector(nc_,0.);
        }
        std::size_t steps=0;
        // Tighten only the algebraic target after the *true* FV physics
        // rejects an otherwise converged candidate. A failed state is never
        // published. No production or historical threshold is relaxed.
        for (;;) {
        while(steps<max_krylov_iterations) {
            auto Ax=apply_exact_schur(x);
            std::vector<double> r(n,0.);
            for(std::size_t i=0;i<n;++i) r[i]=b[i]-Ax[i];
            const double beta=dom_schur_norm_v44(r);
            result.true_schur_relative=dom_schur_max_relative_v44(r,b);
            if(result.true_schur_relative<=gate) break;
            if(!(beta>0.) || !std::isfinite(beta))
                throw std::runtime_error("DOM Schur V44: invalid true residual");
            const std::size_t m=std::min(restart,max_krylov_iterations-steps);
            std::vector<std::vector<double>> V(m+1,std::vector<double>(n,0.));
            std::vector<std::vector<double>> Z(m,std::vector<double>(n,0.));
            std::vector<std::vector<double>> H(m+1,std::vector<double>(m,0.));
            std::vector<double> cs(m,0.),sn(m,0.),g(m+1,0.);
            for(std::size_t i=0;i<n;++i) V[0][i]=r[i]/beta;
            g[0]=beta;
            std::size_t used=0;
            for(std::size_t j=0;j<m;++j) {
                Z[j]=custom_preconditioner ? custom_preconditioner(V[j]) :
                    (use_diffusion_preconditioner ?
                        precondition(V[j],diffusion,diffusion_mix) : V[j]);
                if(Z[j].size()!=n) throw std::invalid_argument("DOM Schur V44: custom precondition dimension");
                for(double value:Z[j]) if(!std::isfinite(value))
                    throw std::invalid_argument("DOM Schur V44: custom precondition NaN");
                if(use_diffusion_preconditioner) ++result.diffusion_solves;
                auto w=apply_exact_schur(Z[j]);
                // Reorthogonalized modified Gram-Schmidt for ill-conditioned
                // high-albedo FV Schur operators.
                for(int pass=0;pass<2;++pass)
                    for(std::size_t i=0;i<=j;++i) {
                        const double h=dom_schur_dot_v44(w,V[i]);
                        H[i][j]+=h;
                        for(std::size_t k=0;k<n;++k) w[k]-=h*V[i][k];
                    }
                H[j+1][j]=dom_schur_norm_v44(w);
                const bool breakdown=H[j+1][j]<=1e-14*rhs_l2;
                if(!breakdown)
                    for(std::size_t k=0;k<n;++k) V[j+1][k]=w[k]/H[j+1][j];
                for(std::size_t i=0;i<j;++i) {
                    const double h=cs[i]*H[i][j]+sn[i]*H[i+1][j];
                    H[i+1][j]=-sn[i]*H[i][j]+cs[i]*H[i+1][j];
                    H[i][j]=h;
                }
                const double rho=std::hypot(H[j][j],H[j+1][j]);
                if(!(rho>0.) || !std::isfinite(rho))
                    throw std::runtime_error("DOM Schur V44: singular Arnoldi basis");
                cs[j]=H[j][j]/rho;sn[j]=H[j+1][j]/rho;
                H[j][j]=rho;H[j+1][j]=0.;
                g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
                used=j+1;++steps;
                if(std::abs(g[j+1])<=gate*rhs_l2 || breakdown) break;
            }
            std::vector<double> y(used,0.);
            for(std::size_t rev=used;rev-- >0;) {
                double s=g[rev];
                for(std::size_t k=rev+1;k<used;++k) s-=H[rev][k]*y[k];
                if(!(std::abs(H[rev][rev])>1e-25))
                    throw std::runtime_error("DOM Schur V44: Arnoldi breakdown");
                y[rev]=s/H[rev][rev];
            }
            for(std::size_t i=0;i<used;++i)
                for(std::size_t k=0;k<n;++k) x[k]+=y[i]*Z[i][k];
        }
        result.krylov_iterations=steps;
        const auto Ax=apply_exact_schur(x);
        std::vector<double> residual(n,0.);
        for(std::size_t i=0;i<n;++i) residual[i]=b[i]-Ax[i];
        result.true_schur_relative=dom_schur_max_relative_v44(residual,b);
        const bool algebraic_gate_met=result.true_schur_relative<=gate;
        if(algebraic_gate_met) {
            std::vector<CellField> angular;
            const auto mapped=project(x,true,&angular);
            std::vector<double> fixed_defect(n,0.);
            for(std::size_t i=0;i<n;++i) fixed_defect[i]=mapped[i]-x[i];
            result.true_schur_relative=std::max(result.true_schur_relative,
                dom_schur_max_relative_v44(fixed_defect,b));
            if(result.true_schur_relative<=controls_.tolerance) {
                // The V42 audit uses the actual diffuse-gray wall fields.
                const auto audit=audit_dom_accepted_state_v42(
                    mesh_,geo_,temperature_,angular,directions_,controls_,wall_bcs_);
                result.local_transport_relative=audit.max_cell_relative;
                result.global_energy_relative=audit.relative_global_energy;
                if(audit.max_cell_relative<=controls_.tolerance &&
                   audit.relative_global_energy<=controls_.tolerance) {
                    bool positive=true;
                    for(const auto& a:angular)
                        for(std::size_t c=0;c<nc_;++c)
                            if(a(c)<-1e-10 || !std::isfinite(a(c))) positive=false;
                    for(double v:x) if(v<0. || !std::isfinite(v)) positive=false;
                    if(positive) {
                        const bool external_ok=!physical_gate_v46 ||
                            physical_gate_v46(x,angular);
                        if(external_ok) {
                            result.converged=true;
                            result.scalar_mean_intensity.reserve(nc_);
                            result.reflected_face_mean_intensity.reserve(nreflect_);
                            for(std::size_t c=0;c<nc_;++c)
                                result.scalar_mean_intensity.push_back(x[c]);
                            for(std::size_t f=0;f<nreflect_;++f)
                                result.reflected_face_mean_intensity.push_back(x[nc_+f]);
                        }
                    }
                }
            }
        }
        if(result.converged || !physical_gate_v46 ||
           !algebraic_gate_met || steps>=max_krylov_iterations ||
           result.physical_refinements_v46>=7 ||
           gate<=std::max(1.e-15,controls_.tolerance*1.e-5)) break;
        ++result.physical_refinements_v46;
        ++result.physical_rejections_v46;
        gate=std::max(1.e-15,gate*0.1);
        }
        result.transport_maps=transport_maps_;
        result.seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-t0).count();
        return result;
    }

private:
    const cfdx::core::Mesh& mesh_;
    const FvGeometry& geo_;
    const CellField& temperature_;
    const std::vector<DiscreteDirection>& directions_;
    const DomWallBoundaryConditions& wall_bcs_;
    RadiationTransportControls controls_;
    std::vector<double> absorption_cell_,scattering_cell_;
    bool variable_optics_=false;
    std::size_t nc_=0,nreflect_=0,transport_maps_=0;
    std::vector<double> ib_;
    std::vector<WallFace> faces_;
    std::vector<std::size_t> face_to_wall_;
    std::vector<ScalarEquation> angular_;
};

} // namespace cfdx::physics
