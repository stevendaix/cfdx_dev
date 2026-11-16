#pragma once
// V45 opt-in diagnostic. All currents use the actual oriented faces and
// angular UPWIND state of V44; the low-order drift is an exact identity at
// a supplied high-order state, NOT a proof of diffusion synthetic acceleration.
#include "cfdx/physics/dom_moment_schur_v44.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct DomFaceMomentAuditV45 {
    std::vector<double> high_order_face_current;
    std::vector<double> diffusive_face_current;
    std::vector<double> closure_face_current;
    std::vector<double> cell_mean_intensity;
    double max_transport_balance_relative=0.;
    double max_physical_balance_relative=0.;
    double max_low_order_identity_relative=0.;
    double global_energy_relative=0.;
    double max_moment_mismatch_relative=0.;
};

// For any candidate (J,R), project() gives the true angular FV solution.
// The current identity including the scattering lag is:
//   div(q_HO) + 4*pi*V*(kt*J_HO-ks*J_candidate-ka*Ib) = 0.
// If J_candidate == J_HO, the physical zero-th moment equation follows.
// The *same* q_HO is shared by owner and neighbour with opposite signs.
inline DomFaceMomentAuditV45 audit_dom_face_moments_v45(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geo,
    const DomMomentSchurSystemV44::CellField& temperature,
    const std::vector<DiscreteDirection>& directions,
    const DomWallBoundaryConditions& walls,
    const RadiationTransportControls& controls,
    const DomMomentSchurSystemV44& system,
    const std::vector<double>& candidate,
    const std::vector<DomMomentSchurSystemV44::CellField>& intensity)
{
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(candidate.size()!=system.dimension() || intensity.size()!=directions.size() ||
       temperature.size()!=nc || geo.cell_volumes.size()!=nc ||
       geo.cell_centres.size()!=nc || geo.face_centres.size()!=nf ||
       geo.face_area_vectors.size()!=nf)
        throw std::invalid_argument("DOM V45: inconsistent face-current dimensions");
    validate_discrete_directions(directions);
    validate_dom_wall_conditions(mesh,walls);
    const auto& reference=system.transport_controls();
    const auto& reference_angles=system.directions();
    if(reference.absorption!=controls.absorption ||
       reference.scattering!=controls.scattering ||
       reference_angles.size()!=directions.size())
        throw std::invalid_argument("DOM V45: transport configuration differs from Schur");
    for(std::size_t m=0;m<directions.size();++m) {
        const auto& a=reference_angles[m],&b=directions[m];
        if(a.dx!=b.dx || a.dy!=b.dy || a.dz!=b.dz || a.weight!=b.weight)
            throw std::invalid_argument("DOM V45: angular quadrature differs from Schur");
    }
    for(double v:candidate) if(!std::isfinite(v))
        throw std::invalid_argument("DOM V45: nonfinite moment input");
    const double pi=std::acos(-1.),fourpi=4.*pi;
    DomFaceMomentAuditV45 report;
    report.high_order_face_current.assign(nf,0.);
    report.diffusive_face_current.assign(nf,0.);
    report.closure_face_current.assign(nf,0.);
    report.cell_mean_intensity.assign(nc,0.);
    std::vector<const DomMomentSchurSystemV44::WallFace*> boundary(nf,nullptr);
    for(const auto& w:system.wall_faces()) {
        if(w.face>=nf || boundary[w.face])
            throw std::invalid_argument("DOM V45: invalid wall partition");
        boundary[w.face]=&w;
    }
    for(std::size_t m=0;m<directions.size();++m) {
        const auto& d=directions[m];
        if(intensity[m].size()!=nc)
            throw std::invalid_argument("DOM V45: intensity cell count");
        for(std::size_t c=0;c<nc;++c) {
            if(!std::isfinite(intensity[m](c)))
                throw std::invalid_argument("DOM V45: nonfinite angular intensity");
            report.cell_mean_intensity[c]+=d.weight*intensity[m](c)/fourpi;
        }
        for(std::size_t f=0;f<nf;++f) {
            const auto& sf=geo.face_area_vectors[f];
            const double normal_flux=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
            const auto owner=mesh.ownership().owner(f);
            const auto neighbour=mesh.ownership().neighbour(f);
            if(owner>=nc || !std::isfinite(normal_flux))
                throw std::invalid_argument("DOM V45: invalid oriented face");
            double upwind=0.;
            if(normal_flux>=0.) upwind=intensity[m](owner);
            else if(neighbour>=0) {
                if(static_cast<std::size_t>(neighbour)>=nc)
                    throw std::invalid_argument("DOM V45: invalid face neighbour");
                upwind=intensity[m](static_cast<std::size_t>(neighbour));
            } else {
                const auto* w=boundary[f];
                if(!w || !(w->hminus>0.) || !(w->hplus>0.))
                    throw std::invalid_argument("DOM V45: unassigned boundary face");
                const double reflect=w->reflect_index<system.n_reflecting_faces()
                    ? (1.-w->emissivity)*w->hplus/w->hminus*
                      candidate[nc+w->reflect_index] : 0.;
                upwind=w->emissivity*w->ib+reflect;
            }
            report.high_order_face_current[f]+=d.weight*normal_flux*upwind;
        }
    }
    std::vector<double> div(nc,0.),div_diff(nc,0.),div_close(nc,0.);
    const double D=1./(3.*(controls.absorption+controls.scattering));
    for(std::size_t f=0;f<nf;++f) {
        const auto& sf=geo.face_area_vectors[f];
        const double area=sf.mag();
        const auto o=mesh.ownership().owner(f);
        const auto nei=mesh.ownership().neighbour(f);
        if(!(area>0.) || !std::isfinite(area))
            throw std::invalid_argument("DOM V45: invalid face area");
        const auto nvec=sf*(1./area);
        double current_diff=0.;
        if(nei>=0) {
            const auto n=static_cast<std::size_t>(nei);
            const double dist=(geo.cell_centres[n]-geo.cell_centres[o]).dot(nvec);
            if(!(dist>0.)) throw std::invalid_argument("DOM V45: inverted internal face");
            current_diff=4.*pi*D*area*(candidate[o]-candidate[n])/dist;
        } else {
            const auto* w=boundary[f];
            if(!w) throw std::invalid_argument("DOM V45: missing wall face");
            const double dist=(geo.face_centres[f]-geo.cell_centres[o]).dot(nvec);
            if(!(dist>0.)) throw std::invalid_argument("DOM V45: inverted wall face");
            if(w->emissivity>0.) {
                const double extrap=2.*(2.-w->emissivity)/w->emissivity;
                current_diff=4.*pi*area*(candidate[o]-w->ib)/(dist/D+extrap);
            }
        }
        const double current=report.high_order_face_current[f];
        report.diffusive_face_current[f]=current_diff;
        report.closure_face_current[f]=current-current_diff;
        div[o]+=current;div_diff[o]+=current_diff;
        div_close[o]+=current-current_diff;
        if(nei>=0) {
            const auto n=static_cast<std::size_t>(nei);
            div[n]-=current;div_diff[n]-=current_diff;
            div_close[n]-=current-current_diff;
        }
    }
    long double total_physical=0.,total_scale=0.;
    double exterior=0.;
    for(std::size_t f=0;f<nf;++f)
        if(mesh.ownership().neighbour(f)<0) exterior+=report.high_order_face_current[f];
    for(std::size_t c=0;c<nc;++c) {
        const double V=geo.cell_volumes[c];
        if(!(V>0.) || !std::isfinite(V) || !std::isfinite(temperature(c)) ||
           !std::isfinite(candidate[c]))
            throw std::invalid_argument("DOM V45: invalid volume, state, temperature");
        const double ib=blackbody_intensity(temperature(c));
        const double j=report.cell_mean_intensity[c];
        const double base=fourpi*V*controls.absorption*(j-ib);
        const double lag=fourpi*V*controls.scattering*(j-candidate[c]);
        const double transport=div[c]+base+lag;
        const double physical=div[c]+base;
        const double low_order=div_diff[c]+div_close[c]+base;
        const double scale=std::max({1.,std::abs(div[c]),std::abs(base),
            fourpi*V*(controls.absorption+controls.scattering)*std::abs(j)});
        report.max_transport_balance_relative=std::max(
            report.max_transport_balance_relative,std::abs(transport)/scale);
        report.max_physical_balance_relative=std::max(
            report.max_physical_balance_relative,std::abs(physical)/scale);
        report.max_low_order_identity_relative=std::max(
            report.max_low_order_identity_relative,std::abs(low_order-physical)/scale);
        report.max_moment_mismatch_relative=std::max(
            report.max_moment_mismatch_relative,std::abs(j-candidate[c])/
            std::max(1.,std::abs(j)));
        total_physical+=base;
        total_scale+=std::abs(base);
    }
    total_physical+=exterior;
    total_scale+=std::abs(exterior)+1.;
    report.global_energy_relative=std::abs(static_cast<double>(total_physical)) /
        static_cast<double>(total_scale);
    return report;
}

// A separate fail-closed physical gate for this face-current oracle. It is
// deliberately not loosened to accommodate a high-scattering fixture, whose
// independent global normalization can be more demanding than the V42 one.
inline bool dom_face_moment_qualified_v45(const DomFaceMomentAuditV45& r,
                                         double tolerance) {
    return std::isfinite(tolerance) && tolerance>0. &&
        std::isfinite(r.max_transport_balance_relative) &&
        std::isfinite(r.max_physical_balance_relative) &&
        std::isfinite(r.max_low_order_identity_relative) &&
        std::isfinite(r.global_energy_relative) &&
        std::isfinite(r.max_moment_mismatch_relative) &&
        r.max_transport_balance_relative<=tolerance &&
        r.max_physical_balance_relative<=tolerance &&
        r.max_low_order_identity_relative<=tolerance &&
        r.global_energy_relative<=tolerance &&
        r.max_moment_mismatch_relative<=tolerance;
}

// Small-problem diagnostic: probe the EXACT V44 Schur columns and invert
// independent (cell + its reflecting faces) diagonal blocks. This samples
// actual transport/gray reflections rather than imposing continuous P1 wall
// coefficients. O(n) angular maps and O(n*n) storage are intentionally capped;
// it is NOT a scalable production DSA preconditioner.
class DomLocalSchurBlocksV45 {
public:
    struct Block {
        std::vector<std::size_t> ids;
        std::vector<double> lu;
        std::vector<std::size_t> swaps;
    };
    explicit DomLocalSchurBlocksV45(DomMomentSchurSystemV44& system,
                                    std::size_t max_dimension=128)
        : n_(system.dimension()) {
        if(n_==0 || n_>max_dimension || max_dimension>256)
            throw std::invalid_argument("DOM V45: diagnostic dense Schur cap exceeded");
        const auto start=std::chrono::steady_clock::now();
        blocks_.resize(system.n_cells());
        for(std::size_t c=0;c<system.n_cells();++c)
            blocks_[c].ids.push_back(c);
        for(const auto& f:system.wall_faces())
            if(f.reflect_index<system.n_reflecting_faces())
                blocks_.at(f.owner).ids.push_back(system.n_cells()+f.reflect_index);
        dense_.assign(n_*n_,0.);
        std::vector<double> e(n_,0.);
        for(std::size_t col=0;col<n_;++col) {
            e[col]=1.;
            const auto output=system.apply_exact_schur(e);
            e[col]=0.;
            for(std::size_t row=0;row<n_;++row)
                dense_[row*n_+col]=output[row];
        }
        probe_maps_=n_;
        for(auto& block:blocks_) {
            const auto p=block.ids.size();
            block.lu.resize(p*p);
            block.swaps.resize(p);
            for(std::size_t i=0;i<p;++i)
                for(std::size_t j=0;j<p;++j)
                    block.lu[i*p+j]=dense_[block.ids[i]*n_+block.ids[j]];
            for(std::size_t k=0;k<p;++k) {
                std::size_t pivot=k;
                for(std::size_t r=k+1;r<p;++r)
                    if(std::abs(block.lu[r*p+k])>std::abs(block.lu[pivot*p+k])) pivot=r;
                const double mag=std::abs(block.lu[pivot*p+k]);
                if(!(mag>1.e-12) || !std::isfinite(mag))
                    throw std::runtime_error("DOM V45: singular local Schur block");
                block.swaps[k]=pivot;
                if(pivot!=k)
                    for(std::size_t j=0;j<p;++j)
                        std::swap(block.lu[k*p+j],block.lu[pivot*p+j]);
                for(std::size_t i=k+1;i<p;++i) {
                    block.lu[i*p+k]/=block.lu[k*p+k];
                    for(std::size_t j=k+1;j<p;++j)
                        block.lu[i*p+j]-=block.lu[i*p+k]*block.lu[k*p+j];
                }
            }
        }
        setup_seconds_=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
    }
    std::size_t dimension() const {return n_;}
    std::size_t probe_maps() const {return probe_maps_;}
    std::size_t apply_calls() const {return apply_calls_;}
    double setup_seconds() const {return setup_seconds_;}
    const std::vector<double>& dense_schur() const {return dense_;}
    const std::vector<Block>& blocks() const {return blocks_;}
    std::vector<double> apply(const std::vector<double>& r) const {
        if(r.size()!=n_) throw std::invalid_argument("DOM V45: precondition input size");
        std::vector<double> z(n_,0.);
        for(const auto& block:blocks_) {
            const auto p=block.ids.size();
            std::vector<double> b(p);
            for(std::size_t i=0;i<p;++i) {
                b[i]=r[block.ids[i]];
                if(!std::isfinite(b[i]))
                    throw std::invalid_argument("DOM V45: precondition NaN");
            }
            for(std::size_t k=0;k<p;++k)
                if(block.swaps[k]!=k) std::swap(b[k],b[block.swaps[k]]);
            for(std::size_t i=0;i<p;++i)
                for(std::size_t j=0;j<i;++j) b[i]-=block.lu[i*p+j]*b[j];
            for(std::size_t ii=p;ii-- >0;) {
                for(std::size_t j=ii+1;j<p;++j) b[ii]-=block.lu[ii*p+j]*b[j];
                b[ii]/=block.lu[ii*p+ii];
                if(!std::isfinite(b[ii]))
                    throw std::runtime_error("DOM V45: block solve nonfinite");
            }
            for(std::size_t i=0;i<p;++i) z[block.ids[i]]=b[i];
        }
        ++apply_calls_;
        return z;
    }
    // Diagnostics, not an eigenvalue claim for nonnormal operators.
    double max_row_sum_I_minus_preconditioned_H() const {
        std::vector<std::vector<double>> columns(n_);
        for(std::size_t col=0;col<n_;++col) {
            std::vector<double> h(n_);
            for(std::size_t row=0;row<n_;++row) h[row]=dense_[row*n_+col];
            columns[col]=apply(h);
        }
        double max_row=0.;
        for(std::size_t row=0;row<n_;++row) {
            double s=0.;
            for(std::size_t col=0;col<n_;++col)
                s+=std::abs((row==col?1.:0.)-columns[col][row]);
            max_row=std::max(max_row,s);
        }
        return max_row;
    }
private:
    std::size_t n_=0,probe_maps_=0;
    double setup_seconds_=0.;
    mutable std::size_t apply_calls_=0;
    std::vector<double> dense_;
    std::vector<Block> blocks_;
};


// Cheap local approximation using only the diagonal of each actual FV
// angular transport matrix. For each cell we eliminate the angular unknowns
// locally and retain its J and all gray-wall R variables on that cell.
// Unlike the sampled Schur blocks above, setup performs ZERO angular solves.
// Its validity/performance must be checked against the exact Schur residual.
class DomLocalAngularJacobiV45 {
public:
    using Block=DomLocalSchurBlocksV45::Block;
    DomLocalAngularJacobiV45(
        const DomMomentSchurSystemV44& system,
        const FvGeometry& geo,
        const std::vector<DiscreteDirection>& directions,
        const RadiationTransportControls& controls)
        : n_(system.dimension()) {
        const auto start=std::chrono::steady_clock::now();
        if(n_==0 || directions.empty() || geo.cell_volumes.size()!=system.n_cells())
            throw std::invalid_argument("DOM V45: invalid local angular setup");
        const auto& reference=system.transport_controls();
        const auto& reference_dirs=system.directions();
        if(reference.absorption!=controls.absorption ||
           reference.scattering!=controls.scattering ||
           reference_dirs.size()!=directions.size())
            throw std::invalid_argument("DOM V45: local block optical data mismatch");
        for(std::size_t m=0;m<directions.size();++m) {
            const auto& a=reference_dirs[m],&b=directions[m];
            if(a.dx!=b.dx || a.dy!=b.dy || a.dz!=b.dz || a.weight!=b.weight)
                throw std::invalid_argument("DOM V45: local block quadrature mismatch");
        }
        blocks_.resize(system.n_cells());
        std::vector<std::vector<const DomMomentSchurSystemV44::WallFace*>> walls(
            system.n_cells());
        for(const auto& w:system.wall_faces())
            if(w.reflect_index<system.n_reflecting_faces())
                walls.at(w.owner).push_back(&w);
        for(std::size_t c=0;c<system.n_cells();++c) {
            auto& block=blocks_[c];
            block.ids.push_back(c);
            for(const auto* f:walls[c])
                block.ids.push_back(system.n_cells()+f->reflect_index);
            const auto p=block.ids.size();
            block.lu.assign(p*p,0.);
            block.swaps.resize(p);
            for(std::size_t i=0;i<p;++i) block.lu[i*p+i]=1.;
            for(std::size_t m=0;m<directions.size();++m) {
                const auto& d=directions[m];
                const double a=system.angular_cell_diagonal(m,c);
                if(!(a>0.) || !std::isfinite(a))
                    throw std::invalid_argument("DOM V45: nonpositive angular diagonal");
                std::vector<double> response(p,0.);
                response[0]=controls.scattering*geo.cell_volumes[c]/a;
                for(std::size_t j=1;j<p;++j) {
                    const auto* face=walls[c][j-1];
                    const auto& sf=geo.face_area_vectors.at(face->face);
                    const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                    response[j]=mu<0. ? (-mu)*(1.-face->emissivity)*
                        face->hplus/face->hminus/a : 0.;
                }
                for(std::size_t j=0;j<p;++j)
                    block.lu[j]-=d.weight*response[j]/(4.*M_PI);
                for(std::size_t i=1;i<p;++i) {
                    const auto* face=walls[c][i-1];
                    const auto& sf=geo.face_area_vectors.at(face->face);
                    const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                    if(mu<=0.) continue;
                    const double weight=d.weight*(mu/face->area)/face->hplus;
                    for(std::size_t j=0;j<p;++j)
                        block.lu[i*p+j]-=weight*response[j];
                }
            }
            // Partial-pivot LU for the cell+wall moment block.
            for(std::size_t k=0;k<p;++k) {
                std::size_t pivot=k;
                for(std::size_t r=k+1;r<p;++r)
                    if(std::abs(block.lu[r*p+k])>std::abs(block.lu[pivot*p+k]))
                        pivot=r;
                const double mag=std::abs(block.lu[pivot*p+k]);
                if(!(mag>1e-12) || !std::isfinite(mag))
                    throw std::runtime_error("DOM V45: local angular Jacobi singular");
                block.swaps[k]=pivot;
                if(pivot!=k)
                    for(std::size_t j=0;j<p;++j)
                        std::swap(block.lu[k*p+j],block.lu[pivot*p+j]);
                for(std::size_t i=k+1;i<p;++i) {
                    block.lu[i*p+k]/=block.lu[k*p+k];
                    for(std::size_t j=k+1;j<p;++j)
                        block.lu[i*p+j]-=block.lu[i*p+k]*block.lu[k*p+j];
                }
            }
        }
        setup_seconds_=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
    }
    std::vector<double> apply(const std::vector<double>& r) const {
        if(r.size()!=n_) throw std::invalid_argument("DOM V45: cheap block dimension");
        std::vector<double> z(n_,0.);
        for(const auto& block:blocks_) {
            const auto p=block.ids.size();
            std::vector<double> b(p);
            for(std::size_t i=0;i<p;++i) {
                b[i]=r[block.ids[i]];
                if(!std::isfinite(b[i]))
                    throw std::invalid_argument("DOM V45: cheap block nonfinite RHS");
            }
            for(std::size_t k=0;k<p;++k)
                if(block.swaps[k]!=k) std::swap(b[k],b[block.swaps[k]]);
            for(std::size_t i=0;i<p;++i)
                for(std::size_t j=0;j<i;++j) b[i]-=block.lu[i*p+j]*b[j];
            for(std::size_t ii=p;ii-- >0;) {
                for(std::size_t j=ii+1;j<p;++j) b[ii]-=block.lu[ii*p+j]*b[j];
                b[ii]/=block.lu[ii*p+ii];
                if(!std::isfinite(b[ii]))
                    throw std::runtime_error("DOM V45: cheap block solve invalid");
            }
            for(std::size_t i=0;i<p;++i) z[block.ids[i]]=b[i];
        }
        ++apply_calls_;
        return z;
    }
    double setup_seconds() const {return setup_seconds_;}
    std::size_t apply_calls() const {return apply_calls_;}
private:
    std::size_t n_=0;
    double setup_seconds_=0.;
    mutable std::size_t apply_calls_=0;
    std::vector<Block> blocks_;
};

} // namespace cfdx::physics
