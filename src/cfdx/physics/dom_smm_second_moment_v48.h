#pragma once
// CFDX V48: research-only FV second-moment operator on the actual mesh graph.
// With P_ij = E_ij J, the low-order first-moment closure is
//     Q_f = -4*pi * Sf . div(P)_f / kappa_t,f .
// Here E is the *refreshed* HO Eddington tensor and the spatial derivative
// operates on P=E*J (NOT merely E*grad J). This is essential when E varies.
// Boundary currents must be supplied explicitly by a wall oracle (HO DOM);
// an independently discretized Marshak condition is NOT silently substituted.
// This matrix-free, gradient-corrected FV operator is not yet a proved
// coercive SMM solver or an industrial DSA. Production is unchanged.
#include "cfdx/physics/dom_holo_face_closure_v47.h"
#include "cfdx/physics/tensor_heat_flux_v37.h"
#include "cfdx/core/fvm/least_squares_gradient.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct DomSmmFluxV48 {
    // Each internal face has exactly one oriented current, positive outward
    // from its owner. Cell summation subtracts the same number at neighbour.
    std::vector<double> face_watts;
    std::vector<double> divergence_watts;
    double internal_cancellation_watts=0.;
};
struct DomSmmConsistencyV48 {
    DomSmmFluxV48 second_moment;
    std::vector<double> face_ho_correction;
    std::vector<double> reconstructed_ho_current;
    std::vector<double> reconstructed_divergence;
    double max_face_identity_relative=0.;
    double max_cell_identity_relative=0.;
};

class DomSmmSecondMomentV48 {
public:
    using Tensor=std::array<double,6>; // xx,yy,zz,xy,xz,yz
    DomSmmSecondMomentV48(const cfdx::core::Mesh& mesh,const FvGeometry& geo)
        :mesh_(mesh),geo_(geo) {
        const auto nc=mesh_.n_cells(),nf=mesh_.n_faces();
        if(nc<8 || !nf || geo_.cell_centres.size()!=nc ||
           geo_.cell_volumes.size()!=nc || geo_.face_centres.size()!=nf ||
           geo_.face_area_vectors.size()!=nf)
            throw std::invalid_argument("V48 SMM: invalid three-dimensional FV mesh");
        neighbour_cells_.resize(nc);
        for(std::size_t c=0;c<nc;++c) {
            if(!(geo_.cell_volumes[c]>0.)|| !finite(geo_.cell_centres[c]))
                throw std::invalid_argument("V48 SMM: invalid cell geometry");
        }
        for(std::size_t f=0;f<nf;++f) {
            const auto o=mesh_.ownership().owner(f);
            const auto n=mesh_.ownership().neighbour(f);
            if(o>=nc || !finite(geo_.face_centres[f]) ||
               !finite(geo_.face_area_vectors[f]) ||
               !(geo_.face_area_vectors[f].mag()>0.))
                throw std::invalid_argument("V48 SMM: invalid oriented face");
            if(n>=0) {
                const auto j=static_cast<std::size_t>(n);
                if(j>=nc || j==o || !((geo_.cell_centres[j]-geo_.cell_centres[o])
                                             .dot(geo_.face_area_vectors[f])>0.))
                    throw std::invalid_argument("V48 SMM: inverted or collapsed cell face");
                neighbour_cells_[o].push_back(j);
                neighbour_cells_[j].push_back(o);
            } else {
                const auto d=geo_.face_centres[f]-geo_.cell_centres[o];
                if(!(d.dot(geo_.face_area_vectors[f])>0.))
                    throw std::invalid_argument("V48 SMM: inverted boundary face");
            }
        }
        // Preserve the one-ring stencil wherever it is already full-rank.
        // Boundary tetrahedra may have only one or two cell neighbours;
        // use the topological second/third ring only for those cells. No
        // artificial boundary samples or rank-2 fallbacks are introduced.
        const auto direct_neighbours=neighbour_cells_;
        for(std::size_t c=0;c<nc;++c) {
            auto& nbr=neighbour_cells_[c];
            std::sort(nbr.begin(),nbr.end());
            nbr.erase(std::unique(nbr.begin(),nbr.end()),nbr.end());
            bool full_rank=false;
            for(unsigned ring=0;ring<3;++ring) {
                if(nbr.size()>=3) {
                    std::vector<cfdx::core::Vec3> xyz;
                    std::vector<double> constant;
                    xyz.reserve(nbr.size());constant.reserve(nbr.size());
                    for(const auto j:nbr) {
                        xyz.push_back(geo_.cell_centres[j]);constant.push_back(1.);
                    }
                    cfdx::core::StencilQuality quality;
                    (void)cfdx::core::weighted_least_squares_gradient_robust(
                        geo_.cell_centres[c],1.,xyz,constant,
                        cfdx::core::GradientWeighting::INVERSE_DISTANCE_SQUARED,
                        0,&quality);
                    if(quality.rank==3) {full_rank=true;break;}
                }
                if(ring==2)break;
                auto extended=nbr;
                for(const auto j:nbr)
                    for(const auto k:direct_neighbours[j])
                        if(k!=c)extended.push_back(k);
                std::sort(extended.begin(),extended.end());
                extended.erase(std::unique(extended.begin(),extended.end()),extended.end());
                if(extended==nbr)break;
                nbr=std::move(extended);
            }
            if(!full_rank)
                throw std::invalid_argument("V52 SMM: no full-rank 3D polyhedral WLS stencil");
        }

    }

    // Distinct from an HO heterogeneous-gas solve: variable kappa_t is
    // supported only in THIS LO operator. In an HOLO call the caller must
    // separately establish that the underlying angular transport matches it.
    DomSmmFluxV48 evaluate(const std::vector<double>& J,
                           const std::vector<Tensor>& E,
                           const std::vector<double>& kappa_t,
                           const std::vector<double>& boundary_face_watts) const {
        return compute(J,E,kappa_t,boundary_face_watts,true);
    }

    // Matrix-free action of the *frozen-E* second-moment operator. Signed
    // Krylov directions are valid; they are NOT treated as physical intensities.
    // Homogeneous boundary currents remove the fixed-wall affine RHS.
    DomSmmFluxV48 linearized_face_action(const std::vector<double>& direction,
                                         const std::vector<Tensor>& E,
                                         const std::vector<double>& kappa_t) const {
        return compute(direction,E,kappa_t,
                       std::vector<double>(mesh_.n_faces(),0.),false);
    }
    std::vector<double> elliptic_action(const std::vector<double>& direction,
                                        const std::vector<Tensor>& E,
                                        const std::vector<double>& kappa_t,
                                        const std::vector<double>& kappa_a) const {
        if(kappa_a.size()!=mesh_.n_cells())
            throw std::invalid_argument("V48 SMM: absorption dimension mismatch");
        auto flux=linearized_face_action(direction,E,kappa_t);
        constexpr double fourpi=12.56637061435917295385;
        for(std::size_t c=0;c<mesh_.n_cells();++c) {
            if(!std::isfinite(kappa_a[c]) || kappa_a[c]<0. ||
               kappa_a[c]>kappa_t[c])
                throw std::invalid_argument("V48 SMM: absorption exceeds extinction");
            flux.divergence_watts[c]+=fourpi*kappa_a[c]*
                geo_.cell_volumes[c]*direction[c];
        }
        return flux.divergence_watts;
    }
private:
    DomSmmFluxV48 compute(const std::vector<double>& J,
                           const std::vector<Tensor>& E,
                           const std::vector<double>& kappa_t,
                           const std::vector<double>& boundary_face_watts,
                           bool positive) const {
        const auto nc=mesh_.n_cells(),nf=mesh_.n_faces();
        if(J.size()!=nc || E.size()!=nc || kappa_t.size()!=nc ||
           boundary_face_watts.size()!=nf)
            throw std::invalid_argument("V48 SMM: input field dimensions mismatch");
        std::array<std::vector<double>,6> P;
        for(auto& p:P)p.resize(nc);
        for(std::size_t c=0;c<nc;++c) {
            if(!std::isfinite(J[c]) || (positive && J[c]<0.) ||
               !(kappa_t[c]>0.) || !std::isfinite(kappa_t[c]))
                throw std::invalid_argument("V48 SMM: invalid moment or extinction");
            for(double value:E[c]) if(!std::isfinite(value))
                throw std::invalid_argument("V48 SMM: nonfinite tensor");
            const double trace=E[c][0]+E[c][1]+E[c][2];
            if(std::abs(trace-1.)>1.e-9 ||
               !v37::tensor_is_spd({E[c][0],E[c][1],E[c][2],
                                      E[c][3],E[c][4],E[c][5]}))
                throw std::invalid_argument("V48 SMM: nonphysical Eddington tensor");
            for(std::size_t k=0;k<6;++k)P[k][c]=E[c][k]*J[c];
        }
        // The WLS gradient acts on P_ij itself, retaining J*grad(E).
        std::array<std::vector<cfdx::core::Vec3>,6> gradP;
        for(auto& g:gradP)g.resize(nc);
        for(std::size_t c=0;c<nc;++c) {
            std::vector<cfdx::core::Vec3> pts;
            pts.reserve(neighbour_cells_[c].size());
            for(auto j:neighbour_cells_[c])pts.push_back(geo_.cell_centres[j]);
            for(std::size_t k=0;k<6;++k) {
                std::vector<double> values;
                values.reserve(neighbour_cells_[c].size());
                for(auto j:neighbour_cells_[c])values.push_back(P[k][j]);
                gradP[k][c]=cfdx::core::weighted_least_squares_gradient_robust(
                    geo_.cell_centres[c],P[k][c],pts,values,
                    cfdx::core::GradientWeighting::INVERSE_DISTANCE_SQUARED,3);
            }
        }
        std::vector<cfdx::core::Vec3> divP(nc);
        for(std::size_t c=0;c<nc;++c) {
            divP[c]={gradP[0][c].x+gradP[3][c].y+gradP[4][c].z,
                     gradP[3][c].x+gradP[1][c].y+gradP[5][c].z,
                     gradP[4][c].x+gradP[5][c].y+gradP[2][c].z};
            if(!finite(divP[c]))throw std::runtime_error("V48 SMM: nonfinite div(EJ)");
        }
        DomSmmFluxV48 out;
        out.face_watts.resize(nf);
        out.divergence_watts.assign(nc,0.);
        constexpr double fourpi=12.56637061435917295385;
        long double internal_owner=0.,internal_neighbour=0.;
        for(std::size_t f=0;f<nf;++f) {
            const auto o=mesh_.ownership().owner(f);
            const auto n=mesh_.ownership().neighbour(f);
            double q;
            if(n>=0) {
                const auto j=static_cast<std::size_t>(n);
                const double extinction=2./(1./kappa_t[o]+1./kappa_t[j]);
                const double do_distance=(geo_.face_centres[f]-geo_.cell_centres[o]).mag();
                const double dn_distance=(geo_.face_centres[f]-geo_.cell_centres[j]).mag();
                const double length=do_distance+dn_distance;
                if(!(length>0.) || !std::isfinite(length))
                    throw std::invalid_argument("V48 SMM: invalid face interpolation");
                const auto grad=(divP[o]*dn_distance+divP[j]*do_distance)*(1./length);
                q=-fourpi/extinction*geo_.face_area_vectors[f].dot(grad);
                // One conservative face current, accumulated on both cells.
                out.divergence_watts[j]-=q;
                internal_owner+=static_cast<long double>(q);
                internal_neighbour-=static_cast<long double>(q);
            } else {
                q=boundary_face_watts[f]; // external HO oracle or prescribed current
            }
            if(!std::isfinite(q))throw std::invalid_argument("V48 SMM: invalid face current");
            out.face_watts[f]=q;
            out.divergence_watts[o]+=q;
        }
        out.internal_cancellation_watts=static_cast<double>(
            internal_owner+internal_neighbour);
        return out;
    }

public:
    // Algebraically conservative HO/LO face correction. Recalculating it
    // whenever E or J changes is MANDATORY; freezing it is not a proof of
    // a converged second-moment method.
    DomSmmConsistencyV48 reconcile(const DomSmmFluxV48& lo,
                                   const std::vector<double>& ho_face_watts) const {
        if(lo.face_watts.size()!=mesh_.n_faces() ||
           lo.divergence_watts.size()!=mesh_.n_cells() ||
           ho_face_watts.size()!=mesh_.n_faces())
            throw std::invalid_argument("V48 SMM: reconciliation dimensions");
        DomSmmConsistencyV48 r;
        r.second_moment=lo;
        r.face_ho_correction.resize(mesh_.n_faces());
        r.reconstructed_ho_current.resize(mesh_.n_faces());
        r.reconstructed_divergence.assign(mesh_.n_cells(),0.);
        std::vector<double> ho_divergence(mesh_.n_cells(),0.);
        for(std::size_t f=0;f<mesh_.n_faces();++f) {
            if(!std::isfinite(ho_face_watts[f]) || !std::isfinite(lo.face_watts[f]))
                throw std::invalid_argument("V48 SMM: nonfinite HO/LO current");
            const double delta=ho_face_watts[f]-lo.face_watts[f];
            const double restored=lo.face_watts[f]+delta;
            r.face_ho_correction[f]=delta;
            r.reconstructed_ho_current[f]=restored;
            r.max_face_identity_relative=std::max(r.max_face_identity_relative,
                std::abs(restored-ho_face_watts[f])/
                std::max({1.,std::abs(restored),std::abs(ho_face_watts[f])}));
            const auto o=mesh_.ownership().owner(f);
            r.reconstructed_divergence[o]+=restored;
            ho_divergence[o]+=ho_face_watts[f];
            const auto n=mesh_.ownership().neighbour(f);
            if(n>=0) {
                const auto j=static_cast<std::size_t>(n);
                r.reconstructed_divergence[j]-=restored;
                ho_divergence[j]-=ho_face_watts[f];
            }
        }
        for(std::size_t c=0;c<mesh_.n_cells();++c)
            r.max_cell_identity_relative=std::max(r.max_cell_identity_relative,
                std::abs(r.reconstructed_divergence[c]-ho_divergence[c])/
                std::max({1.,std::abs(ho_divergence[c]),
                             std::abs(r.reconstructed_divergence[c])}));
        return r;
    }
private:
    static bool finite(const cfdx::core::Vec3& p){
        return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z);
    }
    const cfdx::core::Mesh& mesh_;
    const FvGeometry& geo_;
    std::vector<std::vector<std::size_t>> neighbour_cells_;
};
} // namespace cfdx::physics
