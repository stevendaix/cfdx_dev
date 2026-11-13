#pragma once
// V54 opt-in experimental numerically linear WLS action. Same geometrical
// second-moment FV flux as V48, with stencil weights precomputed once and
// long-double products/face ledger. This is NOT a replacement for the
// historical operator unless equivalence AND both physical gates are shown.
#include "cfdx/physics/dom_smm_gray_wall_v49.h"
#include "cfdx/physics/dom_hilo_arithmetic_v55.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace cfdx::physics {
struct DomV54StableFlux {
    std::vector<double> face_watts,divergence_watts;
};
// Diagnostic extended representation: never implicitly publish 80-bit J
// through CFDX double-valued fields. This is not a production backend.
struct DomV54ExtendedFlux {
    std::vector<long double> face_watts,divergence_watts;
};
class DomStableWlsV54 {
public:
    using Tensor=DomSmmSecondMomentV48::Tensor;
    DomStableWlsV54(const cfdx::core::Mesh& mesh,const FvGeometry& geo,
                    const DomSmmGrayWallV49& wall,
                    std::vector<Tensor> E,std::vector<double> extinction,
                    std::vector<double> absorption)
     :mesh_(mesh),geo_(geo),wall_(wall),E_(std::move(E)),
      extinction_(std::move(extinction)),absorption_(std::move(absorption)) {
        if(!wall_.frozen_operator_matches_v52(
                mesh_,geo_,E_,extinction_,absorption_))
            throw std::invalid_argument("V55 stable WLS: stale/mismatched frozen data");
        const auto n=mesh.n_cells();
        if(n<8 || n!=wall.dimension() || E_.size()!=n ||
           extinction_.size()!=n || absorption_.size()!=n)
            throw std::invalid_argument("V54 stable WLS: field sizes");
        std::vector<std::vector<std::size_t>> direct(n);
        for(std::size_t f=0;f<mesh_.n_faces();++f){
            const auto o=mesh_.ownership().owner(f);
            const auto j=mesh_.ownership().neighbour(f);
            if(j>=0){
                direct[o].push_back(static_cast<std::size_t>(j));
                direct[static_cast<std::size_t>(j)].push_back(o);
            }
        }
        stencils_.resize(n);
        for(std::size_t c=0;c<n;++c){
            auto neighbours=direct[c];
            std::sort(neighbours.begin(),neighbours.end());
            neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
            bool complete=false;
            for(unsigned ring=0;ring<3;++ring){
                if(neighbours.size()>=3){
                    std::vector<cfdx::core::Vec3> pts;
                    std::vector<double> constant;
                    for(const auto j:neighbours){
                        pts.push_back(geo.cell_centres[j]);constant.push_back(1.);
                    }
                    cfdx::core::StencilQuality q;
                    (void)cfdx::core::weighted_least_squares_gradient_robust(
                        geo.cell_centres[c],1.,pts,constant,
                        cfdx::core::GradientWeighting::INVERSE_DISTANCE_SQUARED,0,&q);
                    if(q.rank==3){complete=true;break;}
                }
                if(ring==2)break;
                auto grown=neighbours;
                for(const auto j:neighbours)
                    for(const auto k:direct[j])if(k!=c)grown.push_back(k);
                std::sort(grown.begin(),grown.end());
                grown.erase(std::unique(grown.begin(),grown.end()),grown.end());
                if(grown==neighbours)break;
                neighbours=std::move(grown);
            }
            if(!complete)throw std::invalid_argument("V54 stable WLS: rank deficient");
            std::vector<cfdx::core::Vec3> xyz;
            for(const auto j:neighbours)xyz.push_back(geo.cell_centres[j]);
            for(std::size_t q=0;q<neighbours.size();++q){
                std::vector<double> basis(neighbours.size(),0.);
                basis[q]=1.;
                const auto g=cfdx::core::weighted_least_squares_gradient_robust(
                    geo.cell_centres[c],0.,xyz,basis,
                    cfdx::core::GradientWeighting::INVERSE_DISTANCE_SQUARED,3);
                stencils_[c].push_back({neighbours[q],g});
            }
        }
        for(std::size_t c=0;c<n;++c){
            if(!(extinction_[c]>0.) || !(absorption_[c]>0.) ||
               absorption_[c]>extinction_[c])
                throw std::invalid_argument("V54 stable WLS: invalid optics");
        }
    }
    // Full 80-bit state, gradient differences, fluxes and cell ledger.
    // This has the same mathematical WLS/FV stencil as V48, but it is
    // deliberately a separate research operator: it is not bitwise-equal
    // to a reassembled double-precision production residual.
    DomV54ExtendedFlux evaluate_extended(const std::vector<long double>& J)const{
        const auto nc=mesh_.n_cells(),nf=mesh_.n_faces();
        if(J.size()!=nc)throw std::invalid_argument("V54 extended WLS: moment size");
        for(const auto j:J)if(!std::isfinite(j))
            throw std::invalid_argument("V54 extended WLS: nonfinite input");
        using V=std::array<long double,3>;
        std::vector<V> divP(nc);
        for(std::size_t c=0;c<nc;++c){
            std::array<V,6> gradient{};
            for(const auto& w:stencils_[c]){
                for(std::size_t k=0;k<6;++k){
                    const long double dc=static_cast<long double>(E_[w.cell][k])*J[w.cell]-
                        static_cast<long double>(E_[c][k])*J[c];
                    gradient[k][0]+=static_cast<long double>(w.coeff.x)*dc;
                    gradient[k][1]+=static_cast<long double>(w.coeff.y)*dc;
                    gradient[k][2]+=static_cast<long double>(w.coeff.z)*dc;
                }
            }
            divP[c]={gradient[0][0]+gradient[3][1]+gradient[4][2],
                     gradient[3][0]+gradient[1][1]+gradient[5][2],
                     gradient[4][0]+gradient[5][1]+gradient[2][2]};
        }
        std::vector<long double> cells(nc,0.);
        DomV54ExtendedFlux r;
        r.face_watts.resize(nf);
        constexpr long double fourpi=12.566370614359172953850573533118L;
        for(std::size_t f=0;f<nf;++f){
            const std::size_t o=mesh_.ownership().owner(f);
            const auto nb=mesh_.ownership().neighbour(f);
            long double watts=0.;
            if(nb>=0){
                const std::size_t j=static_cast<std::size_t>(nb);
                const auto& sf=geo_.face_area_vectors[f];
                const auto& fc=geo_.face_centres[f];
                const double a=(fc-geo_.cell_centres[o]).mag();
                const double b=(fc-geo_.cell_centres[j]).mag();
                const long double inv=static_cast<long double>(1./extinction_[o])+
                                      1./extinction_[j];
                const long double factor=-(fourpi/2.)*inv/(a+b);
                const long double dot=static_cast<long double>(sf.x)*
                    (divP[o][0]*b+divP[j][0]*a)+
                    static_cast<long double>(sf.y)*
                    (divP[o][1]*b+divP[j][1]*a)+
                    static_cast<long double>(sf.z)*
                    (divP[o][2]*b+divP[j][2]*a);
                watts=factor*dot;
                cells[j]-=watts;
            }
            cells[o]+=watts;
            r.face_watts[f]=watts;
        }
        for(const auto& wall:wall_.walls())
            cells[wall.owner]+=static_cast<long double>(wall.conductance)*J[wall.owner];
        r.divergence_watts.resize(nc);
        for(std::size_t c=0;c<nc;++c){
            cells[c]+=fourpi*absorption_[c]*geo_.cell_volumes[c]*J[c];
            r.divergence_watts[c]=cells[c];
            if(!std::isfinite(cells[c]))
                throw std::runtime_error("V54 extended WLS: nonfinite flux");
        }
        return r;
    }
    std::vector<long double> apply_extended(const std::vector<long double>& x)const {
        return evaluate_extended(x).divergence_watts;
    }
    // V55 experimental: the SAME precomputed face/stencil geometry as V54,
    // but every E*J product, WLS dot, face power and cell accumulation is
    // evaluated with a portable two-double field. Importantly, this does NOT
    // make the historic V48 double reconstruction bitwise identical.
    struct SplitFlux {
        std::vector<DomHiLoV55> face_watts,divergence_watts;
    };
    SplitFlux evaluate_hilo_v55(const std::vector<DomHiLoV55>& J) const {
        using H=DomHiLoV55;
        using V=std::array<H,3>;
        const auto nc=mesh_.n_cells(), nf=mesh_.n_faces();
        if(J.size()!=nc)
            throw std::invalid_argument("V55 WLS: split moment dimension");
        for(const auto& j:J)
            if(!dom_finite_v55(j))
                throw std::invalid_argument("V55 WLS: nonfinite split state");
        std::vector<V> divP(nc);
        for(std::size_t c=0;c<nc;++c){
            std::array<V,6> grad{};
            for(const auto& w:stencils_[c]){
                for(std::size_t k=0;k<6;++k){
                    const H diff=H(E_[w.cell][k])*J[w.cell]-H(E_[c][k])*J[c];
                    grad[k][0]+=H(w.coeff.x)*diff;
                    grad[k][1]+=H(w.coeff.y)*diff;
                    grad[k][2]+=H(w.coeff.z)*diff;
                }
            }
            divP[c]={grad[0][0]+grad[3][1]+grad[4][2],
                     grad[3][0]+grad[1][1]+grad[5][2],
                     grad[4][0]+grad[5][1]+grad[2][2]};
        }
        SplitFlux result;
        result.face_watts.resize(nf);
        result.divergence_watts.resize(nc);
        const H fourpi=dom_four_pi_v55();
        for(std::size_t f=0;f<nf;++f){
            const std::size_t o=mesh_.ownership().owner(f);
            const auto nb=mesh_.ownership().neighbour(f);
            H watts;
            if(nb>=0){
                const std::size_t j=static_cast<std::size_t>(nb);
                const auto& sf=geo_.face_area_vectors[f];
                const auto& fc=geo_.face_centres[f];
                const double a=(fc-geo_.cell_centres[o]).mag();
                const double b=(fc-geo_.cell_centres[j]).mag();
                if(!(a+b>0.)||!std::isfinite(a+b))
                    throw std::invalid_argument("V55 WLS: invalid face distance");
                const H inv=H(1.)/H(extinction_[o])+H(1.)/H(extinction_[j]);
                const H factor=-(fourpi/H(2.))*inv/H(a+b);
                const H dot=H(sf.x)*(divP[o][0]*H(b)+divP[j][0]*H(a))+
                    H(sf.y)*(divP[o][1]*H(b)+divP[j][1]*H(a))+
                    H(sf.z)*(divP[o][2]*H(b)+divP[j][2]*H(a));
                watts=factor*dot;
                // Exactly the same split face current is used by owner and
                // neighbour; do not separately re-evaluate either side.
                result.divergence_watts[j]-=watts;
            }
            result.divergence_watts[o]+=watts;
            result.face_watts[f]=watts;
        }
        for(const auto& w:wall_.walls())
            result.divergence_watts[w.owner]+=H(w.conductance)*J[w.owner];
        for(std::size_t c=0;c<nc;++c){
            result.divergence_watts[c]+=
                fourpi*H(absorption_[c])*H(geo_.cell_volumes[c])*J[c];
            if(!dom_finite_v55(result.divergence_watts[c]))
                throw std::runtime_error("V55 WLS: nonfinite split face balance");
        }
        return result;
    }
    std::vector<DomHiLoV55> apply_hilo_v55(
        const std::vector<DomHiLoV55>& J)const{
        return evaluate_hilo_v55(J).divergence_watts;
    }
    DomV54StableFlux evaluate(const std::vector<double>& J)const{
        std::vector<long double> x(J.begin(),J.end());
        const auto out=evaluate_extended(x);
        DomV54StableFlux r;
        r.face_watts.reserve(out.face_watts.size());
        r.divergence_watts.reserve(out.divergence_watts.size());
        for(auto q:out.face_watts)r.face_watts.push_back(static_cast<double>(q));
        for(auto q:out.divergence_watts)r.divergence_watts.push_back(static_cast<double>(q));
        return r;
    }
    std::vector<double> apply(const std::vector<double>& x)const {
        return evaluate(x).divergence_watts;
    }

private:
    struct Coeff{std::size_t cell; cfdx::core::Vec3 coeff;};
    const cfdx::core::Mesh& mesh_;
    const FvGeometry& geo_;
    const DomSmmGrayWallV49& wall_;
    std::vector<Tensor> E_;
    std::vector<double> extinction_,absorption_;
    std::vector<std::vector<Coeff>> stencils_;
};
}
