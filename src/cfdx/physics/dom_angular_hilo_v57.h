#pragma once
// CFDX V57: CPU research-only double-double angular FV solve. Actual upwind
// angular CSR is solved in hi+lo BEFORE forming any HO moment/face flux.
// No API of the production DOM solver is switched to this backend.
#include "cfdx/physics/dom_hilo_holo_v56.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

inline DomHiLoV55 dom_abs_v57(DomHiLoV55 x) noexcept {
    return (x.hi<0. || (x.hi==0. && x.lo<0.)) ? -x : x;
}

struct DomV57AngularSolution {
    std::vector<std::vector<DomHiLoV55>> intensities;
    std::vector<DomHiLoV55> moment, face_current;
    std::vector<DomSmmSecondMomentV48::Tensor> eddington;
    std::vector<double> gas_ib;
    double max_true_csr_relative=0.;
    double max_noncanonical_low=0.;
};

// An upwind ordinate is a triangular graph on the tested grids: create a
// topological solve order from the ACTUAL CSR, not from guessed geometric
// sorting. Reject cycles rather than quietly asserting a sweep works on
// arbitrary recirculating/unstructured meshes. Sparse fallback is future work.
class DomAngularHiloTransportV57 {
public:
    explicit DomAngularHiloTransportV57(const DomMomentSchurSystemV44& ho)
       :ho_(ho) {
        const auto n=ho_.n_cells();
        if(n<1)throw std::invalid_argument("V57: empty HO grid");
        orders_.reserve(ho_.directions().size());
        for(std::size_t m=0;m<ho_.directions().size();++m){
            const auto& A=ho_.angular_matrix_research_v57(m);
            if(A.n_rows()!=n || A.n_cols()!=n || !A.is_consistent())
                throw std::invalid_argument("V57: invalid angular CSR");
            std::vector<std::size_t> indegree(n,0);
            std::vector<std::vector<std::size_t>> downstream(n);
            const auto* rows=A.row_offsets_data();
            const auto* cols=A.columns_data();
            const auto* vals=A.values_data();
            for(std::size_t i=0;i<n;++i){
                bool diagonal=false;
                for(auto k=rows[i];k<rows[i+1];++k){
                    const std::size_t j=cols[k];
                    if(j==i){
                        if(!(vals[k]>0.) || !std::isfinite(vals[k]))
                            throw std::invalid_argument("V57: nonpositive angular diagonal");
                        diagonal=true;
                    } else if(vals[k]!=0.) {
                        if(!std::isfinite(vals[k]))
                            throw std::invalid_argument("V57: nonfinite angular stencil");
                        ++indegree[i];
                        downstream[j].push_back(i);
                    }
                }
                if(!diagonal)throw std::invalid_argument("V57: missing angular diagonal");
            }
            std::queue<std::size_t> available;
            for(std::size_t i=0;i<n;++i)if(indegree[i]==0)available.push(i);
            std::vector<std::size_t> order;
            order.reserve(n);
            while(!available.empty()){
                const auto i=available.front();available.pop();
                order.push_back(i);
                for(const auto j:downstream[i])
                    if(--indegree[j]==0)available.push(j);
            }
            if(order.size()!=n)
                throw std::runtime_error("V57: angular CSR contains directed cycles; sweep unqualified");
            orders_.push_back(std::move(order));
        }
    }

    DomV57AngularSolution solve(const std::vector<DomHiLoV55>& x) const {
        using H=DomHiLoV55;
        const auto& mesh=ho_.source_mesh();
        const auto& geo=ho_.source_geometry();
        const auto& dirs=ho_.directions();
        const std::size_t nc=ho_.n_cells(),nf=mesh.n_faces();
        if(x.size()!=ho_.dimension())
            throw std::invalid_argument("V57: moment dimension");
        for(const auto& a:x){
            if(!dom_finite_v55(a) || (a.hi<0.) ||
               (a.hi==0. && a.lo<0.))
                throw std::invalid_argument("V57: invalid split moment");
        }
        DomV57AngularSolution s;
        s.intensities.resize(dirs.size());
        s.moment.assign(x.size(),H(0.));
        s.face_current.assign(nf,H(0.));
        s.eddington.resize(nc);
        s.gas_ib.resize(nc);
        for(std::size_t c=0;c<nc;++c)
            s.gas_ib[c]=blackbody_intensity(ho_.source_temperature()(c));
        std::vector<const DomMomentSchurSystemV44::WallFace*> walls(nf,nullptr);
        for(const auto& w:ho_.wall_faces())walls[w.face]=&w;
        std::vector<std::array<H,6>> second(nc);
        const auto& ka=ho_.cell_absorption();
        const auto& ks=ho_.cell_scattering();
        for(std::size_t m=0;m<dirs.size();++m){
            const auto& d=dirs[m];
            const auto& A=ho_.angular_matrix_research_v57(m);
            const auto* rows=A.row_offsets_data();
            const auto* cols=A.columns_data();
            const auto* vals=A.values_data();
            std::vector<H> rhs(nc),sol(nc);
            for(std::size_t c=0;c<nc;++c)
                rhs[c]=H(geo.cell_volumes[c])*(H(ks[c])*x[c]+H(ka[c])*H(s.gas_ib[c]));
            for(const auto& wall:ho_.wall_faces()){
                const auto& sf=geo.face_area_vectors[wall.face];
                const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                if(mu>=0.)continue;
                const H incoming=H(wall.emissivity)*H(wall.ib)+
                    (wall.reflect_index<ho_.n_reflecting_faces() ?
                     H((1.-wall.emissivity)*(wall.hplus/wall.hminus))*
                        x[nc+wall.reflect_index]:H(0.));
                rhs[wall.owner]+=H(-mu)*incoming;
            }
            for(const auto i:orders_[m]){
                H remainder=rhs[i];
                double diagonal=0.;
                for(auto k=rows[i];k<rows[i+1];++k){
                    if(cols[k]==i)diagonal=vals[k];
                    else remainder-=H(vals[k])*sol[cols[k]];
                }
                sol[i]=remainder/H(diagonal);
                if(!dom_finite_v55(sol[i]) || sol[i].hi<0.)
                    throw std::runtime_error("V57: negative/nonfinite angular intensity");
                // Actually re-evaluate b-A*I on the stored two-double fields.
                H defect=rhs[i];
                for(auto k=rows[i];k<rows[i+1];++k)
                    defect-=H(vals[k])*sol[cols[k]];
                s.max_true_csr_relative=std::max(s.max_true_csr_relative,
                    static_cast<double>(std::abs(dom_to_long_double_v55(defect)) /
                      std::max(1.L,std::abs(dom_to_long_double_v55(rhs[i])))));
                s.max_noncanonical_low=std::max(s.max_noncanonical_low,
                    std::abs(sol[i].lo)/std::max(1.,std::abs(sol[i].hi)));
                const H weighted=H(d.weight)*sol[i];
                s.moment[i]+=weighted/dom_four_pi_v55();
                const H components[6]={H(d.dx)*H(d.dx),H(d.dy)*H(d.dy),
                    H(d.dz)*H(d.dz),H(d.dx)*H(d.dy),
                    H(d.dx)*H(d.dz),H(d.dy)*H(d.dz)};
                for(unsigned k=0;k<6;++k)second[i][k]+=weighted*components[k];
            }
            for(const auto& wall:ho_.wall_faces())
                if(wall.reflect_index<ho_.n_reflecting_faces()){
                    const auto& sf=geo.face_area_vectors[wall.face];
                    const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                    if(mu>0.)s.moment[nc+wall.reflect_index]+=
                        (H(d.weight)*H(mu)/H(wall.area)/H(wall.hplus))*sol[wall.owner];
                }
            for(std::size_t f=0;f<nf;++f){
                const auto& sf=geo.face_area_vectors[f];
                const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
                const std::size_t owner=mesh.ownership().owner(f);
                const auto nb=mesh.ownership().neighbour(f);
                H upwind(0.);
                if(mu>=0.)upwind=sol[owner];
                else if(nb>=0)upwind=sol[static_cast<std::size_t>(nb)];
                else {
                    const auto* w=walls[f];
                    if(!w)throw std::runtime_error("V57: missing wall face");
                    upwind=H(w->emissivity)*H(w->ib)+
                        (w->reflect_index<ho_.n_reflecting_faces()?
                         H((1.-w->emissivity)*(w->hplus/w->hminus))*
                            x[nc+w->reflect_index]:H(0.));
                }
                s.face_current[f]+=H(d.weight)*H(mu)*upwind;
            }
            s.intensities[m]=std::move(sol);
        }
        for(std::size_t c=0;c<nc;++c){
            const auto denom=dom_four_pi_v55()*s.moment[c];
            if(denom.hi<=0.)throw std::runtime_error("V57: nonpositive moment");
            for(unsigned k=0;k<6;++k){
                const auto fraction=second[c][k]/denom;
                s.eddington[c][k]=fraction.hi+fraction.lo;
            }
        }
        return s;
    }
private:
    const DomMomentSchurSystemV44& ho_;
    std::vector<std::vector<std::size_t>> orders_;
};

struct DomV57IndependentAudit {
    bool accepted=false;
    double max_ordinate_equation=0.,max_wall_reflection=0.;
    double max_cell_energy=0.,global_energy=0.,max_fixed_moment=0.;
    double max_face_recompute=0.,csr_residual=0.;
};

// Independent *face*-reassembly: uses geometric directions, signed upwind
// intensities, sources and wall currents. It NEVER calls the CSR matvec,
// corrected LO RHS or HO-LO identity and uses separate summation order.
inline DomV57IndependentAudit audit_dom_angular_hilo_v57(
    const DomMomentSchurSystemV44& ho,
    const std::vector<DomHiLoV55>& x,
    const DomV57AngularSolution& s,
    double tol) {
    using H=DomHiLoV55;
    const auto& mesh=ho.source_mesh();
    const auto& geo=ho.source_geometry();
    const auto& dirs=ho.directions();
    const std::size_t nc=ho.n_cells(),nf=mesh.n_faces();
    if(!(tol>0.) || tol>=1. || x.size()!=ho.dimension() ||
       s.intensities.size()!=dirs.size() || s.face_current.size()!=nf)
        throw std::invalid_argument("V57 audit: inconsistent state");
    std::vector<const DomMomentSchurSystemV44::WallFace*> walls(nf,nullptr);
    for(const auto& w:ho.wall_faces())walls[w.face]=&w;
    std::vector<H> face_current(nf),cell_moment(nc),outgoing(nf);
    DomV57IndependentAudit a;
    a.csr_residual=s.max_true_csr_relative;
    const auto& ka=ho.cell_absorption();
    const auto& ks=ho.cell_scattering();
    for(std::size_t m=0;m<dirs.size();++m){
        const auto& d=dirs[m];
        const auto& I=s.intensities[m];
        if(I.size()!=nc)throw std::invalid_argument("V57 audit: ordinate size");
        std::vector<H> div(nc),absolute(nc);
        for(std::size_t f=0;f<nf;++f){
            const auto& sf=geo.face_area_vectors[f];
            // Each face is evaluated here without reference to cached HO q.
            const H mu=H(d.dx)*H(sf.x)+H(d.dy)*H(sf.y)+H(d.dz)*H(sf.z);
            const std::size_t owner=mesh.ownership().owner(f);
            const auto nb=mesh.ownership().neighbour(f);
            H upwind(0.);
            if(mu.hi>=0.){
                upwind=I[owner];
                if(nb<0)outgoing[f]+=H(d.weight)*mu*upwind;
            }else if(nb>=0)upwind=I[static_cast<std::size_t>(nb)];
            else {
                const auto* w=walls[f];
                if(!w)throw std::runtime_error("V57 audit: missing wall");
                upwind=H(w->emissivity)*H(w->ib)+
                    (w->reflect_index<ho.n_reflecting_faces()?
                     H((1.-w->emissivity)*(w->hplus/w->hminus))*
                         x[nc+w->reflect_index]:H(0.));
            }
            const H q=H(d.weight)*mu*upwind;
            face_current[f]+=q;
            div[owner]+=q;
            absolute[owner]+=dom_abs_v57(q);
            if(nb>=0){
                const auto j=static_cast<std::size_t>(nb);
                div[j]-=q;
                absolute[j]+=dom_abs_v57(q);
            }
        }
        for(std::size_t c=0;c<nc;++c){
            const H emission=H(geo.cell_volumes[c])*H(d.weight)*
                ((H(ka[c])+H(ks[c]))*I[c]-H(ka[c])*H(s.gas_ib[c])-H(ks[c])*x[c]);
            const H r=div[c]+emission;
            const long double norm=std::max({1.L,
                std::abs(dom_to_long_double_v55(absolute[c])),
                std::abs(dom_to_long_double_v55(emission))});
            a.max_ordinate_equation=std::max(a.max_ordinate_equation,
                static_cast<double>(std::abs(dom_to_long_double_v55(r))/norm));
            cell_moment[c]+=H(d.weight)*I[c]/dom_four_pi_v55();
        }
    }
    std::vector<H> cell_flux(nc),abs_flux(nc);
    H boundary_total(0.),vol_total(0.),total_abs(1.);
    for(std::size_t f=0;f<nf;++f){
        const auto owner=mesh.ownership().owner(f);
        const auto nb=mesh.ownership().neighbour(f);
        const H q=face_current[f];
        cell_flux[owner]+=q;
        abs_flux[owner]+=dom_abs_v57(q);
        if(nb>=0){
            const auto j=static_cast<std::size_t>(nb);
            cell_flux[j]-=q;
            abs_flux[j]+=dom_abs_v57(q);
        }else {
            boundary_total+=q;
            total_abs+=dom_abs_v57(q);
            const auto* w=walls[f];
            if(!w)throw std::runtime_error("V57 audit: unclassified boundary");
            if(w->reflect_index<ho.n_reflecting_faces()){
                const H Rout=outgoing[f]/H(w->hplus*w->area);
                const auto residual=x[nc+w->reflect_index]-Rout;
                a.max_wall_reflection=std::max(a.max_wall_reflection,
                    static_cast<double>(std::abs(dom_to_long_double_v55(residual)) /
                    std::max({1.L,
                        std::abs(dom_to_long_double_v55(Rout)),
                        std::abs(dom_to_long_double_v55(x[nc+w->reflect_index]))})));
            }
        }
    }
    for(std::size_t c=0;c<nc;++c){
        const auto moment=cell_moment[c];
        const H emission=dom_four_pi_v55()*H(ka[c])*H(geo.cell_volumes[c])*
            (moment-H(s.gas_ib[c]));
        const H r=cell_flux[c]+emission;
        const long double norm=std::max({1.L,
            std::abs(dom_to_long_double_v55(abs_flux[c])),
            std::abs(dom_to_long_double_v55(emission))});
        a.max_cell_energy=std::max(a.max_cell_energy,
            static_cast<double>(std::abs(dom_to_long_double_v55(r))/norm));
        const auto diff=x[c]-moment;
        a.max_fixed_moment=std::max(a.max_fixed_moment,
            static_cast<double>(std::abs(dom_to_long_double_v55(diff))/
            std::max({1.L,std::abs(dom_to_long_double_v55(x[c])),
                     std::abs(dom_to_long_double_v55(moment))})));
        vol_total+=emission;
        total_abs+=dom_abs_v57(emission);
    }
    a.global_energy=static_cast<double>(
        std::abs(dom_to_long_double_v55(boundary_total+vol_total)) /
        dom_to_long_double_v55(total_abs));
    for(std::size_t f=0;f<nf;++f){
        const auto diff=face_current[f]-s.face_current[f];
        a.max_face_recompute=std::max(a.max_face_recompute,
            static_cast<double>(std::abs(dom_to_long_double_v55(diff)) /
            std::max(1.L,std::abs(dom_to_long_double_v55(face_current[f])))));
    }
    a.accepted=a.max_ordinate_equation<=tol && a.max_wall_reflection<=tol &&
        a.max_cell_energy<=tol && a.global_energy<=tol &&
        a.max_fixed_moment<=tol && a.max_face_recompute<=tol &&
        a.csr_residual<=1e-13;
    return a;
}

} // namespace cfdx::physics
