#pragma once
// Read-only V56 differential audit on the *production* cell/face topology.
// The reconstruction and limiter value function is also used by the actual
// FV assembly and flux replay. This audit differentiates the convective part
// F_face=phi_f*psi_f (not diffusion, Rhie-Chow, or the full six-field Fphys).
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

using CellScalarV56=cfdx::core::Field<double,cfdx::core::Location::CELL>;
using FaceScalarV56=cfdx::core::Field<double,cfdx::core::Location::FACE>;

struct CellBoundsJetV56 { FaceJetV56 low,high; };

inline ScalarBoundaryConditions homogeneous_boundary_v56(const ScalarBoundaryConditions& bc) {
    auto out=bc;
    for(auto& kv:out){kv.second.value=0.0;kv.second.gradient=0.0;}
    return out;
}
inline ScalarBoundaryFaceValues homogeneous_face_values_v56(const ScalarBoundaryFaceValues* bc) {
    ScalarBoundaryFaceValues out;
    if(!bc) return out;
    out=*bc;
    // Non-finite unused positions must remain unselected by the resolver.
    for(auto& kv:out.values) for(auto& v:kv.second)
        if(std::isfinite(v)) v=0.0;
    return out;
}
inline ScalarBoundaryFaceConditions homogeneous_face_conditions_v56(
        const ScalarBoundaryFaceConditions* bc,std::size_t nf) {
    ScalarBoundaryFaceConditions out(nf);
    if(bc){out=*bc;for(auto& item:out.conditions)item.value=0.0;}
    return out;
}

inline std::vector<CellBoundsJetV56> scalar_stencil_bounds_jets_v56(
        const cfdx::core::Mesh& mesh,const FvGeometry& geo,
        const CellScalarV56& field,const CellScalarV56& direction,
        const ScalarBoundaryConditions& bcs,
        const ScalarBoundaryFaceValues* face_values=nullptr,
        const ScalarBoundaryFaceConditions* face_conditions=nullptr) {
    const auto nc=mesh.n_cells();
    if(field.size()!=nc || direction.size()!=nc ||
       field.dimension()!=1 || direction.dimension()!=1)
        throw std::invalid_argument("V56 stencil field shape mismatch");
    std::vector<CellBoundsJetV56> result(nc);
    for(std::size_t c=0;c<nc;++c){
        auto lo=face_jet_v56(field(c),direction(c)),hi=lo;
        for(auto k=mesh.cells().offsets_data()[c];k<mesh.cells().offsets_data()[c+1];++k){
            const auto f=static_cast<std::size_t>(mesh.cells().faces_data()[k]);
            const auto nr=mesh.ownership().neighbour(f);
            FaceJetV56 sample;
            if(nr>=0){
                const auto other=mesh.ownership().owner(f)==c
                    ? static_cast<std::size_t>(nr) : mesh.ownership().owner(f);
                sample=face_jet_v56(field(other),direction(other));
            } else {
                const auto patch=geo.face_patch[f];
                if(patch<mesh.boundary().n_patches() &&
                   mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
                    continue;
                const auto bc=resolve_scalar_boundary_condition(mesh,geo,f,bcs,
                    face_values,face_conditions);
                if(bc.type==ScalarBoundaryType::FIXED_VALUE)
                    sample=face_jet_v56(bc.value,0.0);
                else if(bc.type==ScalarBoundaryType::FIXED_GRADIENT)
                    sample=face_jet_v56(field(c)+bc.gradient*boundary_normal_distance(geo,f,c),
                                        direction(c));
                else sample=face_jet_v56(field(c),direction(c));
            }
            lo=face_min_v56(lo,sample);
            hi=face_max_v56(hi,sample);
        }
        result[c]={lo,hi};
    }
    return result;
}

inline std::vector<FaceFluxTangentV56> audit_scalar_convective_faces_v56(
        const cfdx::core::Mesh& mesh,const FvGeometry& geo,
        const CellScalarV56& field,const CellScalarV56& direction,
        const FaceScalarV56& face_flux,const FaceScalarV56& flux_direction,
        ConvectionScheme scheme,const ScalarBoundaryConditions& bcs,
        const ScalarBoundaryFaceValues* face_values=nullptr,
        const ScalarBoundaryFaceConditions* face_conditions=nullptr,
        double high_order_fraction=1.0) {
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(field.size()!=nc || direction.size()!=nc || face_flux.size()!=nf ||
       flux_direction.size()!=nf || field.dimension()!=1 || direction.dimension()!=1 ||
       face_flux.dimension()!=1 || flux_direction.dimension()!=1 ||
       geo.cell_centres.size()!=nc || geo.face_centres.size()!=nf ||
       geo.face_patch.size()!=nf || !(high_order_fraction>=0.0 && high_order_fraction<=1.0))
        throw std::invalid_argument("V56 face tangent shapes / homotopy invalid");
    if(scheme!=ConvectionScheme::UPWIND) (void)convection_kind_v56(scheme);
    CellScalarV56 grad,dgrad;
    const bool high=scheme!=ConvectionScheme::UPWIND;
    if(high){
        grad=compute_boundary_aware_wls_gradient(field,mesh,geo,bcs,
            face_values,face_conditions);
        const auto homogeneous=homogeneous_boundary_v56(bcs);
        const auto dhfv=homogeneous_face_values_v56(face_values);
        const auto dhfc=homogeneous_face_conditions_v56(face_conditions,nf);
        dgrad=compute_boundary_aware_wls_gradient(direction,mesh,geo,
            homogeneous,face_values?&dhfv:nullptr,face_conditions?&dhfc:nullptr);
    }
    const auto bounds=high?scalar_stencil_bounds_jets_v56(
        mesh,geo,field,direction,bcs,face_values,face_conditions)
        :std::vector<CellBoundsJetV56>{};
    std::vector<FaceFluxTangentV56> out(nf);
    for(std::size_t f=0;f<nf;++f){
        const auto o=mesh.ownership().owner(f);
        const auto nr=mesh.ownership().neighbour(f);
        const auto F=face_jet_v56(face_flux(f),flux_direction(f));
        if(nr<0){
            const auto patch=geo.face_patch[f];
            if(patch<mesh.boundary().n_patches() &&
               mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
                continue;
            const auto bc=resolve_scalar_boundary_condition(mesh,geo,f,bcs,
                face_values,face_conditions);
            auto psi=face_jet_v56(field(o),direction(o));
            if(bc.type==ScalarBoundaryType::FIXED_VALUE && F.value<0.0)
                psi=face_jet_v56(bc.value);
            if(bc.type==ScalarBoundaryType::FIXED_GRADIENT)
                psi=psi+face_jet_v56(bc.gradient*boundary_normal_distance(geo,f,o));
            out[f]=limited_convective_flux_v56(F,psi,psi,1.0);
            continue;
        }
        const auto n=static_cast<std::size_t>(nr);
        const bool owner_up=F.value>=0.0;
        const auto up=owner_up?o:n,down=owner_up?n:o;
        const auto u=face_jet_v56(field(up),direction(up));
        const auto d=face_jet_v56(field(down),direction(down));
        auto reconstructed=u;
        if(high){
            const auto& cu=geo.cell_centres[up];
            const auto& cf=geo.face_centres[f];
            const auto& cd=geo.cell_centres[down];
            const auto ginc=cf-cu,gfull=cd-cu;
            const auto value_increment=grad(up,0)*ginc.x+grad(up,1)*ginc.y+grad(up,2)*ginc.z;
            const auto directional_increment=dgrad(up,0)*ginc.x+dgrad(up,1)*ginc.y+dgrad(up,2)*ginc.z;
            const auto value_full=grad(up,0)*gfull.x+grad(up,1)*gfull.y+grad(up,2)*gfull.z;
            const auto directional_full=dgrad(up,0)*gfull.x+dgrad(up,1)*gfull.y+dgrad(up,2)*gfull.z;
            auto low=face_min_v56(u,d),high_bound=face_max_v56(u,d);
            if(scheme==ConvectionScheme::MUSCL_CELL_LIMITED){
                low=bounds[up].low;high_bound=bounds[up].high;
            } else if(scheme==ConvectionScheme::TVD ||
                    scheme==ConvectionScheme::TVD_VANLEER ||
                    scheme==ConvectionScheme::TVD_VANALBADA){
                const auto edge=geo.cell_centres[n]-geo.cell_centres[o];
                const auto d2=edge.mag2();
                if(!(d2>0.0))throw std::invalid_argument("V56 degenerate face tangent geometry");
                const auto offset=geo.face_centres[f]-geo.cell_centres[o];
                const auto skew=offset-edge*(offset.dot(edge)/d2);
                if(skew.mag2()>1.e-24*d2){
                    low=face_min_v56(bounds[o].low,bounds[n].low);
                    high_bound=face_max_v56(bounds[o].high,bounds[n].high);
                }
            }
            reconstructed=limited_high_face_v56(convection_kind_v56(scheme),u,d,
                face_jet_v56(value_increment,directional_increment),
                face_jet_v56(value_full,directional_full),low,high_bound);
        }
        out[f]=limited_convective_flux_v56(F,u,reconstructed,high_order_fraction);
    }
    return out;
}

// The bounded steady equation contains -div(phi)*psi. Its derivative is
// essential when face flux changes with the state; a face-only tangent is not
// enough for the assembled momentum/transport residual.
struct CellConvectiveResidualV56 {
    std::vector<double> value, tangent;
    std::size_t nonsmooth_faces=0;
};
inline CellConvectiveResidualV56 scatter_bounded_convection_v56(
        const cfdx::core::Mesh& mesh,const FvGeometry& geo,
        const CellScalarV56& field,const CellScalarV56& direction,
        const FaceScalarV56& phi,const FaceScalarV56& dphi,
        const std::vector<FaceFluxTangentV56>& faces,
        bool bounded=true) {
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(field.size()!=nc || direction.size()!=nc || phi.size()!=nf ||
       dphi.size()!=nf || faces.size()!=nf)
        throw std::invalid_argument("V56 bounded residual shape mismatch");
    CellConvectiveResidualV56 out;
    out.value.assign(nc,0.0);out.tangent.assign(nc,0.0);
    std::vector<double> div(nc,0.0),ddiv(nc,0.0);
    for(std::size_t f=0;f<nf;++f){
        const auto o=mesh.ownership().owner(f);
        const auto nr=mesh.ownership().neighbour(f);
        if(nr<0){
            const auto patch=geo.face_patch[f];
            if(patch<mesh.boundary().n_patches() &&
               mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
                continue;
        }
        if(faces[f].nonsmooth)++out.nonsmooth_faces;
        out.value[o]+=faces[f].flux;
        out.tangent[o]+=faces[f].tangent;
        div[o]+=phi(f);ddiv[o]+=dphi(f);
        if(nr>=0){
            const auto n=static_cast<std::size_t>(nr);
            out.value[n]-=faces[f].flux;
            out.tangent[n]-=faces[f].tangent;
            div[n]-=phi(f);ddiv[n]-=dphi(f);
        }
    }
    if(bounded) for(std::size_t c=0;c<nc;++c){
        out.value[c]-=div[c]*field(c);
        out.tangent[c]-=ddiv[c]*field(c)+div[c]*direction(c);
    }
    return out;
}
} // namespace cfdx::physics
