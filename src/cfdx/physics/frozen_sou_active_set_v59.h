#pragma once
// CFDX V59 -- opt-in, read-only frozen generalized SOU face derivative.
// Freeze each active branch ONCE at the physical base state, before Krylov.
// A fixed selection gives a linear face action even at SOU triple ties.
// IMPORTANT: a frozen generalized selection is NOT the two-sided derivative
// of the nonsmooth primal and is NOT a complete 6N Jacobian. V57 physical
// fidelity, true-Krylov and post-step merit vetoes MUST remain in force.
#include "cfdx/physics/convective_face_audit_v56.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

enum class FrozenSouBranchV59 : unsigned char {
    Boundary, Upstream, Downstream, Extrapolated
};

struct FrozenSouFaceV59 {
    std::size_t upstream=0, downstream=0;
    FrozenSouBranchV59 branch=FrozenSouBranchV59::Boundary;
    bool nonunique=false;
    double primal_face_value=0.;
};

struct FrozenSouActiveSetV59 {
    std::vector<FrozenSouFaceV59> faces;
    std::vector<double> base_field, base_flux;
    double high_order_fraction=1.;
    std::size_t nonunique_internal_faces=0;
    std::size_t triple_ties=0;
};

inline FrozenSouActiveSetV59 freeze_sou_active_set_v59(
    const cfdx::core::Mesh& mesh,const FvGeometry& geo,
    const CellScalarV56& field,const FaceScalarV56& phi,
    const ScalarBoundaryConditions& bc,double high_order_fraction=1.,
    const ScalarBoundaryFaceValues* face_values=nullptr,
    const ScalarBoundaryFaceConditions* face_conditions=nullptr) {
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(field.size()!=nc || phi.size()!=nf || field.dimension()!=1 ||
       phi.dimension()!=1 || !(high_order_fraction>=0.) ||
       !(high_order_fraction<=1.) || !std::isfinite(high_order_fraction) ||
       geo.cell_centres.size()!=nc || geo.face_centres.size()!=nf)
        throw std::invalid_argument("V59 invalid frozen SOU mesh, fields, or fraction");
    FrozenSouActiveSetV59 frozen;
    frozen.faces.resize(nf);
    frozen.high_order_fraction=high_order_fraction;
    frozen.base_field.resize(nc);
    frozen.base_flux.resize(nf);
    for(std::size_t c=0;c<nc;++c) {
        if(!std::isfinite(field(c)))
            throw std::invalid_argument("V59 nonfinite frozen velocity");
        frozen.base_field[c]=field(c);
    }
    for(std::size_t f=0;f<nf;++f) {
        if(!std::isfinite(phi(f)))
            throw std::invalid_argument("V59 nonfinite frozen face flux");
        frozen.base_flux[f]=phi(f);
    }
    const auto grad=compute_boundary_aware_wls_gradient(
        field,mesh,geo,bc,face_values,face_conditions);
    for(std::size_t f=0;f<nf;++f) {
        auto& chosen=frozen.faces[f];
        const auto nr=mesh.ownership().neighbour(f);
        if(nr<0) continue;
        const auto o=mesh.ownership().owner(f);
        const auto n=static_cast<std::size_t>(nr);
        chosen.upstream=phi(f)>=0.?o:n;
        chosen.downstream=phi(f)>=0.?n:o;
        const double u=field(chosen.upstream),d=field(chosen.downstream);
        const auto distance=geo.face_centres[f]-geo.cell_centres[chosen.upstream];
        const double inc=grad(chosen.upstream,0)*distance.x+
            grad(chosen.upstream,1)*distance.y+
            grad(chosen.upstream,2)*distance.z;
        const double low=std::min(u,d),high=std::max(u,d),raw=u+inc;
        chosen.primal_face_value=std::clamp(raw,low,high);
        // Explicit, deterministic generalized-Jacobian selection at ties.
        // Do not make this choice depend on the Krylov *direction*.
        if(u==d) {
            chosen.branch=FrozenSouBranchV59::Upstream;
            chosen.nonunique=true;
            ++frozen.triple_ties;
        }else if(raw<=low) {
            chosen.branch=(u<d)?FrozenSouBranchV59::Upstream:
                                   FrozenSouBranchV59::Downstream;
            chosen.nonunique=(raw==low);
        }else if(raw>=high) {
            chosen.branch=(u>d)?FrozenSouBranchV59::Upstream:
                                   FrozenSouBranchV59::Downstream;
            chosen.nonunique=(raw==high);
        }else {
            chosen.branch=FrozenSouBranchV59::Extrapolated;
        }
        chosen.nonunique=chosen.nonunique || phi(f)==0.;
        frozen.nonunique_internal_faces+=chosen.nonunique?1u:0u;
    }
    return frozen;
}

// The direction is evaluated with the BASE-state upwind index, face bounds
// and reconstruction selection. Reapplying V56 is only for the unchanged
// primal flux and boundary-face derivative; EVERY internal SOU tangent is
// replaced with the frozen (direction-independent) branch action.
inline std::vector<FaceFluxTangentV56> apply_frozen_sou_v59(
    const FrozenSouActiveSetV59& frozen,
    const cfdx::core::Mesh& mesh,const FvGeometry& geo,
    const CellScalarV56& field,const CellScalarV56& direction,
    const FaceScalarV56& phi,const FaceScalarV56& dphi,
    const ScalarBoundaryConditions& bc,
    const ScalarBoundaryFaceValues* face_values=nullptr,
    const ScalarBoundaryFaceConditions* face_conditions=nullptr) {
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(frozen.base_field.size()!=nc || frozen.base_flux.size()!=nf ||
       frozen.faces.size()!=nf || field.size()!=nc || direction.size()!=nc ||
       phi.size()!=nf || dphi.size()!=nf ||
       field.dimension()!=1 || direction.dimension()!=1 ||
       phi.dimension()!=1 || dphi.dimension()!=1)
        throw std::invalid_argument("V59 invalid frozen action dimensions");
    for(std::size_t c=0;c<nc;++c)
        if(field(c)!=frozen.base_field[c] || !std::isfinite(direction(c)))
            throw std::invalid_argument("V59 stale base cell state or nonfinite direction");
    for(std::size_t f=0;f<nf;++f)
        if(phi(f)!=frozen.base_flux[f] || !std::isfinite(dphi(f)))
            throw std::invalid_argument("V59 stale base face flux or nonfinite direction");
    auto faces=audit_scalar_convective_faces_v56(mesh,geo,field,direction,
        phi,dphi,ConvectionScheme::SECOND_ORDER_UPWIND,bc,
        face_values,face_conditions,frozen.high_order_fraction);
    const auto homogeneous=homogeneous_boundary_v56(bc);
    const auto hv=homogeneous_face_values_v56(face_values);
    const auto hc=homogeneous_face_conditions_v56(face_conditions,nf);
    const auto dgrad=compute_boundary_aware_wls_gradient(direction,mesh,geo,
        homogeneous,face_values?&hv:nullptr,face_conditions?&hc:nullptr);
    for(std::size_t f=0;f<nf;++f) {
        const auto& branch=frozen.faces[f];
        if(branch.branch==FrozenSouBranchV59::Boundary)continue;
        const auto up=branch.upstream,down=branch.downstream;
        double directional_face=direction(up);
        if(branch.branch==FrozenSouBranchV59::Downstream)
            directional_face=direction(down);
        else if(branch.branch==FrozenSouBranchV59::Extrapolated) {
            const auto offset=geo.face_centres[f]-geo.cell_centres[up];
            directional_face+=dgrad(up,0)*offset.x+
                dgrad(up,1)*offset.y+dgrad(up,2)*offset.z;
        }
        const auto expected=(1.-frozen.high_order_fraction)*field(up)+
            frozen.high_order_fraction*branch.primal_face_value;
        const double tolerance=32.*std::numeric_limits<double>::epsilon()*
            std::max({1.,std::abs(expected),std::abs(faces[f].face_value)});
        if(std::abs(expected-faces[f].face_value)>tolerance)
            throw std::runtime_error("V59 frozen selection disagrees with production primal");
        const double selected=(1.-frozen.high_order_fraction)*direction(up)+
            frozen.high_order_fraction*directional_face;
        faces[f].tangent=dphi(f)*faces[f].face_value+phi(f)*selected;
        faces[f].nonsmooth=branch.nonunique; // KEEP existing Newton veto.
        if(!std::isfinite(faces[f].tangent))
            throw std::runtime_error("V59 nonfinite generalized convective tangent");
    }
    return faces;
}
} // namespace cfdx::physics
