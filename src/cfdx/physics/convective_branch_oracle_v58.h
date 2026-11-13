#pragma once
// CFDX V58: read-only, branch-resolved DIFFERENCE audit of the exact V56
// primal convective face flux. This is not a replacement Jacobian. It rejects
// certifying a bilateral derivative where the two one-sided slopes disagree.
// The face mass-flux perturbation is linearized at the candidate; for a full
// six-field derivative, closure of Rhie-Chow/rAU must still be audited by V57.
#include "cfdx/physics/convective_face_audit_v56.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct ConvectiveFaceBranchV58 {
    std::size_t face=0;
    std::size_t owner=0;
    long long neighbour=-1;
    double mass_flux=0.0, dmass_flux=0.0;
    double right=0.0, left=0.0, analytic=0.0;
    double sided_relative=0.0, halving_relative=0.0;
    bool sign_crossing=false;
    bool active_set_ambiguous=false;
    bool effective_kink=false;
    bool unresolved=false;
    bool inert_branch_flag=false;
    bool below_resolution=false;
};
struct ConvectiveBranchEvidenceV58 {
    std::vector<ConvectiveFaceBranchV58> faces;
    std::size_t effective_kinks=0, unresolved=0, inert_branch_flags=0;
    std::size_t flux_crossings=0, ambiguous_active_sets=0;
    std::size_t roundoff_limited=0;
    double global_sided_relative=0.0, global_halving_relative=0.0;
    double maximum_smooth_tangent_error=0.0;
    bool bilateral_certifiable=false;
};


// The exact active *base* inputs of the production reconstruction on a
// selected internal face. In particular, delta=down-up identifies the TVD
// |delta|<=1e-14 fallback; it cannot be inferred from sign(phi) alone.
struct ConvectiveFaceInputsV58 {
    std::size_t up=0,down=0;
    double up_value=0.,down_value=0.,increment=0.,full=0.;
    double lower=0.,upper=0.;
    double reconstructed=0.,central_ratio=0.;
    bool tvd_small_delta_fallback=false;
};

inline ConvectiveFaceInputsV58 inspect_convective_face_v58(
    const cfdx::core::Mesh& mesh,const FvGeometry& geo,
    const CellScalarV56& field,const FaceScalarV56& phi,
    std::size_t face, ConvectionScheme scheme,
    const ScalarBoundaryConditions& bc,
    const ScalarBoundaryFaceValues* face_values=nullptr,
    const ScalarBoundaryFaceConditions* face_conditions=nullptr) {
    if(face>=mesh.n_faces() || mesh.ownership().neighbour(face)<0 ||
       field.size()!=mesh.n_cells() || phi.size()!=mesh.n_faces())
        throw std::invalid_argument("V58 face state needs interior topology");
    ConvectiveFaceInputsV58 info;
    const auto owner=mesh.ownership().owner(face);
    const auto neighbour=static_cast<std::size_t>(mesh.ownership().neighbour(face));
    info.up=phi(face)>=0.0?owner:neighbour;
    info.down=phi(face)>=0.0?neighbour:owner;
    info.up_value=field(info.up);
    info.down_value=field(info.down);
    info.lower=std::min(info.up_value,info.down_value);
    info.upper=std::max(info.up_value,info.down_value);
    info.reconstructed=info.up_value;
    if(scheme==ConvectionScheme::UPWIND) return info;
    const auto grad=compute_boundary_aware_wls_gradient(
        field,mesh,geo,bc,face_values,face_conditions);
    const auto inc=geo.face_centres[face]-geo.cell_centres[info.up];
    const auto full=geo.cell_centres[info.down]-geo.cell_centres[info.up];
    info.increment=grad(info.up,0)*inc.x+grad(info.up,1)*inc.y+
                   grad(info.up,2)*inc.z;
    info.full=grad(info.up,0)*full.x+grad(info.up,1)*full.y+
              grad(info.up,2)*full.z;
    const auto kind=convection_kind_v56(scheme);
    if(scheme==ConvectionScheme::MUSCL_CELL_LIMITED ||
       scheme==ConvectionScheme::TVD ||
       scheme==ConvectionScheme::TVD_VANLEER ||
       scheme==ConvectionScheme::TVD_VANALBADA){
        CellScalarV56 zeros(mesh.n_cells(),"dir_zero_v58","",1);
        zeros.fill(0.0);
        const auto bounds=scalar_stencil_bounds_jets_v56(
            mesh,geo,field,zeros,bc,face_values,face_conditions);
        if(scheme==ConvectionScheme::MUSCL_CELL_LIMITED){
            info.lower=bounds[info.up].low.value;
            info.upper=bounds[info.up].high.value;
        }else{
            const auto edge=geo.cell_centres[neighbour]-geo.cell_centres[owner];
            const auto d2=edge.mag2();
            if(!(d2>0.))throw std::invalid_argument("V58 degenerate face distance");
            const auto offset=geo.face_centres[face]-geo.cell_centres[owner];
            const auto skew=offset-edge*(offset.dot(edge)/d2);
            if(skew.mag2()>1e-24*d2){
                info.lower=std::min(bounds[owner].low.value,bounds[neighbour].low.value);
                info.upper=std::max(bounds[owner].high.value,bounds[neighbour].high.value);
            }
        }
    }
    const auto delta=info.down_value-info.up_value;
    info.tvd_small_delta_fallback=(scheme==ConvectionScheme::TVD ||
        scheme==ConvectionScheme::TVD_VANLEER ||
        scheme==ConvectionScheme::TVD_VANALBADA) &&
        std::abs(delta)<=1e-14;
    if(!info.tvd_small_delta_fallback && delta!=0.)
        info.central_ratio=(2.*info.full-delta)/delta;
    const auto limited=limited_high_face_v56(kind,
        face_jet_v56(info.up_value),face_jet_v56(info.down_value),
        face_jet_v56(info.increment),face_jet_v56(info.full),
        face_jet_v56(info.lower),face_jet_v56(info.upper));
    info.reconstructed=limited.value;
    return info;
}


// Experimental symmetric generalized selection at exact SOU clipping ties.
// It changes ONLY a derivative; the production primal flux remains bitwise
// untouched. For up==down and a nonzero extrapolation, the SOU reconstruction
// clamps to max(up,down) or min(up,down). The mean of their active derivatives
// is a fixed linear selection from the convex generalized derivative set.
// This does not claim a unique two-sided derivative. The caller must still
// compare the selected Jv to central Fphys and the true Arnoldi residual.
inline std::size_t apply_symmetric_sou_ties_v58(
    const cfdx::core::Mesh& mesh,const FvGeometry& geo,
    const CellScalarV56& field,const CellScalarV56& direction,
    const FaceScalarV56& phi,const FaceScalarV56& dphi,
    const ScalarBoundaryConditions& bc,double high_order_fraction,
    std::vector<FaceFluxTangentV56>& faces,
    const ScalarBoundaryFaceValues* face_values=nullptr,
    const ScalarBoundaryFaceConditions* face_conditions=nullptr) {
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    if(field.size()!=nc || direction.size()!=nc || phi.size()!=nf ||
       dphi.size()!=nf || faces.size()!=nf ||
       !(high_order_fraction>=0. && high_order_fraction<=1.) ||
       !std::isfinite(high_order_fraction))
        throw std::invalid_argument("V58 generalized SOU shape/fraction invalid");
    if(high_order_fraction==0.)return 0;
    const auto grad=compute_boundary_aware_wls_gradient(
        field,mesh,geo,bc,face_values,face_conditions);
    std::size_t changed=0;
    for(std::size_t f=0;f<nf;++f){
        const auto nr=mesh.ownership().neighbour(f);
        if(nr<0 || !(std::abs(phi(f))>0.))continue;
        const auto o=mesh.ownership().owner(f),n=static_cast<std::size_t>(nr);
        if(field(o)!=field(n))continue; // no exact active tie
        const auto up=phi(f)>=0.?o:n;
        const auto dr=geo.face_centres[f]-geo.cell_centres[up];
        const double extrap=grad(up,0)*dr.x+grad(up,1)*dr.y+grad(up,2)*dr.z;
        if(!(std::abs(extrap)>32.*std::numeric_limits<double>::epsilon()*
                                    std::max(1.,std::abs(field(up)))))
            continue; // zero gradient has no forced clipping branch
        if(std::abs(faces[f].face_value-field(up))>
           1e-12*std::max(1.,std::abs(field(up))))
            throw std::runtime_error("V58 generalized branch altered primal face flux");
        const double selected=0.5*(direction(o)+direction(n));
        const double mixed=(1.-high_order_fraction)*direction(up)+
                           high_order_fraction*selected;
        faces[f].tangent=dphi(f)*field(up)+phi(f)*mixed;
        if(!std::isfinite(faces[f].tangent))
            throw std::runtime_error("V58 nonfinite generalized tangent");
        // Retain the truth of nonsmoothness. The caller decides whether its
        // separate semismooth+Krylov gate permits using this candidate.
        faces[f].nonsmooth=true;
        ++changed;
    }
    return changed;
}

inline ConvectiveBranchEvidenceV58 audit_convective_branches_v58(
    const cfdx::core::Mesh& mesh,const FvGeometry& geo,
    const CellScalarV56& field,const CellScalarV56& direction,
    const FaceScalarV56& mass_flux,const FaceScalarV56& dmass_flux,
    ConvectionScheme scheme,const ScalarBoundaryConditions& bc,
    double h=1e-6,double kink_relative=1e-2,
    const ScalarBoundaryFaceValues* face_values=nullptr,
    const ScalarBoundaryFaceConditions* face_conditions=nullptr,
    double high_order_fraction=1.0) {
    const auto nf=mesh.n_faces(),nc=mesh.n_cells();
    if(field.size()!=nc || direction.size()!=nc || mass_flux.size()!=nf ||
       dmass_flux.size()!=nf || field.dimension()!=1 || direction.dimension()!=1 ||
       mass_flux.dimension()!=1 || dmass_flux.dimension()!=1 ||
       !(h>=1e-9 && h<=1e-3) || !std::isfinite(h) ||
       !(kink_relative>0.0 && kink_relative<0.5) ||
       !std::isfinite(kink_relative))
        throw std::invalid_argument("V58 invalid state, face shapes or branch-step");
    for(std::size_t c=0;c<nc;++c)
        if(!std::isfinite(field(c)) || !std::isfinite(direction(c)))
            throw std::invalid_argument("V58 nonfinite candidate cell");
    for(std::size_t f=0;f<nf;++f)
        if(!std::isfinite(mass_flux(f)) || !std::isfinite(dmass_flux(f)))
            throw std::invalid_argument("V58 nonfinite candidate flux");
    const auto base=audit_scalar_convective_faces_v56(mesh,geo,field,direction,
        mass_flux,dmass_flux,scheme,bc,face_values,face_conditions,
        high_order_fraction);
    const auto primal=[&](double step){
        auto x=field; auto phi=mass_flux;
        for(std::size_t c=0;c<nc;++c)x(c)+=step*direction(c);
        for(std::size_t f=0;f<nf;++f)phi(f)+=step*dmass_flux(f);
        // Diffusion is explicitly zero: only the production primal convective
        // reconstruction is returned (same boundary conditions and limiter).
        return reconstruct_scalar_transport_flux(mesh,geo,phi,x,0.,bc,true,
            scheme,face_conditions,nullptr,face_values,nullptr,
            high_order_fraction);
    };
    const auto fp=primal(h),fm=primal(-h);
    const auto hp=primal(0.5*h),hm=primal(-0.5*h);
    ConvectiveBranchEvidenceV58 out;
    out.faces.reserve(nf);
    long double sided_sq=0.0,halving_sq=0.0,reference_sq=0.0;
    for(std::size_t f=0;f<nf;++f){
        ConvectiveFaceBranchV58 v;
        v.face=f;
        v.owner=mesh.ownership().owner(f);
        v.neighbour=mesh.ownership().neighbour(f);
        v.mass_flux=mass_flux(f);v.dmass_flux=dmass_flux(f);
        v.analytic=base[f].tangent;
        const double r=(fp(f)-base[f].flux)/h;
        const double l=(base[f].flux-fm(f))/h;
        v.right=2.*(hp(f)-base[f].flux)/h;
        v.left=2.*(base[f].flux-hm(f))/h;
        const double scale=std::max({std::abs(v.right),std::abs(v.left),1e-12});
        v.sided_relative=std::abs(v.right-v.left)/scale;
        v.halving_relative=std::max(std::abs(r-v.right),std::abs(l-v.left))/scale;
        v.active_set_ambiguous=base[f].nonsmooth;
        v.sign_crossing=(mass_flux(f)-h*dmass_flux(f))*
                         (mass_flux(f)+h*dmass_flux(f))<0.0;
        // A genuine kink must be reproducible on halving. If slopes are
        // still moving, the audit is inconclusive rather than "smooth".
        // Subtracting two nearly identical primal face fluxes amplifies
        // absolute roundoff by 1/h. A large *local relative* mismatch on an
        // almost-zero derivative must not be called a physical kink. The
        // floor is grounded in the floating point precision and the primal
        // flux scale, not in a relaxed CFD convergence tolerance.
        const double primal_scale=std::max({std::abs(base[f].flux),
                                            std::abs(fp(f)),std::abs(fm(f)),
                                            std::abs(hp(f)),std::abs(hm(f))});
        const double resolution=128.*std::numeric_limits<double>::epsilon()*
                                primal_scale/h+1e-13;
        const double h_abs=std::max(std::abs(r-v.right),std::abs(l-v.left));
        const double s_abs=std::abs(v.right-v.left);
        v.unresolved=h_abs>kink_relative*scale && h_abs>resolution;
        v.effective_kink=!v.unresolved && s_abs>kink_relative*scale &&
                         s_abs>resolution;
        v.below_resolution=(h_abs>kink_relative*scale ||
                            s_abs>kink_relative*scale) &&
                           h_abs<=resolution && s_abs<=resolution;
        v.inert_branch_flag=v.active_set_ambiguous && !v.effective_kink &&
                            !v.unresolved;
        if(v.active_set_ambiguous)++out.ambiguous_active_sets;
        if(v.sign_crossing)++out.flux_crossings;
        if(v.effective_kink)++out.effective_kinks;
        if(v.unresolved)++out.unresolved;
        if(v.inert_branch_flag)++out.inert_branch_flags;
        if(v.below_resolution)++out.roundoff_limited;
        if(!v.effective_kink && !v.unresolved)
            out.maximum_smooth_tangent_error=std::max(
                out.maximum_smooth_tangent_error,
                std::abs(v.analytic-0.5*(v.right+v.left)));
        sided_sq+=static_cast<long double>(v.right-v.left)*(v.right-v.left);
        halving_sq+=static_cast<long double>(v.right-r)*(v.right-r)+
                    static_cast<long double>(v.left-l)*(v.left-l);
        reference_sq+=static_cast<long double>(v.right)*v.right+
                      static_cast<long double>(v.left)*v.left;
        out.faces.push_back(v);
    }
    const double denom=std::max(std::sqrt(static_cast<double>(reference_sq)),1e-12);
    out.global_sided_relative=std::sqrt(static_cast<double>(sided_sq))/denom;
    out.global_halving_relative=std::sqrt(static_cast<double>(halving_sq))/denom;
    out.bilateral_certifiable=(out.effective_kinks==0 && out.unresolved==0 &&
        out.global_sided_relative<=kink_relative &&
        out.global_halving_relative<=kink_relative);
    return out;
}
} // namespace cfdx::physics
