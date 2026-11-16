#pragma once
// CFDX V57: bounded directional reconstruction of the actual candidate-only
// six-field operator and its intermediate physical coefficients. This is a
// NUMERICAL chain diagnostic (not a claimed analytic RC/viscosity tangent).
// All stage derivatives are recomputed from X +/- h*v; no accepted-state phi,
// rAU, wall coefficients, nu_t, or production is read or written.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct PhysicalChainTangentV57 {
    std::vector<double> signed_raw;
    std::vector<double> dphi, drau, dmu, dstrain, df2, dwall_mu, dwall_production;
    std::array<double,6> row_halving{};
    double global_halving=0.;
    double global_sided_nonsmoothness=0.;
    std::array<double,6> block_sided{};
    std::size_t flux_sign_crossings=0,wall_active_faces=0;
    bool smooth_both_sides=true;
    std::size_t evaluations=0;
};

inline double v57_norm2(const std::vector<double>& x) {
    long double total=0.;
    for(const auto a:x) {
        if(!std::isfinite(a)) throw std::invalid_argument("V57 nonfinite derivative");
        total+=static_cast<long double>(a)*a;
    }
    return std::sqrt(static_cast<double>(total));
}
inline std::vector<double> v57_difference(const std::vector<double>& a,
                                          const std::vector<double>& b,double h) {
    if(a.size()!=b.size() || !(h>0.) || !std::isfinite(h))
        throw std::invalid_argument("V57 inconsistent stage derivative dimensions");
    std::vector<double> out(a.size());
    for(std::size_t i=0;i<a.size();++i){
        out[i]=(a[i]-b[i])/(2.*h);
        if(!std::isfinite(out[i]))throw std::runtime_error("V57 nonfinite stage derivative");
    }
    return out;
}

template<class Evaluate>
PhysicalChainTangentV57 audit_physical_chain_v57(
    const PhysicalCandidateV51& base,const std::vector<double>& direction,
    Evaluate&& callback,double h,double k_min=0.,double omega_min=0.) {
    validate_physical_candidate_v51(base,k_min,omega_min);
    if(direction.size()!=base.fields.size() || !(h>=1e-9 && h<=1e-3) ||
       !std::isfinite(h) || !(v57_norm2(direction)>0.))
        throw std::invalid_argument("V57 invalid stage tangent direction / step");
    PhysicalChainTangentV57 info;
    auto eval=[&](double a) {
        auto x=base;
        for(std::size_t i=0;i<x.fields.size();++i)
            x.fields[i]+=a*direction[i];
        auto result=evaluate_physical_candidate_v51(x,callback,k_min,omega_min);
        if(!result.stages_v57)
            throw std::runtime_error("V57 evaluator did not reassemble physical stages");
        ++info.evaluations;
        return result;
    };
    // No use of an old normalized merit or a lagged accepted face flux.
    const auto center=eval(0.);
    const auto plus=eval(h),minus=eval(-h);
    const auto smallplus=eval(0.5*h),smallminus=eval(-0.5*h);
    const auto& cp=*center.stages_v57;
    const auto& bp=*smallplus.stages_v57;
    const auto& bm=*smallminus.stages_v57;
    info.signed_raw=v57_difference(smallplus.signed_raw,smallminus.signed_raw,0.5*h);
    info.dphi=v57_difference(smallplus.face_flux,smallminus.face_flux,0.5*h);
    info.drau=v57_difference(bp.inverse_momentum_diagonal,bm.inverse_momentum_diagonal,0.5*h);
    info.dmu=v57_difference(bp.effective_dynamic_viscosity,bm.effective_dynamic_viscosity,0.5*h);
    info.dstrain=v57_difference(bp.strain_rate,bm.strain_rate,0.5*h);
    info.df2=v57_difference(bp.sst_blend_f2,bm.sst_blend_f2,0.5*h);
    info.dwall_mu=v57_difference(bp.wall_momentum_face_mu,bm.wall_momentum_face_mu,0.5*h);
    info.dwall_production=v57_difference(bp.wall_production,bm.wall_production,0.5*h);
    const auto coarse=v57_difference(plus.signed_raw,minus.signed_raw,h);
    std::vector<double> difference(coarse.size());
    for(std::size_t block=0;block<6;++block){
        long double e=0.,denom=0.;
        for(std::size_t i=block*base.cells;i<(block+1)*base.cells;++i){
            const double delta=info.signed_raw[i]-coarse[i];
            difference[i]=delta;
            e+=static_cast<long double>(delta)*delta;
            denom+=static_cast<long double>(info.signed_raw[i])*info.signed_raw[i];
        }
        info.row_halving[block]=std::sqrt(static_cast<double>(e))/
            std::max(std::sqrt(static_cast<double>(denom)),1e-14);
    }
    info.global_halving=v57_norm2(difference)/std::max(v57_norm2(info.signed_raw),1e-14);
    std::vector<double> sided(difference.size());
    for(std::size_t block=0;block<6;++block){
        long double side=0.,scale=0.;
        for(std::size_t i=block*base.cells;i<(block+1)*base.cells;++i){
            // This is (F+ - 2 F0 + F-) / h; it distinguishes an
            // actual branch kink from mere two-step FD consistency.
            sided[i]=(plus.signed_raw[i]-2.*center.signed_raw[i]+
                      minus.signed_raw[i])/h;
            side+=static_cast<long double>(sided[i])*sided[i];
            scale+=static_cast<long double>(info.signed_raw[i])*info.signed_raw[i];
        }
        info.block_sided[block]=std::sqrt(static_cast<double>(side))/
            std::max(std::sqrt(static_cast<double>(scale)),1e-14);
    }
    info.global_sided_nonsmoothness=v57_norm2(sided)/
        std::max(v57_norm2(info.signed_raw),1e-14);
    for(std::size_t f=0;f<plus.face_flux.size();++f)
        if(plus.face_flux[f]*minus.face_flux[f]<0.) ++info.flux_sign_crossings;
    for(double mu:cp.wall_momentum_face_mu)if(mu>0.)++info.wall_active_faces;
    // A two-sided difference can be evaluated across a kink; don't interpret
    // sign crossings as proof of physical non-differentiability by themselves.
    info.smooth_both_sides=(info.global_halving<=0.05 &&
                            info.global_sided_nonsmoothness<=0.05);
    return info;
}
} // namespace cfdx::physics
