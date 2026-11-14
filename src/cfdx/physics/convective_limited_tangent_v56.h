#pragma once
// V56: value and active-set directional tangent of the very same limited
// high-order face reconstruction used by scalar FV assembly and flux replay.
// No branch smoothing; nondifferentiable points are reported, not certified.
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

enum class LimitedFaceKindV56 { SOU, TVD_MINMOD, TVD_VANLEER, TVD_VANALBADA, MUSCL_CELL };

struct FaceJetV56 {
    double value = 0.0;
    double tangent = 0.0;
    bool nonsmooth = false;
};

inline FaceJetV56 face_jet_v56(double a, double da = 0.0) {
    if (!std::isfinite(a) || !std::isfinite(da))
        throw std::invalid_argument("V56 non-finite face primal or tangent");
    return {a, da, false};
}
inline FaceJetV56 operator+(FaceJetV56 a, FaceJetV56 b) {
    return {a.value+b.value, a.tangent+b.tangent, a.nonsmooth||b.nonsmooth};
}
inline FaceJetV56 operator-(FaceJetV56 a, FaceJetV56 b) {
    return {a.value-b.value, a.tangent-b.tangent, a.nonsmooth||b.nonsmooth};
}
inline FaceJetV56 operator*(FaceJetV56 a, FaceJetV56 b) {
    return {a.value*b.value, a.tangent*b.value+a.value*b.tangent,
            a.nonsmooth||b.nonsmooth};
}
inline FaceJetV56 operator/(FaceJetV56 a, FaceJetV56 b) {
    if (!(std::abs(b.value)>0.0)) throw std::invalid_argument("V56 zero face divisor");
    return {a.value/b.value,
            (a.tangent*b.value-a.value*b.tangent)/(b.value*b.value),
            a.nonsmooth||b.nonsmooth};
}
inline FaceJetV56 face_abs_v56(FaceJetV56 a) {
    if (a.value<0.0) return {-a.value,-a.tangent,a.nonsmooth};
    return {a.value,a.tangent,a.nonsmooth||(a.value==0.0 && a.tangent!=0.0)};
}
inline FaceJetV56 face_min_v56(FaceJetV56 a,FaceJetV56 b) {
    const bool ambiguous=(a.value==b.value && a.tangent!=b.tangent);
    const FaceJetV56 chosen=a.value<=b.value ? a : b;
    return {std::min(a.value,b.value),chosen.tangent,
            a.nonsmooth||b.nonsmooth||ambiguous};
}
inline FaceJetV56 face_max_v56(FaceJetV56 a,FaceJetV56 b) {
    const bool ambiguous=(a.value==b.value && a.tangent!=b.tangent);
    const FaceJetV56 chosen=a.value>=b.value ? a : b;
    return {std::max(a.value,b.value),chosen.tangent,
            a.nonsmooth||b.nonsmooth||ambiguous};
}
inline FaceJetV56 face_clamp_v56(FaceJetV56 x, FaceJetV56 lower, FaceJetV56 upper) {
    if (lower.value>upper.value) throw std::invalid_argument("V56 inverted face envelope");
    // The primal operation is std::clamp(x,lower,upper); retain its ties.
    return face_min_v56(face_max_v56(x,lower),upper);
}

inline FaceJetV56 face_limiter_v56(FaceJetV56 r,LimitedFaceKindV56 kind) {
    const auto zero=face_jet_v56(0.0), one=face_jet_v56(1.0);
    switch(kind) {
        case LimitedFaceKindV56::TVD_MINMOD:
            return face_max_v56(zero,face_min_v56(one,r));
        case LimitedFaceKindV56::TVD_VANLEER: {
            const auto ar=face_abs_v56(r);
            return face_max_v56(zero,(r+ar)/(one+ar));
        }
        case LimitedFaceKindV56::TVD_VANALBADA:
            return face_max_v56(zero,(r*r+r)/(r*r+one));
        default: throw std::invalid_argument("V56 non-TVD limiter");
    }
}

// `low` and `high` are the exact owner-neighbour or one-ring extrema used by
// production; their tangents carry the selected active sample's derivative.
// `increment` is grad(upwind) dot (Cf-Cu), `full` is grad(upwind) dot(Cd-Cu).
inline FaceJetV56 limited_high_face_v56(
        LimitedFaceKindV56 kind,FaceJetV56 up,FaceJetV56 down,
        FaceJetV56 increment,FaceJetV56 full,
        FaceJetV56 low,FaceJetV56 high) {
    FaceJetV56 result=up;
    if(kind==LimitedFaceKindV56::MUSCL_CELL) {
        auto theta=face_jet_v56(1.0);
        if(increment.value>0.0)
            theta=face_clamp_v56((high-up)/increment,
                                  face_jet_v56(0.0),face_jet_v56(1.0));
        else if(increment.value<0.0)
            theta=face_clamp_v56((low-up)/increment,
                                  face_jet_v56(0.0),face_jet_v56(1.0));
        else if(increment.tangent!=0.0) theta.nonsmooth=true;
        result=up+theta*increment;
    } else if(kind==LimitedFaceKindV56::SOU) {
        result=face_clamp_v56(up+increment,face_min_v56(up,down),
                              face_max_v56(up,down));
    } else {
        const auto delta=down-up;
        FaceJetV56 ratio;
        if(std::abs(delta.value)>1e-14)
            ratio=(face_jet_v56(2.0)*full-delta)/delta;
        else {
            ratio=face_jet_v56(std::abs(full.value)<=1e-14 ? 1.0 : 0.0);
            if (delta.tangent!=0.0 || full.tangent!=0.0) ratio.nonsmooth=true;
        }
        result=face_clamp_v56(up+face_limiter_v56(ratio,kind)*increment,
                              low,high);
    }
    if(!std::isfinite(result.value) || !std::isfinite(result.tangent))
        throw std::runtime_error("V56 non-finite limited face result");
    return result;
}

struct FaceFluxTangentV56 {
    double flux = 0.0;
    double tangent = 0.0;
    double face_value = 0.0;
    bool nonsmooth = false;
};
inline FaceFluxTangentV56 limited_convective_flux_v56(
        FaceJetV56 mass_flux, FaceJetV56 up, FaceJetV56 high,
        double high_order_fraction=1.0) {
    if (!(high_order_fraction>=0.0 && high_order_fraction<=1.0) ||
        !std::isfinite(high_order_fraction))
        throw std::invalid_argument("V56 invalid high-order fraction");
    auto psi=up+face_jet_v56(high_order_fraction)*(high-up);
    auto flux=mass_flux*psi;
    if(!std::isfinite(flux.value) || !std::isfinite(flux.tangent))
        throw std::runtime_error("V56 non-finite convective face flux");
    return {flux.value,flux.tangent,psi.value,flux.nonsmooth ||
        (mass_flux.value==0.0 && mass_flux.tangent!=0.0)};
}
} // namespace cfdx::physics
