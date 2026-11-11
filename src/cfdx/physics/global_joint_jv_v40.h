#pragma once
// V40 contract for a future production full-FV Newton-Krylov operator.
// The supplied evaluator MUST recompute conservative face flux (including
// Rhie-Chow), momentum, viscosity and its SST dependence, wall terms, and
// k/omega rows at ONE supplied state. This utility cannot certify those
// physics unless the evaluator is wired to the real production assembler.
// No state is ever modified, and an active min/max branch is identified by
// inconsistency of one-sided derivatives rather than regularized silently.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
using GlobalJointResidualEvaluatorV40 =
    std::function<std::vector<double>(const std::vector<double>&)>;
struct GlobalJointJvV40 {
    std::vector<double> scaled_jv;
    double max_relative_one_sided_gap = std::numeric_limits<double>::infinity();
    double normalized_step = 0.0;
    bool smooth = false;
    bool residual_finite = false;
    std::size_t function_evaluations = 0;
};

// Components are contiguous [Ux,Uy,Uz,p,k,omega], each of length n.
// The perturbation uses affine U/p and logarithmic k/omega. The direction
// belongs to normalized coordinates (U/p scaled by variable_steps and
// k/omega in logarithmic units). Pressure gauge can be pinned explicitly.
inline GlobalJointJvV40 audit_joint_jv_v40(
    const std::vector<double>& physical,
    const std::vector<double>& direction,
    const std::vector<double>& variable_steps,
    const std::vector<double>& residual_scales,
    std::size_t n,
    const GlobalJointResidualEvaluatorV40& evaluate,
    double k_floor,
    double omega_floor,
    std::size_t pressure_gauge = std::numeric_limits<std::size_t>::max(),
    double branch_disagreement_max = 0.03)
{
    if (!n || n>std::numeric_limits<std::size_t>::max()/6 ||
        physical.size()!=6*n || direction.size()!=6*n ||
        variable_steps.size()!=6*n || residual_scales.size()!=6*n || !evaluate ||
        !std::isfinite(k_floor) || !(k_floor>0.0) ||
        !std::isfinite(omega_floor) || !(omega_floor>0.0) ||
        !std::isfinite(branch_disagreement_max) ||
        !(branch_disagreement_max>0.0 && branch_disagreement_max<1.0))
        throw std::invalid_argument("V40 global Jv: incomplete physical residual contract");
    if (pressure_gauge<n && direction[3*n+pressure_gauge]!=0.0)
        throw std::invalid_argument("V40 global Jv: pressure reference must not drift");
    double max_direction=0.0;
    for (std::size_t i=0;i<6*n;++i) {
        if (!std::isfinite(physical[i]) || !std::isfinite(direction[i]) ||
            !(variable_steps[i]>0.0) || !std::isfinite(variable_steps[i]) ||
            !(residual_scales[i]>0.0) || !std::isfinite(residual_scales[i]))
            throw std::invalid_argument("V40 global Jv: nonfinite state/scale");
        if (i>=4*n && physical[i] <=
            (i<5*n?k_floor:omega_floor)*(1.0+1e-4))
            throw std::invalid_argument("V40 global Jv: active SST positivity floor");
        max_direction=std::max(max_direction,std::abs(direction[i]));
    }
    if (!(max_direction>0.0)) throw std::invalid_argument("V40 global Jv: zero direction");
    const double h=1e-5/max_direction;
    auto perturb=[&](double sign) {
        auto x=physical;
        for (std::size_t i=0;i<6*n;++i) {
            if (i<4*n) x[i]+=sign*h*direction[i]*variable_steps[i];
            else x[i]*=std::exp(sign*h*direction[i]);
            if (!std::isfinite(x[i]) || (i>=4*n &&
                x[i] <= (i<5*n?k_floor:omega_floor)))
                throw std::domain_error("V40 global Jv: out-of-domain trial");
        }
        return x;
    };
    const auto sample=[&](const std::vector<double>& x) {
        auto v=evaluate(x);
        if (v.size()!=6*n)
            throw std::runtime_error("V40 global Jv: missing physical residual rows");
        for (double z:v)
            if (!std::isfinite(z))
                throw std::runtime_error("V40 global Jv: nonfinite reassembled FV residual");
        return v;
    };
    GlobalJointJvV40 out;
    const auto f=sample(physical);++out.function_evaluations;
    const auto fp=sample(perturb(1.0));++out.function_evaluations;
    const auto fm=sample(perturb(-1.0));++out.function_evaluations;
    out.scaled_jv.resize(6*n);
    double max_gap=0.0;
    for (std::size_t i=0;i<6*n;++i) {
        const double dp=(fp[i]-f[i])/(h*residual_scales[i]);
        const double dm=(f[i]-fm[i])/(h*residual_scales[i]);
        const double dc=(fp[i]-fm[i])/(2*h*residual_scales[i]);
        const double gap=std::abs(dp-dm)/std::max({1e-10,std::abs(dc),std::abs(dp),std::abs(dm)});
        if (!std::isfinite(dc) || !std::isfinite(gap))
            throw std::runtime_error("V40 global Jv: nonfinite directional derivative");
        out.scaled_jv[i]=dc;
        max_gap=std::max(max_gap,gap);
    }
    out.normalized_step=h;
    out.max_relative_one_sided_gap=max_gap;
    out.smooth=max_gap<=branch_disagreement_max;
    out.residual_finite=true;
    return out;
}
} // namespace cfdx::physics
