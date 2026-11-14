#pragma once
// V57: additive split of the true candidate-only FV residual, not a modified
// physical equation. Differentiates the SAME production V56 limited convection
// analytically on its active branch, and differences the REMAINDER, including
// rebuilt RC/rAU, turbulent/viscous stress and SST. Requires two fresh F(X).
// If the active convective branch is not locally consistent with those two
// fresh evaluations, this diagnostic refuses certification.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct ConvectiveTangentV57 {
    std::vector<double> signed_raw; // 3N, exact active-set convection tangent
    std::size_t nonsmooth_faces=0;
};
struct HybridPhysicalJvV57 {
    std::vector<double> action;
    double convective_relative_discrepancy=0.;
    std::size_t nonsmooth_faces=0;
};

using GlobalPhysicalConvectiveTangentV57=std::function<ConvectiveTangentV57(
    const PhysicalCandidateV51&,const std::vector<double>&,
    const std::vector<double>&,const std::vector<double>&)>;

inline HybridPhysicalJvV57 combine_physical_jv_v57(
    const PhysicalCandidateV51& base,const PhysicalResidualV51& base_result,
    const PhysicalResidualV51& plus,const PhysicalResidualV51& minus,
    const std::vector<double>& physical_direction,
    const std::array<double,6>& state_scale,
    const std::array<double,6>& row_scale,
    double h,const GlobalPhysicalConvectiveTangentV57& convective,
    bool pin_pressure,std::size_t pressure_cell) {
    const std::size_t n=base.cells,len=base.fields.size();
    if(!convective || !base_result.stages_v57 || !plus.stages_v57 ||
       !minus.stages_v57 || physical_direction.size()!=len || len!=6*n ||
       !(h>0.) || !std::isfinite(h) || pressure_cell>=n ||
       plus.signed_raw.size()!=len || minus.signed_raw.size()!=len ||
       plus.face_flux.size()!=base_result.face_flux.size() ||
       minus.face_flux.size()!=base_result.face_flux.size())
        throw std::invalid_argument("V57 hybrid Jv missing candidate-owned stage");
    for(std::size_t b=0;b<6;++b)
        if(!(row_scale[b]>0.) || !(state_scale[b]>0.) ||
           !std::isfinite(row_scale[b]) || !std::isfinite(state_scale[b]))
            throw std::invalid_argument("V57 invalid frozen state/row scaling");
    std::vector<double> dphi(plus.face_flux.size());
    for(std::size_t f=0;f<dphi.size();++f)
        dphi[f]=(plus.face_flux[f]-minus.face_flux[f])/(2.*h);
    auto analytic=convective(base,physical_direction,base_result.face_flux,dphi);
    if(analytic.signed_raw.size()!=3*n)
        throw std::runtime_error("V57 incomplete convective tangent");
    const auto& pc=plus.stages_v57->momentum_convection;
    const auto& mc=minus.stages_v57->momentum_convection;
    if(pc.size()!=3*n || mc.size()!=3*n)
        throw std::runtime_error("V57 incomplete rebuilt convective primal");
    long double e=0.,norm=0.;
    HybridPhysicalJvV57 result;
    result.action.resize(len);
    result.nonsmooth_faces=analytic.nonsmooth_faces;
    for(std::size_t i=0;i<len;++i) {
        double deriv=(plus.signed_raw[i]-minus.signed_raw[i])/(2.*h);
        if(i<3*n){
            const double numerical=(pc[i]-mc[i])/(2.*h);
            const double diff=analytic.signed_raw[i]-numerical;
            e+=static_cast<long double>(diff)*diff;
            norm+=static_cast<long double>(numerical)*numerical;
            // Correctly differentiates the residual remainder F-C, not F
            // itself, otherwise convection would be counted twice.
            deriv+=diff;
        }
        result.action[i]=deriv/row_scale[i/n];
        if(!std::isfinite(result.action[i]))
            throw std::runtime_error("V57 nonfinite hybrid Jacobian action");
    }
    if(pin_pressure)
        result.action[3*n+pressure_cell]=
            physical_direction[3*n+pressure_cell]/state_scale[3];
    result.convective_relative_discrepancy=std::sqrt(static_cast<double>(e))/
        std::max(std::sqrt(static_cast<double>(norm)),1e-12);
    if(!std::isfinite(result.convective_relative_discrepancy))
        throw std::runtime_error("V57 invalid analytic/numeric convective split");
    return result;
}
} // namespace cfdx::physics
