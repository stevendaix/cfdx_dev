#pragma once
// V55 exploratory arithmetic gate. V48 historic apply(double) is the
// authoritative algebraic double-state oracle. No HO/LO or physical gate is
// replaced, and NO state is published by this function.
#include "cfdx/physics/dom_wls_precision_audit_v54.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct DomWlsHiloGateV55 {
    double hilo_l2=std::numeric_limits<double>::infinity();
    double extended_l2=std::numeric_limits<double>::infinity();
    double historic_double_l2=std::numeric_limits<double>::infinity();
    double absolute_face_error_watts=std::numeric_limits<double>::infinity();
    double max_apply_difference=std::numeric_limits<double>::infinity();
    bool hilo_pass=false,original_double_pass=false;
    bool publish=false; // remains false: physical HO and energy gates absent
};
inline DomWlsHiloGateV55 audit_wls_hilo_v55(
        const DomStableWlsV54& stable,
        const DomSmmGrayWallV49& historic,
        const std::vector<double>& rhs,
        const std::vector<DomHiLoV55>& diagnostic_state,
        double strict_tolerance=1e-13) {
    if(diagnostic_state.size()!=rhs.size() || rhs.size()!=historic.dimension() ||
       !(strict_tolerance>0.) || !std::isfinite(strict_tolerance))
        throw std::invalid_argument("V55 gate: invalid input");
    std::vector<double> rounded(diagnostic_state.size());
    std::vector<long double> extended(diagnostic_state.size());
    for(std::size_t i=0;i<rhs.size();++i){
        if(!std::isfinite(rhs[i]) || !dom_finite_v55(diagnostic_state[i]))
            throw std::invalid_argument("V55 gate: nonfinite state/RHS");
        rounded[i]=diagnostic_state[i].hi;
        extended[i]=dom_to_long_double_v55(diagnostic_state[i]);
    }
    const auto split=stable.evaluate_hilo_v55(diagnostic_state);
    const auto long_double=stable.evaluate_extended(extended);
    const auto original=historic.apply(rounded);
    long double split_sq=0.,extended_sq=0.,historic_sq=0.,den=0.;
    long double faces=0.,max_apply=0.;
    for(std::size_t i=0;i<rhs.size();++i){
        const long double b=rhs[i];
        const long double s=dom_to_long_double_v55(split.divergence_watts[i]);
        const long double e=long_double.divergence_watts[i];
        const long double h=original[i];
        if(!std::isfinite(s) || !std::isfinite(e) || !std::isfinite(h))
            throw std::runtime_error("V55 gate: nonfinite evaluated residual");
        split_sq+=(s-b)*(s-b);
        extended_sq+=(e-b)*(e-b);
        historic_sq+=(h-b)*(h-b);
        den+=b*b;
        max_apply=std::max(max_apply,std::abs(s-e)/
                           std::max(1.L,std::abs(e)));
    }
    if(!(den>0.))throw std::invalid_argument("V55 gate: zero RHS norm");
    for(std::size_t i=0;i<split.face_watts.size();++i){
        const long double e=long_double.face_watts[i];
        const long double s=dom_to_long_double_v55(split.face_watts[i]);
        faces=std::max(faces,std::abs(s-e));
    }
    DomWlsHiloGateV55 report;
    report.hilo_l2=static_cast<double>(std::sqrt(split_sq/den));
    report.extended_l2=static_cast<double>(std::sqrt(extended_sq/den));
    report.historic_double_l2=static_cast<double>(std::sqrt(historic_sq/den));
    report.absolute_face_error_watts=static_cast<double>(faces);
    report.max_apply_difference=static_cast<double>(max_apply);
    report.hilo_pass=report.hilo_l2<=strict_tolerance;
    report.original_double_pass=report.historic_double_l2<=strict_tolerance;
    // Deliberately not equivalent to a physical convergence certificate.
    report.publish=false;
    return report;
}
} // namespace cfdx::physics
