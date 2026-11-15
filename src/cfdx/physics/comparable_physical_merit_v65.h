#pragma once
// V65: opt-in, state-only, same-operator merit for the existing V47 backtracking.
// This does NOT define a new physical model, alter the primal F(X), or loosen
// final convergence gates. Both F(X_old) and F(X_trial) must be independently
// reassembled by the same V51 callback at every trial. Scaling is FROZEN from
// X_old for that trial; per-state normalization is not compared directly.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct ComparablePhysicalMeritV65 {
    std::array<double,6> raw_old{},raw_trial{},reference_scale{};
    std::array<double,6> scaled_old{},scaled_trial{},growth{};
    double old_merit=0., trial_merit=0., merit_ratio=0.;
    double largest_growth=0.;
    std::size_t worst_block=0;
    bool baseline=false;
    bool joint_allowed=false;
    bool monotone_allowed=false;
    bool accepted=false;
};

inline ComparablePhysicalMeritV65 evaluate_comparable_merit_v65(
    const PhysicalCandidateV51& old_state,
    const PhysicalCandidateV51& candidate,
    const PhysicalResidualV51& old_residual,
    const PhysicalResidualV51& trial_residual,
    double relative_tolerance, double continuity_tolerance,
    bool baseline=false, bool certified_physical_convergence=false)
{
    validate_physical_candidate_v51(old_state);
    validate_physical_candidate_v51(candidate);
    if (old_state.cells!=candidate.cells ||
        !(relative_tolerance>0.) || !std::isfinite(relative_tolerance) ||
        !(continuity_tolerance>0.) || !std::isfinite(continuity_tolerance))
        throw std::invalid_argument("V65 incompatible physical states or tolerances");
    validate_physical_result_v51(old_state,old_residual);
    validate_physical_result_v51(candidate,trial_residual);
    const std::size_t n=old_state.cells;
    ComparablePhysicalMeritV65 result;
    result.baseline=baseline;
    result.joint_allowed=true;
    constexpr double joint_growth=1.20; // unchanged V13 acceptance envelope
    constexpr double max_merit_ratio=0.995; // unchanged V14 requirement
    for (std::size_t b=0;b<6;++b) {
        double old_normalized=0.;
        for (std::size_t i=0;i<n;++i) {
            const auto j=b*n+i;
            result.raw_old[b]=std::max(result.raw_old[b],std::abs(old_residual.signed_raw[j]));
            result.raw_trial[b]=std::max(result.raw_trial[b],std::abs(trial_residual.signed_raw[j]));
            old_normalized=std::max(old_normalized,std::abs(old_residual.normalized[j]));
        }
        // Match the BASELINE block's physical normalization, then freeze it.
        // Do NOT divide the candidate by its own changed RHS or SST row scale.
        result.reference_scale[b]=(result.raw_old[b]>0. && old_normalized>0.)
            ? result.raw_old[b]/old_normalized : 1.;
        if (!(result.reference_scale[b]>0.) ||
            !std::isfinite(result.reference_scale[b]))
            throw std::runtime_error("V65 invalid frozen physical scale");
        result.scaled_old[b]=result.raw_old[b]/result.reference_scale[b];
        result.scaled_trial[b]=result.raw_trial[b]/result.reference_scale[b];
        if (!std::isfinite(result.scaled_old[b]) || !std::isfinite(result.scaled_trial[b]))
            throw std::runtime_error("V65 nonfinite frozen-scale residual");
        const double tolerance=b==3?continuity_tolerance:relative_tolerance;
        const double denominator=std::max(result.scaled_old[b],10.*tolerance);
        result.growth[b]=result.scaled_trial[b]/denominator;
        if (result.growth[b]>result.largest_growth) {
            result.largest_growth=result.growth[b];
            result.worst_block=b;
        }
        if (result.growth[b]>joint_growth)
            result.joint_allowed=false;
        result.old_merit=std::max(result.old_merit,result.scaled_old[b]);
        result.trial_merit=std::max(result.trial_merit,result.scaled_trial[b]);
    }
    result.merit_ratio=result.trial_merit/std::max(result.old_merit,1e-300);
    if (!std::isfinite(result.merit_ratio))
        throw std::runtime_error("V65 nonfinite comparable merit ratio");
    result.monotone_allowed=certified_physical_convergence ||
        result.merit_ratio<=max_merit_ratio;
    // Preserve V13/V14 first-commit behavior: the first trial establishes a
    // reference, but is never called a physical convergence certificate.
    result.accepted=baseline || (result.joint_allowed && result.monotone_allowed);
    return result;
}

struct ComparableFluxAuditV65 {
    double reconstructed_linf=0.;
    double reconstructed_relative_l1=0.;
    double continuity_linf=0.;
    double authoritative_relative_l1=0.; // only meaningful for previous committed field
};

template<class Owner, class Neighbour>
ComparableFluxAuditV65 audit_comparable_flux_v65(
    const PhysicalCandidateV51& state,
    const PhysicalResidualV51& residual,
    const std::vector<double>& independently_reconstructed,
    Owner&& owner, Neighbour&& neighbour,
    const std::vector<double>* previous_authoritative=nullptr)
{
    validate_physical_result_v51(state,residual);
    if (independently_reconstructed.size()!=residual.face_flux.size() ||
        (previous_authoritative &&
         previous_authoritative->size()!=residual.face_flux.size()))
        throw std::invalid_argument("V65 inconsistent face flux vector sizes");
    std::vector<long double> divergence(state.cells,0.0L);
    ComparableFluxAuditV65 audit;
    long double reference_l1=0.,diff_l1=0.,historical_diff=0.;
    for (std::size_t f=0;f<residual.face_flux.size();++f) {
        const auto phi=residual.face_flux[f],again=independently_reconstructed[f];
        if (!std::isfinite(phi) || !std::isfinite(again))
            throw std::runtime_error("V65 nonfinite reconstructed flux");
        const auto o=static_cast<std::size_t>(owner(f));
        const auto nb=neighbour(f);
        if (o>=state.cells || nb < -1 ||
            (nb>=0 && static_cast<std::size_t>(nb)>=state.cells))
            throw std::runtime_error("V65 invalid face topology");
        divergence[o]+=static_cast<long double>(phi);
        if (nb>=0)divergence[static_cast<std::size_t>(nb)]-=static_cast<long double>(phi);
        reference_l1+=std::abs(static_cast<long double>(phi));
        diff_l1+=std::abs(static_cast<long double>(phi)-static_cast<long double>(again));
        audit.reconstructed_linf=std::max(audit.reconstructed_linf,std::abs(phi-again));
        if (previous_authoritative) {
            if (!std::isfinite((*previous_authoritative)[f]))
                throw std::runtime_error("V65 nonfinite previous authoritative flux");
            historical_diff+=std::abs(static_cast<long double>(phi)-
                                      static_cast<long double>((*previous_authoritative)[f]));
        }
    }
    const auto divisor=std::max(reference_l1,1.0e-300L);
    audit.reconstructed_relative_l1=static_cast<double>(diff_l1/divisor);
    audit.authoritative_relative_l1=static_cast<double>(historical_diff/divisor);
    for (std::size_t i=0;i<state.cells;++i)
        audit.continuity_linf=std::max(audit.continuity_linf,
            std::abs(static_cast<double>(divergence[i])-
                     residual.signed_raw[3*state.cells+i]));
    if (!std::isfinite(audit.reconstructed_relative_l1) ||
        !std::isfinite(audit.continuity_linf) ||
        audit.reconstructed_relative_l1>1e-10 ||
        audit.continuity_linf>1e-10)
        throw std::runtime_error("V65 Rhie-Chow/continuity operator identity not closed");
    return audit;
}

} // namespace cfdx::physics
