@@ -0,0 +1,130 @@
#pragma once
// CFDX V67: READ-ONLY counterfactual block damping study on one native
// six-field state. The actual candidate, residual, SGS/SST transport, RC face
// flux and nonlinear acceptance gates remain untouched. The synthetic mixed
// candidates are NOT SST transport solutions and may not be committed.
#include "cfdx/physics/comparable_physical_merit_v65.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct BfsBlockMixV67 {
    double k_fraction;
    double omega_fraction;
};

// Fraction 1 keeps the converged SST trial; 0 freezes that SST component at
// the old state. Momentum and gauge-fixed pressure always come from the
// same carrier trial, never from another retry.
inline PhysicalCandidateV51 mixed_sst_candidate_v67(
    const PhysicalCandidateV51& old,
    const PhysicalCandidateV51& trial,
    const BfsBlockMixV67 fractions,
    const double k_min, const double omega_min)
{
    validate_physical_candidate_v51(old,k_min,omega_min);
    validate_physical_candidate_v51(trial,k_min,omega_min);
    if (old.cells!=trial.cells ||
        !std::isfinite(fractions.k_fraction) ||
        !std::isfinite(fractions.omega_fraction) ||
        fractions.k_fraction<0.0 || fractions.k_fraction>1.0 ||
        fractions.omega_fraction<0.0 || fractions.omega_fraction>1.0)
        throw std::invalid_argument("V67 invalid block damping fractions or state");
    PhysicalCandidateV51 mixed=trial;
    for (std::size_t i=0;i<trial.cells;++i) {
        const std::size_t k=4*trial.cells+i;
        const std::size_t omega=5*trial.cells+i;
        mixed.fields[k]=old.fields[k]+fractions.k_fraction*
            (trial.fields[k]-old.fields[k]);
        mixed.fields[omega]=old.fields[omega]+fractions.omega_fraction*
            (trial.fields[omega]-old.fields[omega]);
    }
    validate_physical_candidate_v51(mixed,k_min,omega_min);
    return mixed;
}

// Unlike carrier-only V47 damping, this is a true straight segment in
// the six physical unknowns. It is evaluated by F_phys but never committed.
inline PhysicalCandidateV51 mixed_six_candidate_v67(
    const PhysicalCandidateV51& old,
    const PhysicalCandidateV51& trial,
    const double alpha,const double k_min,const double omega_min)
{
    validate_physical_candidate_v51(old,k_min,omega_min);
    validate_physical_candidate_v51(trial,k_min,omega_min);
    if (old.cells!=trial.cells || !(alpha>0.0) || alpha>1.0 ||
        !std::isfinite(alpha))
        throw std::invalid_argument("V67 full-six invalid state or step");
    PhysicalCandidateV51 mixed=old;
    for (std::size_t j=0;j<old.fields.size();++j)
        mixed.fields[j]+=alpha*(trial.fields[j]-old.fields[j]);
    validate_physical_candidate_v51(mixed,k_min,omega_min);
    return mixed;
}

struct BfsBlockProbeResultV67 {
    BfsBlockMixV67 fractions{};
    double merit_ratio=0.0;
    double worst_growth=0.0;
    std::size_t worst_block=0;
    bool joint_allowed=false;
    bool monotone_allowed=false;
    bool would_pass_original_guards=false;
};

template<class Evaluate>
inline BfsBlockProbeResultV67 probe_sst_block_mix_v67(
    const PhysicalCandidateV51& old,
    const PhysicalCandidateV51& trial,
    const PhysicalResidualV51& old_residual,
    const BfsBlockMixV67 fractions,
    Evaluate&& evaluate,
    const double relative_tolerance,
    const double continuity_tolerance,
    const double k_min, const double omega_min)
{
    const auto mixed=mixed_sst_candidate_v67(
        old,trial,fractions,k_min,omega_min);
    const auto physical=evaluate_physical_candidate_v51(
        mixed,evaluate,k_min,omega_min);
    const auto merit=evaluate_comparable_merit_v65(
        old,mixed,old_residual,physical,
        relative_tolerance,continuity_tolerance,false,false);
    return {fractions,merit.merit_ratio,merit.largest_growth,
            merit.worst_block,merit.joint_allowed,merit.monotone_allowed,
            merit.accepted};
}

template<class Evaluate>
inline ComparablePhysicalMeritV65 probe_full_six_line_v67(
    const PhysicalCandidateV51& old,
    const PhysicalCandidateV51& trial,
    const PhysicalResidualV51& old_residual,
    const double alpha,Evaluate&& evaluate,
    const double relative_tolerance,const double continuity_tolerance,
    const double k_min,const double omega_min)
{
    const auto candidate=mixed_six_candidate_v67(
        old,trial,alpha,k_min,omega_min);
    const auto residual=evaluate_physical_candidate_v51(
        candidate,evaluate,k_min,omega_min);
    return evaluate_comparable_merit_v65(
        old,candidate,old_residual,residual,
        relative_tolerance,continuity_tolerance,false,false);
}

// This is a mathematical diagnostic, NOT a V65 acceptance condition. The
// unchanged V65 0.995 is reported alongside an Armijo-style comparison
// with the same full-step fractional decrease (0.005). The latter is scaled
// by the actual carrier step beta and cannot be used to commit a state.
inline double v67_armijo_reference_ratio(const double carrier_beta) {
    if (!(carrier_beta>0.0) || carrier_beta>1.0 ||
        !std::isfinite(carrier_beta))
        throw std::invalid_argument("V67 invalid carrier step");
    return 1.0-0.005*carrier_beta;
}

} // namespace cfdx::physics
