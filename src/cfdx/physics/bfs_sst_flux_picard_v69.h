@@ -0,0 +1,144 @@
#pragma once
// V69. Experimental, strictly read-only SST/RC fixed-point audit.
// The caller supplies the *production* SST transport as a pure transformation
// of a private six-field candidate with a fixed RC face flux. The physical
// evaluator is exactly V51 F_phys, not a new nonlinear merit.
// There is no outer transaction commit or user-visible flow update here.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include "cfdx/physics/comparable_physical_merit_v65.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct BfsSstFluxPicardStepV69 {
    std::size_t cycle=0;
    double omega_before=0.,omega_after=0.,k_before=0.,k_after=0.;
    double omega_ratio=0.,k_ratio=0.;
    double merit_before=0.,merit_after=0.,merit_ratio=0.;
    double continuity_before=0.,continuity_after=0.; // signed FV divergence L1
    double rc_flux_delta_linf=0.; // new F_phys RC phi - phi used for SST transport
    double nut_delta_linf=0.;
    bool rc_closed=false;
    // Counterfactual V65 decision on this same old state with frozen scales.
    // Diagnostic only: no production transaction is committed by this API.
    double v65_merit_ratio=0.,v65_worst_growth=0.;
    std::size_t v65_worst_block=0;
    bool v65_joint_ok=false,v65_monotone_ok=false,v65_would_pass=false;
};

inline double bfs_sst_block_linf_v69(const PhysicalResidualV51& f,
                                   std::size_t cells,std::size_t block) {
    if (f.signed_raw.size()!=6*cells || block>=6)
        throw std::invalid_argument("V69 malformed physical block");
    double worst=0.;
    for(std::size_t i=0;i<cells;++i)
        worst=std::max(worst,std::abs(f.signed_raw[block*cells+i]));
    return worst;
}

inline double bfs_sst_continuity_l1_v69(const PhysicalResidualV51& f,
                                       std::size_t cells) {
    if(f.signed_raw.size()!=6*cells)
        throw std::invalid_argument("V69 malformed continuity block");
    double sum=0.;
    for(std::size_t i=0;i<cells;++i)
        sum+=std::abs(f.signed_raw[3*cells+i]);
    return sum;
}

// Correct must never mutate its input or modify U,p; it returns a new state
// after a genuine converged SST solve using *only* the supplied RC face flux.
// At each cycle the full RC/rAU/mu_t is reassembled by the V51 callback.
// This API deliberately does not evaluate the V65 acceptance rule: it cannot
// accept any candidate or replace the authoritative COUPLED conservative phi.
template<class Evaluate,class Correct>
std::vector<BfsSstFluxPicardStepV69> probe_bfs_sst_flux_picard_v69(
    const PhysicalCandidateV51& start,const PhysicalResidualV51& initial,
    Evaluate&& evaluate,Correct&& correct,
    double k_min,double omega_min,std::size_t cycles=2,
    const PhysicalCandidateV51* old_state=nullptr,
    const PhysicalResidualV51* old_result=nullptr,
    double relative_tolerance=1e-8,double continuity_tolerance=1e-8,
    PhysicalCandidateV51* private_final_state=nullptr,
    PhysicalResidualV51* private_final_residual=nullptr) {
    if(cycles==0 || cycles>4)
        throw std::invalid_argument("V69 requires 1..4 fixed-point cycles");
    if ((old_state==nullptr)!=(old_result==nullptr))
        throw std::invalid_argument("V69 V65 baseline state/result must be supplied together");
    validate_physical_candidate_v51(start,k_min,omega_min);
    validate_physical_result_v51(start,initial);
    auto state=start;
    auto previous=initial;
    std::vector<BfsSstFluxPicardStepV69> evidence;
    evidence.reserve(cycles);
    const auto n=start.cells;
    for(std::size_t step=0;step<cycles;++step) {
        const auto frozen_state=state;
        const auto trial=correct(frozen_state,previous.face_flux);
        if(frozen_state.fields!=state.fields)
            throw std::runtime_error("V69 SST solver mutated an input candidate");
        validate_physical_candidate_v51(trial,k_min,omega_min);
        if(trial.cells!=n || trial.fields.size()!=state.fields.size())
            throw std::runtime_error("V69 SST solver returned incompatible state");
        // Exactly Ux,Uy,Uz,p must be identical: SST is the only free block.
        if(!std::equal(state.fields.begin(),state.fields.begin()+static_cast<std::ptrdiff_t>(4*n),
                       trial.fields.begin()))
            throw std::runtime_error("V69 SST-only corrector changed U or p");
        const auto physical=evaluate_physical_candidate_v51(
            trial,evaluate,k_min,omega_min);
        const auto replay=evaluate_physical_candidate_v51(
            trial,evaluate,k_min,omega_min);
        verify_candidate_replay_v51(physical,replay);
        if(previous.face_flux.size()!=physical.face_flux.size() ||
           previous.eddy_viscosity.size()!=physical.eddy_viscosity.size())
            throw std::runtime_error("V69 RC flux or nu_t shape changed");
        BfsSstFluxPicardStepV69 row;
        row.cycle=step+1;
        row.k_before=bfs_sst_block_linf_v69(previous,n,4);
        row.k_after=bfs_sst_block_linf_v69(physical,n,4);
        row.omega_before=bfs_sst_block_linf_v69(previous,n,5);
        row.omega_after=bfs_sst_block_linf_v69(physical,n,5);
        row.k_ratio=row.k_after/std::max(row.k_before,1e-300);
        row.omega_ratio=row.omega_after/std::max(row.omega_before,1e-300);
        row.continuity_before=bfs_sst_continuity_l1_v69(previous,n);
        row.continuity_after=bfs_sst_continuity_l1_v69(physical,n);
        row.merit_before=previous.merit;
        row.merit_after=physical.merit;
        row.merit_ratio=row.merit_after/std::max(row.merit_before,1e-300);
        for(std::size_t f=0;f<previous.face_flux.size();++f)
            row.rc_flux_delta_linf=std::max(row.rc_flux_delta_linf,
                std::abs(physical.face_flux[f]-previous.face_flux[f]));
        for(std::size_t i=0;i<n;++i)
            row.nut_delta_linf=std::max(row.nut_delta_linf,
                std::abs(physical.eddy_viscosity[i]-previous.eddy_viscosity[i]));
        row.rc_closed=physical.flux_closed;
        if(old_state) {
            const auto gate=evaluate_comparable_merit_v65(
                *old_state,trial,*old_result,physical,
                relative_tolerance,continuity_tolerance,false,false);
            row.v65_merit_ratio=gate.merit_ratio;
            row.v65_worst_growth=gate.largest_growth;
            row.v65_worst_block=gate.worst_block;
            row.v65_joint_ok=gate.joint_allowed;
            row.v65_monotone_ok=gate.monotone_allowed;
            row.v65_would_pass=gate.accepted;
        }
        if(!std::isfinite(row.k_ratio) || !std::isfinite(row.omega_ratio) ||
           !std::isfinite(row.merit_ratio) || !std::isfinite(row.continuity_after) ||
           !std::isfinite(row.rc_flux_delta_linf) || !std::isfinite(row.nut_delta_linf))
            throw std::runtime_error("V69 non-finite SST flux cycle evidence");
        evidence.push_back(row);
        state=trial;
        previous=physical;
    }
    // V70 can request the final PRIVATE state for a subsequent conservative
    // pressure/flux trial. No production carrier or SST state is published.
    if(private_final_state)*private_final_state=state;
    if(private_final_residual)*private_final_residual=previous;
    return evidence;
}
}
