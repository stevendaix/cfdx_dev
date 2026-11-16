#pragma once
// V73: after-SST, same-state diagnostic of the EXISTING V71 Picard direction.
// Does not introduce a new Schur/Newton, alter the primal FV equations, or
// relax V65. All final states are private; a native momentum/SST re-solve
// (which this experiment does NOT perform) is mandatory before any commit.
#include "cfdx/physics/bfs_joint_fixed_point_v72.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct BfsFaceCellDeltaV73 {
    std::size_t cell=0;
    double raux_before=0.,raux_after=0.;
    double nut_before=0.,nut_after=0.;
    double momentum_x_before=0.,momentum_x_after=0.;
    double omega_signed_before=0.,omega_signed_after=0.;
    double k_signed_before=0.,k_signed_after=0.;
};

struct BfsDirectionAttemptV73 {
    double gamma=0.,pressure_fraction=0.;
    double face_flux_before=0.,face_flux_after=0.;
    double mass_before=0.,mass_after=0.;
    double merit_ratio=std::numeric_limits<double>::infinity();
    double largest_growth=std::numeric_limits<double>::infinity();
    std::array<double,6> frozen_scaled_before{},frozen_scaled_after{};
    std::array<BfsFaceCellDeltaV73,2> cells{};
    bool evaluated=false,guards_ok=false,mass_improved=false;
    bool momentum_improved=false,k_improved=false,omega_improved=false;
    bool all_relevant_blocks_improved=false;
    bool final_mass_closed=false,all_six_blocks_closed=false;
    bool eligible_for_commit=false;
    std::string rejection;
};

struct BfsDirectionProbeV73 {
    std::vector<BfsDirectionAttemptV73> attempts;
    std::size_t evaluated_count=0,guard_accepted_count=0;
    std::size_t simultaneous_descent_count=0,commits=0;
    double initial_mass=0.,smallest_mass=std::numeric_limits<double>::infinity();
    bool native_momentum_linear_recertified=false;
    bool native_sst_transport_recertified=false;
    bool commit_authorized=false;
};

// Compare the same physical state and flux after the SST update. The only
// new degrees of freedom are the two damping factors of the existing V71
// diagonal-Picard predictor and V70 pressure response. There is no interpolation
// of k/omega: transport is *solved again* on the updated physical RC flux.
// The previous V72 privately accepted state is the source. One trial cannot
// affect another trial. A candidate with a smaller FV mass residual is not
// declared conservative or eligible for production commit.
template<class Evaluate,class Project,class Transport,class Owner,class Neighbour,
         class Observe>
BfsDirectionProbeV73 bfs_probe_post_sst_directions_v73(
    const PhysicalCandidateV51& old_state,const PhysicalResidualV51& old_f,
    const PhysicalCandidateV51& post_sst,const PhysicalResidualV51& post_sst_f,
    Evaluate&& evaluate,Project&& project,Transport&& transport,
    Owner&& owner,Neighbour&& neighbour,
    std::size_t interface_face,std::size_t cell_a,std::size_t cell_b,
    const std::vector<BfsJointChoiceV72>& choices,Observe&& observe,
    double k_min,double omega_min,double relative_tolerance,
    double continuity_tolerance,double final_mass_tolerance=1e-10) {
    validate_physical_result_v51(old_state,old_f);
    validate_physical_result_v51(post_sst,post_sst_f);
    const auto n=post_sst.cells;
    if(n!=old_state.cells || n>152 || !n || cell_a>=n || cell_b>=n ||
       interface_face>=post_sst_f.face_flux.size() || choices.empty() ||
       choices.size()>24 || !(final_mass_tolerance>0.) ||
       !std::isfinite(final_mass_tolerance) ||
       static_cast<std::size_t>(owner(interface_face))!=cell_a ||
       neighbour(interface_face)!=static_cast<int>(cell_b))
        throw std::invalid_argument("V73 incompatible face or bounded post-SST trial domain");
    const auto repeated=evaluate_physical_candidate_v51(
        post_sst,evaluate,k_min,omega_min);
    verify_candidate_replay_v51(post_sst_f,repeated);
    BfsDirectionProbeV73 report;
    report.initial_mass=bfs_mass_balance_v70(n,post_sst_f.face_flux,
                                              owner,neighbour).divergence_linf;
    const auto source_fields=post_sst.fields;
    const auto source_flux=post_sst_f.face_flux;
    const auto source_rau=post_sst_f.inverse_momentum_diagonal;
    const auto source_nut=post_sst_f.eddy_viscosity;
    const auto old_merit=evaluate_comparable_merit_v65(old_state,post_sst,
        old_f,post_sst_f,relative_tolerance,continuity_tolerance,false,false);
    for(const auto& choice:choices) {
        BfsDirectionAttemptV73 item;
        item.gamma=choice.momentum_fraction;
        item.pressure_fraction=choice.pressure_fraction;
        item.face_flux_before=post_sst_f.face_flux[interface_face];
        item.mass_before=report.initial_mass;
        try {
            PhysicalCandidateV51 candidate;
            PhysicalResidualV51 physical;
            (void)bfs_joint_private_cycle_v71(old_state,old_f,post_sst,
                post_sst_f,evaluate,project,transport,owner,neighbour,
                choice.momentum_fraction,choice.pressure_fraction,
                k_min,omega_min,relative_tolerance,continuity_tolerance,
                final_mass_tolerance,2,&candidate,&physical);
            // Strict isolation: all trials begin at the very same post-SST state.
            if(post_sst.fields!=source_fields ||
               post_sst_f.face_flux!=source_flux ||
               post_sst_f.inverse_momentum_diagonal!=source_rau ||
               post_sst_f.eddy_viscosity!=source_nut)
                throw std::runtime_error("V73 candidate contaminated the source state");
            const auto cert=bfs_recertify_six_field_v72(old_state,old_f,
                candidate,physical,evaluate,owner,neighbour,
                relative_tolerance,continuity_tolerance,
                final_mass_tolerance,k_min,omega_min);
            const auto merit=evaluate_comparable_merit_v65(old_state,
                candidate,old_f,physical,relative_tolerance,
                continuity_tolerance,false,false);
            item.evaluated=true;
            item.mass_after=cert.mass_divergence_linf;
            item.face_flux_after=physical.face_flux[interface_face];
            item.merit_ratio=merit.merit_ratio;
            item.largest_growth=merit.largest_growth;
            item.guards_ok=merit.accepted && cert.original_v65_guards_ok;
            item.mass_improved=item.mass_after<report.initial_mass;
            item.final_mass_closed=cert.conservation_reassembled_closed;
            item.all_six_blocks_closed=cert.all_six_physical_blocks_closed;
            item.eligible_for_commit=cert.eligible_for_production_commit;
            item.frozen_scaled_before=old_merit.scaled_trial;
            item.frozen_scaled_after=merit.scaled_trial;
            for(std::size_t a=0;a<2;++a) {
                auto& c=item.cells[a];c.cell=a==0?cell_a:cell_b;
                const auto i=c.cell;
                c.raux_before=post_sst_f.inverse_momentum_diagonal[i];
                c.raux_after=physical.inverse_momentum_diagonal[i];
                c.nut_before=post_sst_f.eddy_viscosity[i];
                c.nut_after=physical.eddy_viscosity[i];
                c.momentum_x_before=post_sst_f.signed_raw[i];
                c.momentum_x_after=physical.signed_raw[i];
                c.k_signed_before=post_sst_f.signed_raw[4*n+i];
                c.k_signed_after=physical.signed_raw[4*n+i];
                c.omega_signed_before=post_sst_f.signed_raw[5*n+i];
                c.omega_signed_after=physical.signed_raw[5*n+i];
            }
            const double momentum_before=std::max({old_merit.scaled_trial[0],
                old_merit.scaled_trial[1],old_merit.scaled_trial[2]});
            const double momentum_after=std::max({merit.scaled_trial[0],
                merit.scaled_trial[1],merit.scaled_trial[2]});
            item.momentum_improved=momentum_after<momentum_before;
            item.k_improved=merit.scaled_trial[4]<old_merit.scaled_trial[4];
            item.omega_improved=merit.scaled_trial[5]<old_merit.scaled_trial[5];
            item.all_relevant_blocks_improved=item.momentum_improved &&
                item.mass_improved && item.k_improved && item.omega_improved;
            if(item.eligible_for_commit)
                throw std::logic_error("V73 private trial unexpectedly authorized a commit");
            // The native source-term audit must pass BEFORE counting an
            // evaluation as usable; otherwise there would be false evidence.
            observe(candidate,physical,item);
            ++report.evaluated_count;
            report.smallest_mass=std::min(report.smallest_mass,item.mass_after);
            if(item.guards_ok)++report.guard_accepted_count;
            if(item.all_relevant_blocks_improved)
                ++report.simultaneous_descent_count;
        } catch(const std::exception& error) {
            item.rejection=error.what();
            item.evaluated=false;
            // A candidate rejected at any stage is never eligible for commit.
            item.eligible_for_commit=false;
        }
        report.attempts.push_back(std::move(item));
    }
    report.commit_authorized=false;
    return report;
}
} // namespace cfdx::physics