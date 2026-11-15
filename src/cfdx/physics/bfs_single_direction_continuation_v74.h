#pragma once
// V74: repeat ONE V73-certified direction on private, reassembled physical
// states. No grid search, no change to primal equations, no production commit.
// The native linear momentum / SST certificates remain unavailable here.
#include "cfdx/physics/bfs_joint_fixed_point_v72.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct BfsContinuationRowV74 {
    std::size_t cycle=0;
    double mass_before=0.,mass_after=0.;
    double merit_before=0.,merit_after=0.;
    double merit_ratio_to_original=0.,largest_growth=0.;
    double phi3_before=0.,phi3_after=0.;
    double raux0_before=0.,raux0_after=0.;
    double raux12_before=0.,raux12_after=0.;
    double nut0_before=0.,nut0_after=0.;
    double nut12_before=0.,nut12_after=0.;
    std::array<double,6> blocks_before{},blocks_after{};
    bool evaluated=false,guards_ok=false,strict_mass_descent=false;
    bool merit_descends=false,relevant_groups_nonincreasing=false;
    bool selected=false,mass_closed=false,six_blocks_closed=false;
    std::string reason;
};
using BfsContinuationObserverV74=std::function<void(
    std::size_t,const PhysicalCandidateV51&,const PhysicalResidualV51&,
    const PhysicalCandidateV51&,const PhysicalResidualV51&,
    const BfsContinuationRowV74&)>;
struct BfsContinuationV74 {
    std::vector<BfsContinuationRowV74> rows;
    std::size_t accepted_private_cycles=0;
    double initial_mass=0.,final_mass=0.;
    bool stopped=true,final_mass_closed=false;
    bool native_momentum_certificate=false,native_sst_certificate=false;
    bool production_commit_authorized=false;
    std::string status="not_run";
};

template<class Evaluate,class Project,class Transport,class Owner,class Neighbour>
BfsContinuationV74 bfs_continue_single_direction_v74(
    const PhysicalCandidateV51& original,const PhysicalResidualV51& original_f,
    const PhysicalCandidateV51& start,const PhysicalResidualV51& start_f,
    Evaluate&& evaluate,Project&& project,Transport&& transport,
    Owner&& owner,Neighbour&& neighbour,
    std::size_t interface_face,std::size_t cell_a,std::size_t cell_b,
    std::size_t max_cycles,double k_min,double omega_min,
    double relative_tolerance,double continuity_tolerance,
    double mass_tolerance=1e-10,
    PhysicalCandidateV51* final_private=nullptr,
    PhysicalResidualV51* final_result=nullptr,
    const BfsContinuationObserverV74& on_evaluated={}) {
    validate_physical_result_v51(original,original_f);
    validate_physical_result_v51(start,start_f);
    if(!start.cells || start.cells>152 || original.cells!=start.cells ||
       cell_a>=start.cells || cell_b>=start.cells ||
       interface_face>=start_f.face_flux.size() ||
       static_cast<std::size_t>(owner(interface_face))!=cell_a ||
       neighbour(interface_face)!=static_cast<int>(cell_b) ||
       max_cycles==0 || max_cycles>8 || !(mass_tolerance>0.) ||
       !std::isfinite(mass_tolerance))
        throw std::invalid_argument("V74 invalid bounded fixed-direction domain");
    verify_candidate_replay_v51(start_f,
        evaluate_physical_candidate_v51(start,evaluate,k_min,omega_min));
    auto current=start;
    auto current_f=start_f;
    BfsContinuationV74 report;
    report.initial_mass=bfs_mass_balance_v70(start.cells,start_f.face_flux,
        owner,neighbour).divergence_linf;
    report.final_mass=report.initial_mass;
    constexpr double gamma=1./128.,alpha_p=1./16.;
    for(std::size_t cycle=1;cycle<=max_cycles;++cycle) {
        BfsContinuationRowV74 row;row.cycle=cycle;
        row.mass_before=report.final_mass;
        row.phi3_before=current_f.face_flux[interface_face];
        row.raux0_before=current_f.inverse_momentum_diagonal[cell_a];
        row.raux12_before=current_f.inverse_momentum_diagonal[cell_b];
        row.nut0_before=current_f.eddy_viscosity[cell_a];
        row.nut12_before=current_f.eddy_viscosity[cell_b];
        const auto before_merit=evaluate_comparable_merit_v65(
            original,current,original_f,current_f,
            relative_tolerance,continuity_tolerance,false,false);
        row.merit_before=before_merit.trial_merit;
        row.blocks_before=before_merit.scaled_trial;
        const auto snapshot=current.fields;
        const auto previous_phi=current_f.face_flux;
        const auto previous_rau=current_f.inverse_momentum_diagonal;
        const auto previous_nut=current_f.eddy_viscosity;
        try {
            PhysicalCandidateV51 trial;
            PhysicalResidualV51 trial_f;
            (void)bfs_joint_private_cycle_v71(original,original_f,current,current_f,
                evaluate,project,transport,owner,neighbour,
                gamma,alpha_p,k_min,omega_min,relative_tolerance,
                continuity_tolerance,mass_tolerance,cycle,&trial,&trial_f);
            if(current.fields!=snapshot || current_f.face_flux!=previous_phi ||
               current_f.inverse_momentum_diagonal!=previous_rau ||
               current_f.eddy_viscosity!=previous_nut)
                throw std::runtime_error("V74 operator mutated its private input");
            const auto cert=bfs_recertify_six_field_v72(
                original,original_f,trial,trial_f,evaluate,owner,neighbour,
                relative_tolerance,continuity_tolerance,mass_tolerance,k_min,omega_min);
            const auto merit=evaluate_comparable_merit_v65(
                original,trial,original_f,trial_f,
                relative_tolerance,continuity_tolerance,false,false);
            row.evaluated=true;
            row.mass_after=cert.mass_divergence_linf;
            row.phi3_after=trial_f.face_flux[interface_face];
            row.raux0_after=trial_f.inverse_momentum_diagonal[cell_a];
            row.raux12_after=trial_f.inverse_momentum_diagonal[cell_b];
            row.nut0_after=trial_f.eddy_viscosity[cell_a];
            row.nut12_after=trial_f.eddy_viscosity[cell_b];
            row.merit_after=merit.trial_merit;
            row.blocks_after=merit.scaled_trial;
            row.merit_ratio_to_original=merit.merit_ratio;
            row.largest_growth=merit.largest_growth;
            row.guards_ok=merit.accepted && cert.original_v65_guards_ok;
            row.strict_mass_descent=row.mass_after<=mass_tolerance ||
                row.mass_after < row.mass_before*(1.-1e-12);
            row.merit_descends=row.merit_after<row.merit_before*(1.-1e-12);
            // In an extruded 2D case the Uz residual can be at roundoff
            // while Ux/Uy contract. Do NOT require strict decrease of an
            // inactive component below its numerical floor. Preserve the
            // original V65 per-block growth guard for all six components,
            // but require the four PHYSICAL GROUPS to stay nonincreasing.
            const double momentum_before=std::max({
                row.blocks_before[0],row.blocks_before[1],row.blocks_before[2]});
            const double momentum_after=std::max({
                row.blocks_after[0],row.blocks_after[1],row.blocks_after[2]});
            row.relevant_groups_nonincreasing=
                momentum_after<=momentum_before*(1.+1e-10)+1e-13 &&
                row.blocks_after[3]<=row.blocks_before[3]*(1.+1e-10)+1e-13 &&
                row.blocks_after[4]<=row.blocks_before[4]*(1.+1e-10)+1e-13 &&
                row.blocks_after[5]<=row.blocks_before[5]*(1.+1e-10)+1e-13;
            row.mass_closed=cert.conservation_reassembled_closed;
            row.six_blocks_closed=cert.all_six_physical_blocks_closed;
            row.selected=row.guards_ok && row.strict_mass_descent &&
                row.merit_descends && row.relevant_groups_nonincreasing;
            // V75 observer may measure both SST budgets even for a rejected
            // cycle. It cannot modify the const states or authorize a commit.
            if(on_evaluated)on_evaluated(cycle,current,current_f,trial,trial_f,row);
            if(!row.selected) {
                if(!row.guards_ok) row.reason="V65_guards";
                else if(!row.strict_mass_descent)row.reason="reassembled_mass_not_decreasing";
                else if(!row.merit_descends)row.reason="same_scale_merit_not_decreasing";
                else row.reason="one_or_more_physical_groups_increasing";
            } else {
                current=std::move(trial);
                current_f=std::move(trial_f);
                report.final_mass=row.mass_after;
                report.final_mass_closed=row.mass_closed;
                ++report.accepted_private_cycles;
            }
        } catch(const std::exception& ex){
            row.reason=ex.what();row.selected=false;
        }
        const bool chosen=row.selected;
        report.rows.push_back(std::move(row));
        if(!chosen){report.status="failed_closed_no_monotone_common_state_direction";break;}
        if(report.final_mass_closed){
            report.status="physical_mass_closed_no_native_linear_certificates";break;
        }
    }
    if(report.status=="not_run")report.status="bounded_cycle_budget_exhausted";
    if(final_private)*final_private=current;
    if(final_result)*final_result=current_f;
    // No V74 private diagonal predictor can issue a native momentum linear
    // certificate, nor certify the SST transport on a finally committed flux.
    report.production_commit_authorized=false;
    return report;
}
} // namespace cfdx::physics
