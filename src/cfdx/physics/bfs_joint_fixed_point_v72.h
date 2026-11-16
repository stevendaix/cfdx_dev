#pragma once
// V72: bounded, private-only continuation of V71's EXISTING Picard operators.
// No alternative Schur, no alternate F_phys, no altered final tolerance.
// A frozen pressure projection is never accepted as a physical mass certificate.
// No production-commit mechanism is exposed by this audit: the authoritative
// momentum linear and SST transport certificates must be reissued by the
// production driver before any future commit integration.
#include "cfdx/physics/bfs_joint_momentum_sst_v71.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct BfsJointChoiceV72 {
    double momentum_fraction=0.;
    double pressure_fraction=0.;
};

struct BfsSixFieldPhysicalCertificateV72 {
    std::array<double,6> signed_linf{};
    std::array<double,6> scaled_linf{};
    std::array<bool,6> physical_block_closed{};
    double mass_divergence_linf=std::numeric_limits<double>::infinity();
    double mass_boundary_signed=0.;
    double face_owner_neighbour_identity=std::numeric_limits<double>::infinity();
    double signed_continuity_vs_face_gap=std::numeric_limits<double>::infinity();
    double rhie_chow_fixed_point_gap=std::numeric_limits<double>::infinity();
    bool replay_identical=false;
    bool original_v65_guards_ok=false;
    bool all_six_physical_blocks_closed=false;
    bool conservation_reassembled_closed=false;
    // Linear momentum is only diagonal-Picard-predicted in V71, not re-solved.
    // Thus no fresh production linear-momentum certificate is available here.
    bool native_momentum_linear_recertified=false;
    // The real SST transport is called only on a private state. It is not a
    // newly issued authoritative production SST certificate after commit.
    bool native_sst_transport_recertified=false;
    bool eligible_for_production_commit=false;
};

struct BfsJointFixedPointAttemptV72 {
    std::size_t cycle=0;
    double momentum_fraction=0.,pressure_fraction=0.;
    double physical_mass_before=0.,physical_mass_after=0.;
    double physical_omega_after=0.,merit_ratio=0.,worst_growth=0.;
    double face_flux_change_linf=0.,rAU_change_linf=0.,nut_change_linf=0.;
    bool evaluated=false,guards_ok=false,mass_improved=false,selected=false;
    std::string rejection;
};

struct BfsJointFixedPointV72 {
    std::vector<BfsJointFixedPointAttemptV72> attempts;
    std::vector<BfsSixFieldPhysicalCertificateV72> certificates;
    std::size_t attempted_cycles=0,accepted_private_cycles=0;
    bool physical_mass_closed=false,all_six_physical_blocks_closed=false;
    bool stalled=false,commit_authorized=false;
    double initial_mass=0.,final_mass=0.;
    std::string status;
};

// A new certificate is constructed from the SAME final six-field candidate,
// a new independent F_phys evaluation, and an owner/neighbour face audit.
// The normalized entries are diagnostics, never used in a derivative.
// At no point is this the (missing) native momentum linear-solver certificate.
template<class Evaluate,class Owner,class Neighbour>
BfsSixFieldPhysicalCertificateV72 bfs_recertify_six_field_v72(
    const PhysicalCandidateV51& old_state,const PhysicalResidualV51& old_result,
    const PhysicalCandidateV51& candidate,const PhysicalResidualV51& final_result,
    Evaluate&& evaluate,Owner&& owner,Neighbour&& neighbour,
    double relative_tolerance,double continuity_tolerance,
    double absolute_mass_tolerance=1e-10,
    double k_min=0.,double omega_min=0.) {
    if(!(relative_tolerance>0.) || !(continuity_tolerance>0.) ||
       !(absolute_mass_tolerance>0.) || !std::isfinite(absolute_mass_tolerance))
        throw std::invalid_argument("V72 invalid certification tolerances");
    validate_physical_result_v51(candidate,final_result);
    const auto independent=evaluate_physical_candidate_v51(
        candidate,evaluate,k_min,omega_min);
    verify_candidate_replay_v51(final_result,independent);
    BfsSixFieldPhysicalCertificateV72 cert;
    cert.replay_identical=true;
    const auto mass=bfs_mass_balance_v70(candidate.cells,final_result.face_flux,
                                         owner,neighbour);
    cert.mass_divergence_linf=mass.divergence_linf;
    cert.mass_boundary_signed=mass.boundary_signed;
    cert.face_owner_neighbour_identity=mass.identity_error;
    std::vector<long double> div(candidate.cells,0.);
    for(std::size_t f=0;f<final_result.face_flux.size();++f) {
        const auto a=static_cast<std::size_t>(owner(f));
        const auto b=neighbour(f);
        const long double flux=final_result.face_flux[f];
        div[a]+=flux;
        if(b>=0)div[static_cast<std::size_t>(b)]-=flux;
    }
    for(std::size_t i=0;i<candidate.cells;++i)
        cert.signed_continuity_vs_face_gap=std::max(
            i==0?0.:cert.signed_continuity_vs_face_gap,
            std::abs(static_cast<double>(div[i])-final_result.signed_raw[3*candidate.cells+i]));
    if(cert.signed_continuity_vs_face_gap>absolute_mass_tolerance ||
       mass.identity_error>absolute_mass_tolerance)
        throw std::runtime_error("V72 signed FV continuity disagrees with owner/neighbour flux");
    cert.rhie_chow_fixed_point_gap=final_result.flux_fixed_point_defect;
    const auto merit=evaluate_comparable_merit_v65(old_state,candidate,
         old_result,final_result,relative_tolerance,continuity_tolerance,false,false);
    cert.original_v65_guards_ok=merit.accepted;
    cert.all_six_physical_blocks_closed=true;
    for(std::size_t b=0;b<6;++b) {
        for(std::size_t i=0;i<candidate.cells;++i) {
            const std::size_t j=b*candidate.cells+i;
            cert.signed_linf[b]=std::max(cert.signed_linf[b],std::abs(final_result.signed_raw[j]));
            cert.scaled_linf[b]=std::max(cert.scaled_linf[b],std::abs(final_result.normalized[j]));
        }
        cert.physical_block_closed[b]=cert.scaled_linf[b]<=
            (b==3?continuity_tolerance:relative_tolerance);
        cert.all_six_physical_blocks_closed &=cert.physical_block_closed[b];
    }
    cert.conservation_reassembled_closed=
        mass.divergence_linf<=absolute_mass_tolerance &&
        mass.identity_error<=absolute_mass_tolerance &&
        cert.signed_continuity_vs_face_gap<=absolute_mass_tolerance &&
        cert.rhie_chow_fixed_point_gap<=absolute_mass_tolerance;
    // Deliberately fail-closed: this private diagonal Picard probe never owns
    // the native linear momentum/SST certificates or the authoritative phi.
    cert.eligible_for_production_commit=cert.original_v65_guards_ok &&
        cert.conservation_reassembled_closed && cert.all_six_physical_blocks_closed &&
        cert.native_momentum_linear_recertified &&
        cert.native_sst_transport_recertified;
    return cert;
}

// Re-uses the V71 predictor -> pressure+flux -> real SST transport -> F_phys
// cycle. A candidate is retained privately only if (i) the unchanged V65
// two guards pass, and (ii) its REASSEMBLED mass divergence strictly decreases.
// All candidate families start from the SAME private state each outer cycle.
// Intermediate states do not become production states, not even when mass is
// closed: a fresh momentum linear solve and production SST certificate are
// not available inside this read-only audit.
template<class Evaluate,class Project,class Transport,class Owner,class Neighbour>
BfsJointFixedPointV72 bfs_joint_private_fixed_point_v72(
    const PhysicalCandidateV51& old_state,const PhysicalResidualV51& old_f,
    const PhysicalCandidateV51& start,const PhysicalResidualV51& start_f,
    Evaluate&& evaluate,Project&& project,Transport&& transport,
    Owner&& owner,Neighbour&& neighbour,
    const std::vector<BfsJointChoiceV72>& choices,
    std::size_t max_cycles,double k_min,double omega_min,
    double relative_tolerance,double continuity_tolerance,
    double final_mass_tolerance=1e-10,
    PhysicalCandidateV51* last_private_state=nullptr,
    PhysicalResidualV51* last_private_residual=nullptr) {
    if(start.cells!=old_state.cells || start.cells>152 || choices.empty() ||
       choices.size()>8 || max_cycles==0 || max_cycles>16 ||
       !(final_mass_tolerance>0.) || !std::isfinite(final_mass_tolerance))
        throw std::invalid_argument("V72 pilot fixed-point contract invalid");
    validate_physical_result_v51(old_state,old_f);
    validate_physical_result_v51(start,start_f);
    auto current=start;
    auto current_f=start_f;
    BfsJointFixedPointV72 result;
    result.initial_mass=bfs_mass_balance_v70(start.cells,start_f.face_flux,
                          owner,neighbour).divergence_linf;
    result.final_mass=result.initial_mass;
    for(std::size_t cycle=1;cycle<=max_cycles;++cycle) {
        result.attempted_cycles=cycle;
        const auto snapshot=current;
        const auto f_snapshot=current_f;
        bool found=false;
        double best_mass=std::numeric_limits<double>::infinity();
        PhysicalCandidateV51 selected_state;
        PhysicalResidualV51 selected_f;
        BfsSixFieldPhysicalCertificateV72 selected_cert;
        std::size_t selected_index=0;
        for(const auto& choice:choices) {
            BfsJointFixedPointAttemptV72 evidence;
            evidence.cycle=cycle;
            evidence.momentum_fraction=choice.momentum_fraction;
            evidence.pressure_fraction=choice.pressure_fraction;
            evidence.physical_mass_before=result.final_mass;
            try {
                PhysicalCandidateV51 proposed;
                PhysicalResidualV51 reassembled;
                const auto step=bfs_joint_private_cycle_v71(
                    old_state,old_f,current,current_f,evaluate,project,transport,
                    owner,neighbour,choice.momentum_fraction,choice.pressure_fraction,
                    k_min,omega_min,relative_tolerance,continuity_tolerance,
                    final_mass_tolerance,cycle,&proposed,&reassembled);
                if(current.fields!=snapshot.fields ||
                   current_f.face_flux!=f_snapshot.face_flux ||
                   current_f.signed_raw!=f_snapshot.signed_raw)
                    throw std::runtime_error("V72 V71 operator mutated the private source state");
                const auto cert=bfs_recertify_six_field_v72(
                    old_state,old_f,proposed,reassembled,evaluate,owner,neighbour,
                    relative_tolerance,continuity_tolerance,final_mass_tolerance,
                    k_min,omega_min);
                evidence.evaluated=true;
                evidence.physical_mass_after=cert.mass_divergence_linf;
                evidence.physical_omega_after=step.omega_after;
                evidence.merit_ratio=step.merit_ratio;
                evidence.worst_growth=step.worst_growth;
                evidence.guards_ok=step.original_guards_ok && cert.original_v65_guards_ok;
                evidence.mass_improved=cert.mass_divergence_linf<=final_mass_tolerance ||
                    cert.mass_divergence_linf<result.final_mass*(1.-1e-12);
                evidence.rAU_change_linf=step.rAU_change_linf;
                evidence.nut_change_linf=step.turbulent_viscosity_change_linf;
                for(std::size_t f=0;f<current_f.face_flux.size();++f)
                    evidence.face_flux_change_linf=std::max(evidence.face_flux_change_linf,
                        std::abs(reassembled.face_flux[f]-current_f.face_flux[f]));
                if(!evidence.guards_ok) evidence.rejection="V65_original_guards";
                else if(!evidence.mass_improved) evidence.rejection="physical_mass_not_decreasing";
                else if(cert.mass_divergence_linf<best_mass) {
                    found=true;
                    best_mass=cert.mass_divergence_linf;
                    selected_state=std::move(proposed);
                    selected_f=std::move(reassembled);
                    selected_cert=cert;
                    selected_index=result.attempts.size();
                }
            } catch(const std::exception& error) {
                evidence.rejection=error.what();
            }
            result.attempts.push_back(std::move(evidence));
        }
        if(!found) {
            result.stalled=true;
            result.status="stalled_no_trial_satisfies_both_original_guards_and_physical_mass";
            break;
        }
        result.attempts[selected_index].selected=true;
        ++result.accepted_private_cycles;
        current=std::move(selected_state);
        current_f=std::move(selected_f);
        result.certificates.push_back(selected_cert);
        result.final_mass=selected_cert.mass_divergence_linf;
        result.physical_mass_closed=selected_cert.conservation_reassembled_closed;
        result.all_six_physical_blocks_closed=selected_cert.all_six_physical_blocks_closed;
        if(result.physical_mass_closed && result.all_six_physical_blocks_closed) {
            result.status="physical_blocks_closed_but_native_linear_certificates_missing";
            break;
        }
    }
    if(result.status.empty())result.status="cycle_budget_exhausted_without_full_certificate";
    result.commit_authorized=false;
    if(last_private_state)*last_private_state=current;
    if(last_private_residual)*last_private_residual=current_f;
    return result;
}
} // namespace cfdx::physics