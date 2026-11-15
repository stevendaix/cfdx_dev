#pragma once
// V71: strictly private, candidate-owned momentum/pressure/RC/SST co-iteration.
// Uses the production F_phys, V70 pressure response, and production SST transport.
// The diagonal momentum predictor is a Picard *diagnostic*, not a Newton step.
// There is deliberately NO commit path and NO change to physical tolerances.
#include "cfdx/physics/bfs_sst_rc_conservative_v70.h"
#include "cfdx/physics/bfs_sst_flux_picard_v69.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct BfsJointPicardStepV71 {
    std::size_t cycle=0;
    double predictor_fraction=0.;
    double momentum_update_linf=0.;
    double rAU_change_linf=0.;
    double turbulent_viscosity_change_linf=0.;
    double mass_before=0.;
    double mass_after_predictor=0.;
    double mass_after_pressure=0.;
    double mass_after_sst=0.;
    double omega_before=0.;
    double omega_after=0.;
    double merit_ratio=0.;
    double worst_growth=0.;
    std::size_t worst_block=0;
    bool joint_ok=false,monotone_ok=false, original_guards_ok=false;
    bool final_mass_closed=false;
    bool private_only=true;
};

// Consistent with the volume-integrated FV momentum system A*U=b:
// U <- U - gamma * diag(A)^{-1} (A*U - b). The rAU and R used for
// this direction MUST come from the same reassembled state. In a 2-D
// extrusion with Uz identically zero the inactive component is not changed.
inline PhysicalCandidateV51 bfs_momentum_predict_v71(
    const PhysicalCandidateV51& state,const PhysicalResidualV51& physical,
    double gamma,double* update_linf=nullptr) {
    validate_physical_candidate_v51(state);
    validate_physical_result_v51(state,physical);
    const std::size_t n=state.cells;
    if(physical.inverse_momentum_diagonal.size()!=3*n ||
       !(gamma>=0.) || !(gamma<=1.) || !std::isfinite(gamma))
        throw std::invalid_argument("V71 invalid diagonal momentum predictor");
    auto trial=state;
    bool extruded_2d=true;
    for(std::size_t i=0;i<n;++i)
        if(state.fields[2*n+i]!=0.)extruded_2d=false;
    double largest=0.;
    for(std::size_t d=0;d<3;++d) {
        if(d==2 && extruded_2d)continue;
        for(std::size_t i=0;i<n;++i) {
            const auto j=d*n+i;
            const double inverse=physical.inverse_momentum_diagonal[j];
            if(!(inverse>0.) || !std::isfinite(inverse))
                throw std::runtime_error("V71 stale or invalid rAU");
            const double delta=-gamma*inverse*physical.signed_raw[j];
            if(!std::isfinite(delta))
                throw std::runtime_error("V71 nonfinite momentum predictor");
            trial.fields[j]+=delta;
            largest=std::max(largest,std::abs(delta));
        }
    }
    validate_physical_candidate_v51(trial);
    if(update_linf)*update_linf=largest;
    return trial;
}

// One existing-operator co-iteration on a PRIVATE copy:
// momentum diagonal Picard -> candidate RC/rAU -> coupled U/p pressure response
// -> true RC/rAU -> SST transport on THAT flux -> true six-field FV residual.
// No projection on frozen coefficients is ever called a final certificate.
// The V65 merit is computed against the SAME old-state normalization.
// The projected response must preserve the candidate's non-pressure fields
// exactly, except the physical U correction caused by delta(p).
// Returns only evidence. In particular, an eligible candidate is NOT committed.
template<class Evaluate,class Project,class Transport,class Owner,class Neighbour>
BfsJointPicardStepV71 bfs_joint_private_cycle_v71(
    const PhysicalCandidateV51& old_state,const PhysicalResidualV51& old_f,
    const PhysicalCandidateV51& start,const PhysicalResidualV51& start_f,
    Evaluate&& evaluate,Project&& project,Transport&& transport,
    Owner&& owner,Neighbour&& neighbour,
    double gamma,double projection_fraction,double k_min,double omega_min,
    double relative_tolerance,double continuity_tolerance,
    double final_mass_tolerance=1e-10,std::size_t cycle=1,
    PhysicalCandidateV51* private_final_state=nullptr,
    PhysicalResidualV51* private_final_residual=nullptr) {
    validate_physical_result_v51(old_state,old_f);
    validate_physical_result_v51(start,start_f);
    if(old_state.cells!=start.cells || !(projection_fraction>0.) ||
       projection_fraction>1. || !(final_mass_tolerance>0.) ||
       !std::isfinite(final_mass_tolerance))
        throw std::invalid_argument("V71 inconsistent private cycle");
    BfsJointPicardStepV71 row;
    row.cycle=cycle;
    row.predictor_fraction=gamma;
    row.mass_before=bfs_mass_balance_v70(start.cells,start_f.face_flux,owner,neighbour).divergence_linf;
    row.omega_before=bfs_sst_block_linf_v69(start_f,start.cells,5);
    auto predictor=bfs_momentum_predict_v71(start,start_f,gamma,&row.momentum_update_linf);
    const auto f_predict=evaluate_physical_candidate_v51(predictor,evaluate,k_min,omega_min);
    row.mass_after_predictor=bfs_mass_balance_v70(start.cells,f_predict.face_flux,owner,neighbour).divergence_linf;
    const auto projected=project(predictor,f_predict,projection_fraction);
    validate_physical_candidate_v51(projected.state,k_min,omega_min);
    if(projected.state.cells!=start.cells ||
       projected.face_flux.size()!=f_predict.face_flux.size())
        throw std::runtime_error("V71 pressure response incompatible with physical state");
    // The pressure-velocity substep is not allowed to solve turbulence on the
    // side; k/omega must be changed by the subsequent *real* SST transport.
    const auto n=start.cells;
    if(!std::equal(predictor.fields.begin()+static_cast<std::ptrdiff_t>(4*n),
                   predictor.fields.end(),
                   projected.state.fields.begin()+static_cast<std::ptrdiff_t>(4*n)))
        throw std::runtime_error("V71 pressure response modified k/omega");
    // The same V70 operator must close pressure/velocity and face flux.
    const auto fixed=bfs_mass_balance_v70(start.cells,projected.face_flux,owner,neighbour);
    row.mass_after_pressure=fixed.divergence_linf;
    if(fixed.identity_error>final_mass_tolerance ||
       (projected.fully_projected && fixed.divergence_linf>final_mass_tolerance))
        throw std::runtime_error("V71 frozen pressure response does not conserve mass");
    const auto f_pressure=evaluate_physical_candidate_v51(projected.state,evaluate,k_min,omega_min);
    const auto after_sst=transport(projected.state,f_pressure.face_flux);
    validate_physical_candidate_v51(after_sst,k_min,omega_min);
    if(after_sst.cells!=start.cells ||
       !std::equal(projected.state.fields.begin(),
                   projected.state.fields.begin()+static_cast<std::ptrdiff_t>(4*start.cells),
                   after_sst.fields.begin()))
        throw std::runtime_error("V71 SST transport changed U or p");
    const auto final=evaluate_physical_candidate_v51(after_sst,evaluate,k_min,omega_min);
    const auto replay=evaluate_physical_candidate_v51(after_sst,evaluate,k_min,omega_min);
    verify_candidate_replay_v51(final,replay);
    row.mass_after_sst=bfs_mass_balance_v70(start.cells,final.face_flux,owner,neighbour).divergence_linf;
    row.omega_after=bfs_sst_block_linf_v69(final,start.cells,5);
    const auto gate=evaluate_comparable_merit_v65(
        old_state,after_sst,old_f,final,relative_tolerance,continuity_tolerance,false,false);
    row.merit_ratio=gate.merit_ratio;
    row.worst_growth=gate.largest_growth;
    row.worst_block=gate.worst_block;
    row.joint_ok=gate.joint_allowed;
    row.monotone_ok=gate.monotone_allowed;
    row.original_guards_ok=gate.accepted;
    row.final_mass_closed=row.mass_after_sst<=final_mass_tolerance;
    if(start_f.inverse_momentum_diagonal.size()!=final.inverse_momentum_diagonal.size() ||
       start_f.eddy_viscosity.size()!=final.eddy_viscosity.size())
        throw std::runtime_error("V71 momentum or SST coefficients incomplete");
    for(std::size_t i=0;i<start_f.inverse_momentum_diagonal.size();++i)
        row.rAU_change_linf=std::max(row.rAU_change_linf,
            std::abs(final.inverse_momentum_diagonal[i]-start_f.inverse_momentum_diagonal[i]));
    for(std::size_t i=0;i<start_f.eddy_viscosity.size();++i)
        row.turbulent_viscosity_change_linf=std::max(row.turbulent_viscosity_change_linf,
            std::abs(final.eddy_viscosity[i]-start_f.eddy_viscosity[i]));
    if(!std::isfinite(row.merit_ratio) || !std::isfinite(row.worst_growth) ||
       !std::isfinite(row.omega_after) || !std::isfinite(row.mass_after_sst))
        throw std::runtime_error("V71 nonfinite joint evidence");
    if(private_final_state)*private_final_state=after_sst;
    if(private_final_residual)*private_final_residual=final;
    return row;
}
} // namespace cfdx::physics
