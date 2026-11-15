#pragma once
// V64: READ-ONLY 2x2 counterfactual decomposition on the same V51 Fphys.
// The four corners are F(flow_before/after, SST_before/after). Each evaluation
// rebuilds RC, wall terms, mu_t, rAU and FV transport from its own candidate.
// This is NOT a frozen-system solve, a Jacobian, or an accepted nonlinear step.
// Never use the mixed-state results to modify production state or tolerances.

#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

struct BfsCausalBlockV64 {
    // Corner order: 00=before, 10=flow-only, 01=SST-only, 11=after.
    std::array<double,4> raw_linf{};
    double flow_effect_linf=0.0;   // F10-F00, U/p update at SST_before
    double sst_effect_linf=0.0;    // F01-F00, SST update at flow_before
    double interaction_linf=0.0;   // F11-F10-F01+F00
    double interaction_relative=0.0; // anchored to max raw corner amplitude
    std::size_t after_hotspot_cell=0;
    std::size_t interaction_hotspot_cell=0;
};

struct BfsCausalAblationV64 {
    std::size_t cells=0;
    bool flow_changed=false;
    bool sst_changed=false;
    double relative_flow_state_change=0.0;
    double relative_sst_state_change=0.0;
    std::array<BfsCausalBlockV64,6> blocks{};
    std::array<double,4> rc_fixed_point_defect{};
    // Closure identity is an arithmetic cross-check, not physics validation.
    double max_recomposition_error=0.0;
    // Changes in the RE-CLOSED RC mass flux and state-derived turbulent mu.
    double rc_flux_flow_effect=0.0;
    double rc_flux_sst_effect=0.0;
    double rc_flux_interaction=0.0;
    std::size_t rc_flux_interaction_face=0;
    double eddy_viscosity_flow_effect=0.0;
    double eddy_viscosity_sst_effect=0.0;
    double eddy_viscosity_interaction=0.0;
    std::size_t eddy_viscosity_interaction_cell=0;
};

template<class Evaluate>
BfsCausalAblationV64 audit_bfs_causal_ablation_v64(
    const PhysicalCandidateV51& before,
    const PhysicalCandidateV51& after,
    const PhysicalResidualV51& after_residual,
    Evaluate&& physical_evaluator,
    double k_min=0.0, double omega_min=0.0)
{
    validate_physical_candidate_v51(before,k_min,omega_min);
    validate_physical_candidate_v51(after,k_min,omega_min);
    if (before.cells!=after.cells)
        throw std::invalid_argument("V64 before/after cell count mismatch");
    validate_physical_result_v51(after,after_residual);
    // Prove that the physical callback did not silently change since the
    // already measured F11. A stateful/history-dependent residual is invalid.
    verify_candidate_replay_v51(after_residual,
        evaluate_physical_candidate_v51(after,physical_evaluator,k_min,omega_min));
    const std::size_t n=before.cells;
    PhysicalCandidateV51 flow_only=before, sst_only=before;
    BfsCausalAblationV64 report;
    report.cells=n;
    for (std::size_t b=0;b<6;++b) {
        for (std::size_t i=0;i<n;++i) {
            const std::size_t index=b*n+i;
            const double oldv=before.fields[index], newv=after.fields[index];
            if (b<4)flow_only.fields[index]=newv;
            else sst_only.fields[index]=newv;
            if (oldv!=newv) {
                const double relative=std::abs(newv-oldv)/
                    std::max({1.0,std::abs(newv),std::abs(oldv)});
                if (b<4) {
                    report.flow_changed=true;
                    report.relative_flow_state_change=
                        std::max(report.relative_flow_state_change,relative);
                } else {
                    report.sst_changed=true;
                    report.relative_sst_state_change=
                        std::max(report.relative_sst_state_change,relative);
                }
            }
        }
    }
    // The callback MUST be the same V51 physical operator as for `after`.
    const std::array<PhysicalResidualV51,4> residuals{
        evaluate_physical_candidate_v51(before,physical_evaluator,k_min,omega_min),
        evaluate_physical_candidate_v51(flow_only,physical_evaluator,k_min,omega_min),
        evaluate_physical_candidate_v51(sst_only,physical_evaluator,k_min,omega_min),
        after_residual
    };
    for (std::size_t corner=0;corner<4;++corner)
        report.rc_fixed_point_defect[corner]=residuals[corner].flux_fixed_point_defect;
    const auto nfaces=residuals[0].face_flux.size();
    for (std::size_t corner=1;corner<4;++corner)
        if (residuals[corner].face_flux.size()!=nfaces)
            throw std::runtime_error("V64 inconsistent Rhie-Chow face counts");
    for (std::size_t f=0;f<nfaces;++f) {
        const double a=residuals[0].face_flux[f];
        const double u=residuals[1].face_flux[f];
        const double t=residuals[2].face_flux[f];
        const double z=residuals[3].face_flux[f];
        report.rc_flux_flow_effect=std::max(report.rc_flux_flow_effect,std::abs(u-a));
        report.rc_flux_sst_effect=std::max(report.rc_flux_sst_effect,std::abs(t-a));
        const double interaction=std::abs((z-u)-(t-a));
        if(interaction>report.rc_flux_interaction){
            report.rc_flux_interaction=interaction;
            report.rc_flux_interaction_face=f;
        }
    }
    for(std::size_t i=0;i<n;++i){
        const double a=residuals[0].eddy_viscosity[i];
        const double u=residuals[1].eddy_viscosity[i];
        const double t=residuals[2].eddy_viscosity[i];
        const double z=residuals[3].eddy_viscosity[i];
        report.eddy_viscosity_flow_effect=std::max(
            report.eddy_viscosity_flow_effect,std::abs(u-a));
        report.eddy_viscosity_sst_effect=std::max(
            report.eddy_viscosity_sst_effect,std::abs(t-a));
        const double interaction=std::abs((z-u)-(t-a));
        if(interaction>report.eddy_viscosity_interaction){
            report.eddy_viscosity_interaction=interaction;
            report.eddy_viscosity_interaction_cell=i;
        }
    }
    for (std::size_t b=0;b<6;++b) {
        auto& block=report.blocks[b];
        for (std::size_t i=0;i<n;++i) {
            const std::size_t index=b*n+i;
            const double a=residuals[0].signed_raw[index];
            const double f=residuals[1].signed_raw[index];
            const double s=residuals[2].signed_raw[index];
            const double z=residuals[3].signed_raw[index];
            const std::array<double,4> corners{a,f,s,z};
            for (std::size_t c=0;c<4;++c)
                block.raw_linf[c]=std::max(block.raw_linf[c],std::abs(corners[c]));
            block.flow_effect_linf=std::max(block.flow_effect_linf,std::abs(f-a));
            block.sst_effect_linf=std::max(block.sst_effect_linf,std::abs(s-a));
            const double interaction=(z-f)-(s-a);
            if (std::abs(interaction)>block.interaction_linf) {
                block.interaction_linf=std::abs(interaction);
                block.interaction_hotspot_cell=i;
            }
            if (std::abs(z)>std::abs(residuals[3].signed_raw[b*n+block.after_hotspot_cell]))
                block.after_hotspot_cell=i;
            // Compare two equivalent decompositions of the same four corners.
            // This detects implementation errors without claiming FV validity.
            const double recomposed=(a+(f-a)) + ((s-a)+interaction);
            report.max_recomposition_error=std::max(report.max_recomposition_error,
                std::abs(recomposed-z));
        }
        const double reference=std::max({block.raw_linf[0],block.raw_linf[1],
            block.raw_linf[2],block.raw_linf[3],1e-300});
        block.interaction_relative=block.interaction_linf/reference;
        if (!std::isfinite(block.interaction_relative) ||
            !std::isfinite(block.flow_effect_linf) ||
            !std::isfinite(block.sst_effect_linf))
            throw std::runtime_error("V64 nonfinite causal ablation measurement");
    }
    return report;
}

} // namespace cfdx::physics
