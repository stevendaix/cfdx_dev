@@ -0,0 +1,249 @@
#pragma once
// V66: bounded, read-only diagnosis of the V65 BFS152 step rejection.
// No change to F_phys, V13/V14 merit, Rhie-Chow, SST, or nonlinear updates.
// Cell/face IDs refer to the ACTUAL mesh ownership and native state-only flux.
#include "cfdx/physics/comparable_physical_merit_v65.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct V66FluxCarrierGap {
    double relative_l1=0.;
    double max_abs=0.;
    double carrier_divergence_gap_linf=0.;
    std::size_t worst_face=0;
    std::size_t worst_cell=0;
};

struct V66CoefficientGap {
    double relative_linf=0.;
    double absolute_linf=0.;
    std::size_t worst_component=0,worst_cell=0;
};

// Compare two array layouts [component][cell] without rewriting either.
inline V66CoefficientGap v66_coefficient_gap(
    const std::vector<double>& state_only,
    const std::vector<double>& carrier,std::size_t cells)
{
    if (!cells || state_only.size()!=carrier.size() ||
        state_only.size()%cells!=0 || state_only.empty())
        throw std::invalid_argument("V66 coefficient layout mismatch");
    V66CoefficientGap out;
    for(std::size_t j=0;j<carrier.size();++j) {
        const double a=state_only[j],b=carrier[j];
        if (!std::isfinite(a) || !std::isfinite(b) || a<0. || b<0.)
            throw std::invalid_argument("V66 invalid positive coefficient");
        const double abs=std::abs(b-a);
        const double rel=abs/std::max({std::abs(a),std::abs(b),1e-300});
        out.absolute_linf=std::max(out.absolute_linf,abs);
        if (rel>out.relative_linf) {
            out.relative_linf=rel;
            out.worst_component=j/cells;
            out.worst_cell=j%cells;
        }
    }
    return out;
}

struct V66OmegaCell {
    std::size_t cell=0;
    double old_signed=0.;
    double trial_signed=0.;
    double delta_signed=0.;
    double omega_old=0.;
    double omega_trial=0.;
    double k_old=0.;
    double k_trial=0.;
    double nut_old=0.;
    double nut_trial=0.;
};

struct V66BfsCausalAudit {
    V66FluxCarrierGap old_carrier,trial_carrier;
    double rau_relative_linf=0.;
    std::size_t rau_worst_cell=0,rau_worst_component=0;
    double nut_relative_linf=0.;
    std::size_t nut_worst_cell=0;
    V66OmegaCell omega_worst_change,omega_worst_trial,omega_worst_old;
    bool valid=false;
};

// The component-wise relative change reports the normalization explicitly:
// max(|a-b|) / max(max(|a|), max(|b|), 1e-300). It is NOT a merit norm.
inline double v66_relative_change(double a,double b) {
    if (!std::isfinite(a) || !std::isfinite(b))
        throw std::invalid_argument("V66 nonfinite coefficient");
    return std::abs(b-a)/std::max({std::abs(a),std::abs(b),1e-300});
}

template<class Owner,class Neighbour>
V66FluxCarrierGap v66_carrier_gap(
    const PhysicalCandidateV51& state,const PhysicalResidualV51& physical,
    const std::vector<double>& carrier,Owner&& owner,Neighbour&& neighbour)
{
    validate_physical_result_v51(state,physical);
    if (carrier.size()!=physical.face_flux.size())
        throw std::invalid_argument("V66 carrier flux size mismatch");
    std::vector<long double> discrepancy(state.cells,0.L);
    long double l1=0.,ref=0.;
    V66FluxCarrierGap out;
    for (std::size_t f=0;f<carrier.size();++f) {
        const auto o=static_cast<std::size_t>(owner(f));
        const auto nb=neighbour(f);
        if (o>=state.cells || nb < -1 ||
            (nb>=0 && static_cast<std::size_t>(nb)>=state.cells) ||
            !std::isfinite(carrier[f]))
            throw std::invalid_argument("V66 invalid flux or face ownership");
        const double err=carrier[f]-physical.face_flux[f];
        if (std::abs(err)>out.max_abs) {
            out.max_abs=std::abs(err);
            out.worst_face=f;
        }
        l1+=std::abs(static_cast<long double>(err));
        ref+=std::abs(static_cast<long double>(physical.face_flux[f]));
        discrepancy[o]+=static_cast<long double>(err);
        if (nb>=0) discrepancy[static_cast<std::size_t>(nb)]-=static_cast<long double>(err);
    }
    out.relative_l1=static_cast<double>(l1/std::max(1.e-300L,ref));
    for(std::size_t c=0;c<state.cells;++c) {
        const double diff=std::abs(static_cast<double>(discrepancy[c]));
        if (diff>out.carrier_divergence_gap_linf) {
            out.carrier_divergence_gap_linf=diff;
            out.worst_cell=c;
        }
    }
    if (!std::isfinite(out.relative_l1) ||
        !std::isfinite(out.carrier_divergence_gap_linf))
        throw std::runtime_error("V66 nonfinite carrier audit");
    return out;
}

template<class Owner,class Neighbour>
V66BfsCausalAudit audit_bfs_flux_sst_v66(
    const PhysicalCandidateV51& old_state,
    const PhysicalCandidateV51& trial_state,
    const PhysicalResidualV51& old_physical,
    const PhysicalResidualV51& trial_physical,
    const std::vector<double>& old_carrier,
    const std::vector<double>& trial_carrier,
    Owner&& owner,Neighbour&& neighbour)
{
    validate_physical_candidate_v51(old_state);
    validate_physical_candidate_v51(trial_state);
    validate_physical_result_v51(old_state,old_physical);
    validate_physical_result_v51(trial_state,trial_physical);
    const auto n=old_state.cells;
    if (n==0 || trial_state.cells!=n ||
        old_physical.inverse_momentum_diagonal.size()!=3*n ||
        trial_physical.inverse_momentum_diagonal.size()!=3*n ||
        old_carrier.size()!=trial_carrier.size())
        throw std::invalid_argument("V66 incomplete or incompatible candidate rAU / states");
    V66BfsCausalAudit out;
    out.old_carrier=v66_carrier_gap(old_state,old_physical,old_carrier,owner,neighbour);
    out.trial_carrier=v66_carrier_gap(trial_state,trial_physical,trial_carrier,owner,neighbour);
    for (std::size_t j=0;j<3*n;++j) {
        const auto change=v66_relative_change(
            old_physical.inverse_momentum_diagonal[j],
            trial_physical.inverse_momentum_diagonal[j]);
        if (change>out.rau_relative_linf) {
            out.rau_relative_linf=change;
            out.rau_worst_component=j/n;
            out.rau_worst_cell=j%n;
        }
    }
    std::size_t worst_delta=0,worst_trial=0,worst_old=0;
    for(std::size_t c=0;c<n;++c) {
        const auto change=v66_relative_change(
            old_physical.eddy_viscosity[c],trial_physical.eddy_viscosity[c]);
        if (change>out.nut_relative_linf) {
            out.nut_relative_linf=change;
            out.nut_worst_cell=c;
        }
        const auto old_w=old_physical.signed_raw[5*n+c];
        const auto trial_w=trial_physical.signed_raw[5*n+c];
        if (std::abs(trial_w-old_w) > std::abs(
                trial_physical.signed_raw[5*n+worst_delta]-
                old_physical.signed_raw[5*n+worst_delta])) worst_delta=c;
        if (std::abs(trial_w)>std::abs(trial_physical.signed_raw[5*n+worst_trial]))
            worst_trial=c;
        if (std::abs(old_w)>std::abs(old_physical.signed_raw[5*n+worst_old]))
            worst_old=c;
    }
    const auto make_cell=[&](std::size_t c) {
        V66OmegaCell row;
        row.cell=c;
        row.old_signed=old_physical.signed_raw[5*n+c];
        row.trial_signed=trial_physical.signed_raw[5*n+c];
        row.delta_signed=row.trial_signed-row.old_signed;
        row.omega_old=old_state.fields[5*n+c];
        row.omega_trial=trial_state.fields[5*n+c];
        row.k_old=old_state.fields[4*n+c];
        row.k_trial=trial_state.fields[4*n+c];
        row.nut_old=old_physical.eddy_viscosity[c];
        row.nut_trial=trial_physical.eddy_viscosity[c];
        return row;
    };
    out.omega_worst_change=make_cell(worst_delta);
    out.omega_worst_trial=make_cell(worst_trial);
    out.omega_worst_old=make_cell(worst_old);
    out.valid=true;
    return out;
}

// The SST row budget is obtained by reassembling the *same* state-only FV
// transport at a selected cell. These are SIGNED row units, not normalized
// V65 block-merit units. The five terms include explicit cross-diffusion and
// other sources in the production sign split; no physics source is modified.
struct V66OmegaBudgetParts {
    double convection=0., diffusion=0., physical_production=0.;
    double physical_destruction=0., physical_cross=0.;
    double full_residual=0., decomposition_error=0.;
    double source_identity_error=0., wall_flux=0.;
    bool constraint=false;
};
struct V66OmegaBudgetDifference {
    double convection=0.,diffusion=0.,production=0.;
    double destruction=0.,cross=0.,sum=0.,residual=0.,closure=0.;
    bool constraint=false;
};
inline V66OmegaBudgetDifference v66_compare_omega_budgets(
    const V66OmegaBudgetParts& old,const V66OmegaBudgetParts& trial)
{
    const auto valid=[](const V66OmegaBudgetParts& b) {
        const double a[]={b.convection,b.diffusion,b.physical_production,
            b.physical_destruction,b.physical_cross,b.full_residual,
            b.decomposition_error,b.source_identity_error,b.wall_flux};
        for (double v:a) if (!std::isfinite(v))
            throw std::invalid_argument("V66 nonfinite SST row budget");
        const double sum=b.convection+b.diffusion+b.physical_production+
            b.physical_destruction+b.physical_cross;
        // Relative to the terms to account for conservative cancellation.
        const double scale=std::max({1.,std::abs(sum),std::abs(b.full_residual)});
        if (std::abs(sum-b.full_residual)>1e-8*scale ||
            std::abs(b.decomposition_error)>1e-8*scale ||
            std::abs(b.source_identity_error)>1e-8*scale)
            throw std::runtime_error("V66 SST budget is not a closed FV row");
    };
    valid(old);valid(trial);
    if (old.constraint!=trial.constraint)
        throw std::runtime_error("V66 SST owner-row branch changed");
    V66OmegaBudgetDifference d;
    d.convection=trial.convection-old.convection;
    d.diffusion=trial.diffusion-old.diffusion;
    d.production=trial.physical_production-old.physical_production;
    d.destruction=trial.physical_destruction-old.physical_destruction;
    d.cross=trial.physical_cross-old.physical_cross;
    d.sum=d.convection+d.diffusion+d.production+d.destruction+d.cross;
    d.residual=trial.full_residual-old.full_residual;
    d.closure=d.residual-d.sum;
    d.constraint=old.constraint;
    return d;
}
}
