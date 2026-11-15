@@ -0,0 +1,270 @@
#pragma once
// V70: pressure/face-flux conservation for a private six-field SST candidate.
// No new physical residual, turbulence model, Schur approximation or tolerance.
// The supplied face-response must be the SAME Rhie-Chow operator, with U,
// viscosity and candidate-owned rAU frozen during the pressure projection.
// The complete F_phys is re-evaluated afterwards; stale pressure operators
// are never allowed to certify the updated candidate or publish a carrier.
#include "cfdx/physics/comparable_physical_merit_v65.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct BfsMassBalanceV70 {
    double divergence_linf=0.;
    double divergence_l1=0.;
    double boundary_signed=0.;
    double flux_l1=0.;
    double identity_error=0.;
};

template<class Owner,class Neighbour>
BfsMassBalanceV70 bfs_mass_balance_v70(
    std::size_t cells,const std::vector<double>& phi,
    Owner&& owner,Neighbour&& neighbour) {
    if (!cells || phi.empty()) throw std::invalid_argument("V70 empty mass balance");
    std::vector<long double> div(cells,0.);
    long double boundary=0.,sum_flux=0.;
    for(std::size_t f=0;f<phi.size();++f) {
        const auto o=static_cast<std::size_t>(owner(f));
        const auto nb=neighbour(f);
        if(o>=cells || nb < -1 || (nb>=0 && static_cast<std::size_t>(nb)>=cells) ||
           !std::isfinite(phi[f]))
            throw std::invalid_argument("V70 invalid face or flux");
        const auto flux=static_cast<long double>(phi[f]);
        div[o]+=flux;
        if(nb>=0) div[static_cast<std::size_t>(nb)]-=flux;
        else boundary+=flux;
        sum_flux+=std::abs(flux);
    }
    BfsMassBalanceV70 a;
    long double sum_div=0.;
    for(const auto value:div) {
        a.divergence_linf=std::max(a.divergence_linf,
            std::abs(static_cast<double>(value)));
        a.divergence_l1+=std::abs(static_cast<double>(value));
        sum_div+=value;
    }
    a.boundary_signed=static_cast<double>(boundary);
    a.flux_l1=static_cast<double>(sum_flux);
    a.identity_error=std::abs(static_cast<double>(sum_div-boundary));
    if(!std::isfinite(a.divergence_l1) || !std::isfinite(a.flux_l1))
        throw std::runtime_error("V70 nonfinite mass conservation evidence");
    return a;
}

struct BfsPressureResponseV70 {
    PhysicalCandidateV51 state;
    std::vector<double> face_flux;
};

struct BfsPressureProjectionV70 {
    PhysicalCandidateV51 state;
    std::vector<double> face_flux;
    double divergence_before=0.,divergence_after=0.;
    double projection_fraction=1.;
    bool fully_projected=false;
    double pressure_delta_linf=0.,face_delta_linf=0.;
    double rc_frozen_replay_gap=0.;
    double mass_identity_gap=0.;
    std::size_t pressure_unknowns=0;
};

// Small-pilot dense pressure response: the REAL candidate-owned RC response
// provides all columns, including its nonorthogonal/deferred pressure terms.
// The velocity must be corrected with -V*rAU*grad(delta_p), so a
// physical-velocity Rhie-Chow flux gets the correct pressure response.
// This is intentionally limited to the BFS152 diagnostic. Larger meshes must
// use a scalable native pressure solve; this cannot silently become one.
// Reject on singularity, incompressible boundary imbalance, stale rAU/phi,
// nonfinite coefficients, or incomplete conservation. No regularization.
template<class ResponseAtDeltaPressure,class Owner,class Neighbour>
BfsPressureProjectionV70 bfs_project_pressure_v70(
    const PhysicalCandidateV51& candidate,
    const PhysicalResidualV51& residual,
    ResponseAtDeltaPressure&& response_at_delta_pressure,
    Owner&& owner,Neighbour&& neighbour,
    bool fixed_pressure_boundary,std::size_t gauge_cell,
    double absolute_mass_tolerance=1e-10,
    double projection_fraction=1.) {
    validate_physical_candidate_v51(candidate);
    validate_physical_result_v51(candidate,residual);
    const auto n=candidate.cells;
    if(n>152 || gauge_cell>=n || residual.inverse_momentum_diagonal.size()!=3*n ||
       !(absolute_mass_tolerance>0.) || !std::isfinite(absolute_mass_tolerance) ||
       !(projection_fraction>0.) || projection_fraction>1. ||
       !std::isfinite(projection_fraction))
        throw std::invalid_argument("V70 pilot pressure projection domain/coefficients invalid");
    const auto nf=residual.face_flux.size();
    const std::vector<double> zeros(n,0.);
    const auto original=response_at_delta_pressure(zeros);
    if(original.state.cells!=n || original.state.fields!=candidate.fields ||
       original.face_flux.size()!=nf)
        throw std::runtime_error("V70 zero pressure increment changed the candidate");
    const auto& initial_flux=original.face_flux;
    double flux_reference=0.,replay_gap=0.;
    for(std::size_t f=0;f<nf;++f) {
        flux_reference+=std::abs(residual.face_flux[f]);
        replay_gap=std::max(replay_gap,
            std::abs(initial_flux[f]-residual.face_flux[f]));
    }
    // The supplied pressure operator MUST reproduce the same candidate's
    // candidate-owned F_phys face flux. A stale carrier is forbidden.
    if(replay_gap>1e-10*std::max(1.,flux_reference))
        throw std::runtime_error("V70 frozen-pressure Rhie-Chow response does not reproduce F_phys");
    const auto old=bfs_mass_balance_v70(n,residual.face_flux,owner,neighbour);
    const double scale=std::max(1.,flux_reference);
    if(old.identity_error>1e-11*scale)
        throw std::runtime_error("V70 original owner/neighbour divergence inconsistent");
    // Build the exact coupled U/p/phi pressure response from zero and unit
    // pressure increments. The pressure-gradient velocity correction is part
    // of the callback, with U/rAU/viscosity frozen at the same candidate.
    const auto& zero_flux=initial_flux;
    std::vector<std::vector<long double>> matrix(n,std::vector<long double>(n+1,0.));
    // RHS = -div(phi(p))
    for(std::size_t f=0;f<nf;++f) {
        const auto o=static_cast<std::size_t>(owner(f));
        const auto nb=neighbour(f);
        matrix[o][n]-=static_cast<long double>(residual.face_flux[f]);
        if(nb>=0)matrix[static_cast<std::size_t>(nb)][n]+=
            static_cast<long double>(residual.face_flux[f]);
    }
    for(std::size_t j=0;j<n;++j) {
        std::vector<double> unit(n,0.);unit[j]=1.;
        const auto response=response_at_delta_pressure(unit);
        if(response.face_flux.size()!=nf || response.state.cells!=n ||
           response.state.fields.size()!=6*n)
            throw std::runtime_error("V70 unit-pressure response incomplete");
        for(std::size_t f=0;f<nf;++f) {
            const auto o=static_cast<std::size_t>(owner(f));
            const auto nb=neighbour(f);
            const auto value=static_cast<long double>(response.face_flux[f])-
                             static_cast<long double>(zero_flux[f]);
            if(!std::isfinite(static_cast<double>(value)))
                throw std::runtime_error("V70 nonfinite pressure response");
            matrix[o][j]+=value;
            if(nb>=0)matrix[static_cast<std::size_t>(nb)][j]-=value;
        }
    }
    if(!fixed_pressure_boundary) {
        // Pure Neumann: pressure is determined only up to a constant. Do not
        // adjust the physical boundary flux or hide incompatibility by pinning.
        if(std::abs(old.boundary_signed)>absolute_mass_tolerance)
            throw std::runtime_error("V70 incompatible pure-Neumann mass inflow");
        std::fill(matrix[gauge_cell].begin(),matrix[gauge_cell].end(),0.);
        matrix[gauge_cell][gauge_cell]=1.;
    }
    long double max_coeff=0.;
    for(std::size_t i=0;i<n;++i)
        for(std::size_t j=0;j<n;++j)
            max_coeff=std::max(max_coeff,std::abs(matrix[i][j]));
    if(!(max_coeff>0.) || !std::isfinite(static_cast<double>(max_coeff)))
        throw std::runtime_error("V70 zero or nonfinite pressure operator");
    for(std::size_t k=0;k<n;++k) {
        std::size_t pivot=k;
        for(std::size_t i=k+1;i<n;++i)
            if(std::abs(matrix[i][k])>std::abs(matrix[pivot][k]))pivot=i;
        if(std::abs(matrix[pivot][k])<=1e-13L*max_coeff)
            throw std::runtime_error("V70 singular/ill-conditioned pressure projection");
        if(pivot!=k)std::swap(matrix[pivot],matrix[k]);
        for(std::size_t i=k+1;i<n;++i) {
            const long double t=matrix[i][k]/matrix[k][k];
            matrix[i][k]=0.;
            for(std::size_t j=k+1;j<=n;++j)matrix[i][j]-=t*matrix[k][j];
        }
    }
    std::vector<long double> correction(n);
    for(std::size_t ii=n;ii-->0;) {
        long double value=matrix[ii][n];
        for(std::size_t j=ii+1;j<n;++j)value-=matrix[ii][j]*correction[j];
        correction[ii]=value/matrix[ii][ii];
    }
    BfsPressureProjectionV70 out;
    out.divergence_before=old.divergence_linf;
    out.pressure_unknowns=n;
    out.rc_frozen_replay_gap=replay_gap;
    out.projection_fraction=projection_fraction;
    out.fully_projected=projection_fraction==1.;
    std::vector<double> delta(n);
    for(std::size_t i=0;i<n;++i) {
        delta[i]=projection_fraction*static_cast<double>(correction[i]);
        if(!std::isfinite(delta[i]))
            throw std::runtime_error("V70 nonfinite pressure correction");
        out.pressure_delta_linf=std::max(out.pressure_delta_linf,std::abs(delta[i]));
    }
    const auto corrected=response_at_delta_pressure(delta);
    const auto corrected_replay=response_at_delta_pressure(delta);
    if(corrected.state.fields!=corrected_replay.state.fields ||
       corrected.face_flux!=corrected_replay.face_flux ||
       corrected.face_flux.size()!=nf || corrected.state.cells!=n)
        throw std::runtime_error("V70 correction is history-dependent");
    out.state=corrected.state;
    out.face_flux=corrected.face_flux;
    validate_physical_candidate_v51(out.state);
    for(std::size_t i=0;i<n;++i) {
        if(std::abs((out.state.fields[3*n+i]-candidate.fields[3*n+i])-delta[i])>
           1e-10*std::max(1.,std::abs(out.state.fields[3*n+i])))
            throw std::runtime_error("V70 pressure response used the wrong increment");
    }
    for(std::size_t i=4*n;i<6*n;++i)
        if(out.state.fields[i]!=candidate.fields[i])
            throw std::runtime_error("V70 pressure correction modified the SST state");
    for(std::size_t f=0;f<nf;++f)
        out.face_delta_linf=std::max(out.face_delta_linf,
             std::abs(out.face_flux[f]-initial_flux[f]));
    const auto after=bfs_mass_balance_v70(n,out.face_flux,owner,neighbour);
    out.divergence_after=after.divergence_linf;
    out.mass_identity_gap=after.identity_error;
    // A damped fraction is *not* a conservative final projection. Only the
    // fully projected branch may use that term. Every partial step must reduce
    // the frozen mass defect, and none may be published as a closed flux.
    if(after.identity_error>absolute_mass_tolerance ||
       (out.fully_projected && after.divergence_linf>absolute_mass_tolerance) ||
       (!out.fully_projected && !(after.divergence_linf<old.divergence_linf)))
        throw std::runtime_error("V70 pressure/flux step did not meet its mass contract");
    return out;
}

// Only an independently *reassembled* F_phys and the unchanged V65 guards
// can qualify the projected candidate. This is a private staging gate, not
// a production commit; physical coefficients changing with pressure/SST must
// not be concealed behind a frozen operator's zero divergence.
template<class Owner,class Neighbour>
ComparablePhysicalMeritV65 bfs_certify_projected_trial_v70(
    const PhysicalCandidateV51& old_state,
    const PhysicalResidualV51& old_residual,
    const BfsPressureProjectionV70& projected,
    const PhysicalResidualV51& reassembled,
    Owner&& owner,Neighbour&& neighbour,
    double relative_tolerance,double continuity_tolerance,
    double absolute_mass_tolerance=1e-10) {
    validate_physical_result_v51(projected.state,reassembled);
    const auto n=projected.state.cells;
    const auto actual=bfs_mass_balance_v70(n,reassembled.face_flux,owner,neighbour);
    if(actual.divergence_linf>absolute_mass_tolerance ||
       actual.identity_error>absolute_mass_tolerance)
        throw std::runtime_error("V70 reassembled F_phys has non-conservative mass flux");
    for(std::size_t i=0;i<n;++i) {
        double div=0.;
        // Signed FV continuity rows are checked below directly from faces.
        for(std::size_t f=0;f<reassembled.face_flux.size();++f) {
            if(static_cast<std::size_t>(owner(f))==i)div+=reassembled.face_flux[f];
            if(neighbour(f)==static_cast<int>(i))div-=reassembled.face_flux[f];
        }
        if(std::abs(div-reassembled.signed_raw[3*n+i])>absolute_mass_tolerance)
            throw std::runtime_error("V70 F_phys continuity differs from face divergence");
    }
    const auto gate=evaluate_comparable_merit_v65(
        old_state,projected.state,old_residual,reassembled,
        relative_tolerance,continuity_tolerance,false,false);
    if(!gate.accepted)
        throw std::runtime_error("V70 updated candidate fails unchanged V65 merit/growth guards");
    return gate;
}
} // namespace cfdx::physics
