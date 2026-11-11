#pragma once

#include "cfdx/physics/sst_trust_subspace_v27.h"
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

// V28 EXPERIMENTAL: trust-region policy for a *proposed* SST step. The
// projected FV residual evaluator is the final authority. This is NOT PETSc
// SNES and does NOT change the root of the SST equations or its final gates.
// The radius is an inner-iteration state, not a global or static variable.
struct SSTV28AdaptiveTrustState {
    double radius = 0.5;
    static constexpr double minimum_radius = 1.0/64.0;
    static constexpr double maximum_radius = 1.5;
    std::size_t committed = 0;
    std::size_t rejected = 0;
    double last_committed_ratio = std::numeric_limits<double>::quiet_NaN();

    void validate() const {
        if (!std::isfinite(radius) || radius<minimum_radius ||
            radius>maximum_radius)
            throw std::invalid_argument("V28 SST trust: radius outside valid range");
    }
    // Call ONLY if this proposal actually won the final FV trial comparison
    // with alpha=1, and was subsequently installed as the accepted state.
    // Never expand the trust region for a hypothetical uncommitted candidate.
    void commit(double actual_over_predicted,double actual_log_step) {
        validate();
        if (!std::isfinite(actual_over_predicted) ||
            !(actual_over_predicted>1e-3) ||
            !std::isfinite(actual_log_step) || !(actual_log_step>0.) ||
            actual_log_step>radius*(1.+1e-8))
            throw std::invalid_argument("V28 SST trust: invalid accepted step");
        ++committed;
        last_committed_ratio=actual_over_predicted;
        if (actual_over_predicted<0.25)
            radius=std::max(minimum_radius,0.5*radius);
        else if (actual_over_predicted>=0.75 &&
                 actual_log_step>=0.95*radius)
            radius=std::min(maximum_radius,2.*radius);
    }
    // No acceptable V28 proposal: reduce only the MODEL radius. Other solver
    // proposals (e.g. V22 Jacobi) remain available and the CFD state is not
    // changed here.
    void reject() {
        validate();
        ++rejected;
        radius=std::max(minimum_radius,0.5*radius);
    }
};

// Rebuild the expensive spatial Jacobian ONLY after an unsuccessful FV
// globalization step, and at most twice before the original stall protection.
// This makes the rescue path reachable at iteration 22 after a rejection at 21
// without imposing an expensive matrix refresh on every accepted iteration.
inline bool sst_v28_stagnation_probe(bool enabled,std::size_t iteration,
                                     std::size_t consecutive_stalls) noexcept {
    return enabled && iteration>1 && (iteration%10)!=0 &&
           consecutive_stalls>0 && consecutive_stalls<=2;
}

struct SSTV28BlendSearch {
    bool accepted = false;
    std::vector<double> fields;
    double merit = std::numeric_limits<double>::infinity();
    double smooth_l2 = std::numeric_limits<double>::infinity();
    double weight_jacobi = 0.;
    double used_log_radius = 0.;
    double actual_over_predicted = 0.;
    double predicted_reduction = 0.;
    double actual_reduction = 0.;
    std::size_t evaluations = 0;
    std::size_t admissible_trials = 0;
    std::size_t rejected_linf = 0;
    std::size_t rejected_l2 = 0;
    std::size_t rejected_projection = 0;
    std::string reason;
};

// Evaluate logarithmic combinations of the old, separately-certified V22
// Newton-Jacobi and V27 Cauchy candidates. A Cauchy candidate is required;
// isolated Newton-Jacobi is already evaluated by the historical safeguard.
// A ratio rho=(Phi(x)-Phi(x_trial))/[m(0)-m(step)] uses a first-order model
// of the ACTUAL dynamically scaled projected L2 merit, differentiated by FV
// reassembly along each mixture. Neither the model nor L2 is a final gate.
inline SSTV28BlendSearch sst_blend_adaptive_trust_v28(
    const std::vector<double>& initial,
    const std::vector<double>& jacobi_fields,
    const std::vector<double>& cauchy_fields,
    std::size_t n_cells, double k_floor, double omega_floor,
    double initial_linf, double initial_l2, double radius,
    const SSTV25PhysicalEvaluator& evaluate)
{
    if (!n_cells || initial.size()!=2*n_cells ||
        jacobi_fields.size()!=initial.size() ||
        cauchy_fields.size()!=initial.size() || !evaluate ||
        !(k_floor>0.) || !(omega_floor>0.) ||
        !(initial_linf>0.) || !(initial_l2>0.) ||
        !std::isfinite(initial_linf) || !std::isfinite(initial_l2) ||
        !std::isfinite(radius) || radius<SSTV28AdaptiveTrustState::minimum_radius ||
        radius>SSTV28AdaptiveTrustState::maximum_radius)
        throw std::invalid_argument("V28 SST trust: invalid input contract");

    SSTV28BlendSearch out;
    const std::size_t n=initial.size();
    std::vector<double> jacobi_log(n),cauchy_log(n);
    for (std::size_t i=0;i<n;++i) {
        const double floor=(i<n_cells?k_floor:omega_floor);
        if (!std::isfinite(initial[i]) || !(initial[i]>floor) ||
            !std::isfinite(jacobi_fields[i]) || !(jacobi_fields[i]>=floor) ||
            !std::isfinite(cauchy_fields[i]) || !(cauchy_fields[i]>=floor))
            throw std::invalid_argument("V28 SST trust: invalid physical candidate");
        jacobi_log[i]=std::log(jacobi_fields[i]/initial[i]);
        cauchy_log[i]=std::log(cauchy_fields[i]/initial[i]);
        if (!std::isfinite(jacobi_log[i]) || !std::isfinite(cauchy_log[i]))
            throw std::invalid_argument("V28 SST trust: nonfinite log step");
    }

    const auto measured=[&](const std::vector<double>& candidate) {
        auto result=evaluate(candidate);
        ++out.evaluations;
        if (result.fields.size()!=n || !std::isfinite(result.merit) ||
            !(result.merit>=0.) || !std::isfinite(result.smooth_l2) ||
            !(result.smooth_l2>=0.) ||
            !std::isfinite(result.k_local) || !(result.k_local>=0.) ||
            !std::isfinite(result.omega_local) || !(result.omega_local>=0.) ||
            result.k_worst>=n_cells || result.omega_worst>=n_cells ||
            std::abs(result.merit-std::max(result.k_local,result.omega_local))>
                1e-12*std::max(1.,result.merit))
            throw std::domain_error("V28 SST trust: invalid FV residual certificate");
        for (std::size_t i=0;i<n;++i) {
            const double floor=i<n_cells?k_floor:omega_floor;
            if (!std::isfinite(result.fields[i]) || !(result.fields[i]>=floor))
                throw std::domain_error("V28 SST trust: invalid wall-projected field");
        }
        return result;
    };
    try {
        const auto base=measured(initial);
        if (std::abs(base.merit-initial_linf)>
              1e-7*std::max(1.,initial_linf) ||
            std::abs(base.smooth_l2-initial_l2)>
              1e-7*std::max(1.,initial_l2)) {
            out.reason="baseline differs from current FV certificate";
            return out;
        }
    } catch (const std::exception& e) {
        out.reason=std::string("invalid projected baseline: ")+e.what();
        return out;
    }
    const double phi0=0.5*initial_l2*initial_l2;
    for (double weight:{0.25,0.50,0.75}) {
        std::vector<double> log_direction(n);
        double norm_inf=0.;
        for (std::size_t i=0;i<n;++i) {
            log_direction[i]=(1.-weight)*cauchy_log[i]+weight*jacobi_log[i];
            norm_inf=std::max(norm_inf,std::abs(log_direction[i]));
        }
        if (!(norm_inf>1e-14) || !std::isfinite(norm_inf)) continue;
        const auto exponential_candidate=[&](double log_amplitude) {
            auto candidate=initial;
            for (std::size_t i=0;i<n;++i) {
                const double floor=i<n_cells?k_floor:omega_floor;
                candidate[i]*=std::exp(log_amplitude*log_direction[i]/norm_inf);
                if (!std::isfinite(candidate[i]) || !(candidate[i]>floor))
                    throw std::domain_error("V28 SST trust: invalid exponential trial");
            }
            return candidate;
        };
        // A measured directional derivative is preferable to assuming Jv of
        // the frozen-scale residual equals a derivative of the dynamic merit.
        double slope=0.;
        try {
            const double probe=std::min(1e-3,0.1*radius);
            const auto p=measured(exponential_candidate(probe));
            slope=(0.5*p.smooth_l2*p.smooth_l2-phi0)/probe;
        } catch (const std::exception&) {continue;}
        if (!std::isfinite(slope) || !(slope<-1e-15)) continue;
        // Small bounded search; two trial radii per mixture. The original
        // Jacobi and V27 candidates continue to compete independently.
        const double maximum_step=std::min(radius,norm_inf);
        for (double factor:{1.,0.5}) {
            const double step=maximum_step*factor;
            if (!(step>0.)) continue;
            try {
                const auto trial=measured(exponential_candidate(step));
                double projected_size=0.;
                for (std::size_t i=0;i<n;++i)
                    projected_size=std::max(projected_size,
                        std::abs(std::log(trial.fields[i]/initial[i])));
                if (!std::isfinite(projected_size) ||
                    projected_size>radius*(1.+1e-9)) {
                    ++out.rejected_projection;continue;
                }
                if (!(trial.merit<initial_linf*(1.-1e-9*step))) {
                    ++out.rejected_linf;continue;
                }
                if (!(trial.smooth_l2<initial_l2*(1.-1e-9*step))) {
                    ++out.rejected_l2;continue;
                }
                const double predicted=-slope*step;
                const double actual=phi0-0.5*trial.smooth_l2*trial.smooth_l2;
                const double rho=actual/predicted;
                if (!std::isfinite(rho) || !(predicted>0.) ||
                    !(actual>0.) || !(rho>1e-3)) continue;
                ++out.admissible_trials;
                if (!out.accepted || trial.merit<out.merit ||
                    (trial.merit==out.merit && trial.smooth_l2<out.smooth_l2)) {
                    out.accepted=true;
                    out.fields=trial.fields;
                    out.merit=trial.merit;
                    out.smooth_l2=trial.smooth_l2;
                    out.weight_jacobi=weight;
                    out.used_log_radius=projected_size;
                    out.actual_over_predicted=rho;
                    out.predicted_reduction=predicted;
                    out.actual_reduction=actual;
                }
            } catch (const std::exception&) {
                // Trial-only failure is not allowed to mutate caller fields.
            }
        }
    }
    out.reason=out.accepted?
        "blend reduced both real FV merits with positive trust ratio":
        "no log blend passed projected Linf, L2 and model-ratio gates";
    return out;
}

} // namespace cfdx::physics
