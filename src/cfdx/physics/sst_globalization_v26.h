#pragma once

#include "cfdx/physics/sst_globalization_v25.h"
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

// V26 diagnostic of the SAME unprojected nonlinear residual used for the
// Newton-Krylov Jv operator. Here d is normalized in log(k),log(omega), so
// h is a dimensionless LOG perturbation independent of a large Krylov step.
// Wall projection is deliberately NOT applied inside this measurement:
// owner-cell constraints are equations in the Newton operator and must not
// be eliminated by projection before differentiation.
struct SSTV26DirectionalAudit {
    bool available=false;
    std::array<double,3> steps{{1e-3,1e-4,1e-5}};
    std::array<double,3> derivative_norms{{0.,0.,0.}};
    std::array<double,3> one_sided_disagreement{{0.,0.,0.}};
    std::array<double,2> adjacent_scale_disagreement{{0.,0.}};
    double maximum_disagreement=std::numeric_limits<double>::infinity();
    std::size_t evaluations=0;
    std::string reason;
};

inline SSTV26DirectionalAudit sst_audit_log_jvp_v26(
    const std::vector<double>& physical,
    const std::vector<double>& log_direction,
    std::size_t n_cells,double k_floor,double omega_floor,
    const SSTPhysicalResidualV22& physical_residual,
    std::array<double,3> steps={{1e-3,1e-4,1e-5}})
{
    if (!n_cells || physical.size()!=2*n_cells ||
        log_direction.size()!=physical.size() || !physical_residual ||
        !(k_floor>0.) || !(omega_floor>0.) ||
        !(steps[0]>steps[1] && steps[1]>steps[2] && steps[2]>0.) ||
        steps[0]>0.01)
        throw std::invalid_argument("V26 log Jv audit: invalid contract");
    SSTV26DirectionalAudit out;
    out.steps=steps;
    double max_direction=0.;
    for (std::size_t i=0;i<physical.size();++i) {
        const double floor=i<n_cells?k_floor:omega_floor;
        if (!std::isfinite(physical[i]) || !(physical[i]>floor) ||
            !std::isfinite(log_direction[i]))
            throw std::invalid_argument("V26 log Jv audit: invalid state");
        max_direction=std::max(max_direction,std::abs(log_direction[i]));
    }
    if (!(max_direction>0.)) {
        out.reason="zero log direction";
        return out;
    }
    // Consistent central and one-sided differences expose cancellation,
    // source-limiter nonsmoothness and changing wall-function branches.
    try {
        auto base=physical_residual(physical);
        ++out.evaluations;
        if (base.size()!=physical.size())
            throw std::domain_error("incomplete residual");
        for (double x:base) if (!std::isfinite(x))
            throw std::domain_error("non-finite base residual");
        std::array<std::vector<double>,3> derivatives;
        for (std::size_t s=0;s<3;++s) {
            auto plus=physical,minus=physical;
            for (std::size_t i=0;i<physical.size();++i) {
                const double delta=steps[s]*log_direction[i]/max_direction;
                plus[i]*=std::exp(delta);
                minus[i]*=std::exp(-delta);
                const double floor=i<n_cells?k_floor:omega_floor;
                if (!std::isfinite(plus[i]) || !std::isfinite(minus[i]) ||
                    !(plus[i]>floor) || !(minus[i]>floor))
                    throw std::domain_error("perturbation violates positivity");
            }
            const auto fp=physical_residual(plus);
            const auto fm=physical_residual(minus);
            out.evaluations+=2;
            if (fp.size()!=physical.size() || fm.size()!=physical.size())
                throw std::domain_error("incomplete trial residual");
            long double sq=0.,one_sided_sq=0.;
            derivatives[s].resize(base.size());
            for (std::size_t i=0;i<base.size();++i) {
                if (!std::isfinite(fp[i]) || !std::isfinite(fm[i]))
                    throw std::domain_error("non-finite trial residual");
                const double central=(fp[i]-fm[i])/(2.*steps[s]);
                const double forward=(fp[i]-base[i])/steps[s];
                const double backward=(base[i]-fm[i])/steps[s];
                if (!std::isfinite(central) || !std::isfinite(forward) ||
                    !std::isfinite(backward))
                    throw std::domain_error("difference quotient overflow");
                derivatives[s][i]=central;
                sq+=static_cast<long double>(central)*central;
                const auto mismatch=static_cast<long double>(forward)-backward;
                one_sided_sq+=mismatch*mismatch;
            }
            out.derivative_norms[s]=static_cast<double>(std::sqrt(sq));
            const auto denominator=std::max(out.derivative_norms[s],1e-300);
            out.one_sided_disagreement[s]=
                static_cast<double>(std::sqrt(one_sided_sq))/denominator;
        }
        for (std::size_t s=0;s<2;++s) {
            long double sq=0.;
            for (std::size_t i=0;i<base.size();++i) {
                const auto delta=static_cast<long double>(derivatives[s][i])-
                                 derivatives[s+1][i];
                sq+=delta*delta;
            }
            out.adjacent_scale_disagreement[s]=
                static_cast<double>(std::sqrt(sq)) /
                std::max(out.derivative_norms[s+1],1e-300);
        }
        out.maximum_disagreement=0.;
        for (double v:out.one_sided_disagreement)
            out.maximum_disagreement=std::max(out.maximum_disagreement,v);
        for (double v:out.adjacent_scale_disagreement)
            out.maximum_disagreement=std::max(out.maximum_disagreement,v);
        out.available=std::isfinite(out.maximum_disagreement);
        out.reason=out.available?"measured only; not a convergence certificate":
                                   "non-finite discrepancy";
    } catch (const std::exception& e) {
        out.available=false;
        out.reason=e.what();
    }
    return out;
}

// Two independent candidate-selection criteria from ONE set of FV evaluations:
// (a) exact V25 min-Linf among physically decreasing candidates;
// (b) min balanced RMS-L2 among candidates decreasing BOTH L2 and physical Linf.
// The caller must still pass either proposal through the UNCHANGED final
// reassembled Linf SST gate. In particular, an L2-only improvement is never
// an admissible physical step. No PDE, BC or convergence threshold changes.
struct SSTV26SmoothSearch {
    SSTV25GlobalizationResult linf;
    bool smooth_accepted=false;
    std::vector<double> smooth_fields;
    double smooth_alpha=0.;
    double baseline_l2=std::numeric_limits<double>::infinity();
    double best_trial_l2=std::numeric_limits<double>::infinity();
    double selected_l2=std::numeric_limits<double>::infinity();
    double selected_linf=std::numeric_limits<double>::infinity();
    std::size_t smooth_improving=0;
};

inline SSTV26SmoothSearch sst_globalize_log_direction_v26(
    const std::vector<double>& initial,
    const std::vector<double>& log_direction,
    std::size_t n_cells,double k_floor,double omega_floor,
    double baseline_linf,double baseline_l2,
    const SSTV25PhysicalEvaluator& evaluate,
    std::size_t maximum_backtracks=12,double log_trust_radius=1.5)
{
    if (!n_cells || initial.size()!=2*n_cells ||
        log_direction.size()!=initial.size() || !evaluate ||
        !(k_floor>0.) || !(omega_floor>0.) ||
        !(baseline_linf>0.) || !std::isfinite(baseline_linf) ||
        !(baseline_l2>0.) || !std::isfinite(baseline_l2) ||
        !(log_trust_radius>0.) || !std::isfinite(log_trust_radius) ||
        maximum_backtracks>24)
        throw std::invalid_argument("V26 smooth log globalizer: invalid contract");
    SSTV26SmoothSearch out;
    out.baseline_l2=baseline_l2;
    for (std::size_t i=0;i<initial.size();++i) {
        const double floor=i<n_cells?k_floor:omega_floor;
        if (!(initial[i]>floor) || !std::isfinite(initial[i]) ||
            !std::isfinite(log_direction[i]))
            throw std::invalid_argument("V26 smooth log globalizer: invalid state");
        out.linf.maximum_absolute_log_direction=std::max(
            out.linf.maximum_absolute_log_direction,std::abs(log_direction[i]));
    }
    if (!(out.linf.maximum_absolute_log_direction>0.)) {
        out.linf.reason="zero Newton log direction";
        return out;
    }
    const double start_alpha=std::min(1.,log_trust_radius/
                                          out.linf.maximum_absolute_log_direction);
    for (std::size_t step=0;step<=maximum_backtracks;++step) {
        const double alpha=std::ldexp(start_alpha,-static_cast<int>(step));
        ++out.linf.attempted;
        auto candidate=initial;
        bool feasible=true;
        for (std::size_t i=0;i<candidate.size();++i) {
            const auto floor=i<n_cells?k_floor:omega_floor;
            candidate[i]*=std::exp(alpha*log_direction[i]);
            if (!std::isfinite(candidate[i]) || !(candidate[i]>=floor)) {
                feasible=false;break;
            }
        }
        if (!feasible) continue;
        try {
            auto measured=evaluate(candidate);
            if (measured.fields.size()!=initial.size() ||
                !(measured.merit>=0.) || !std::isfinite(measured.merit) ||
                !(measured.smooth_l2>=0.) || !std::isfinite(measured.smooth_l2) ||
                !(measured.k_local>=0.) || !std::isfinite(measured.k_local) ||
                !(measured.omega_local>=0.) || !std::isfinite(measured.omega_local) ||
                std::abs(measured.merit-
                   std::max(measured.k_local,measured.omega_local))>
                   1e-12*std::max(1.,measured.merit) ||
                measured.k_worst>=n_cells || measured.omega_worst>=n_cells)
                throw std::domain_error("incomplete projected physical certificate");
            for (std::size_t i=0;i<measured.fields.size();++i) {
                const auto floor=i<n_cells?k_floor:omega_floor;
                if (!std::isfinite(measured.fields[i]) || measured.fields[i]<floor)
                    throw std::domain_error("invalid projected state");
            }
            ++out.linf.valid;
            if (step==0) out.linf.full_step_merit=measured.merit;
            out.linf.best_trial_merit=std::min(out.linf.best_trial_merit,
                                               measured.merit);
            out.best_trial_l2=std::min(out.best_trial_l2,measured.smooth_l2);
            const bool physical_ok=measured.merit<baseline_linf*(1.-1e-9*alpha);
            if (!physical_ok) continue;
            ++out.linf.improving;
            if (!out.linf.accepted || measured.merit<out.linf.merit) {
                out.linf.accepted=true;
                out.linf.alpha=alpha;
                out.linf.merit=measured.merit;
                out.linf.fields=measured.fields;
                out.linf.k_worst=measured.k_worst;
                out.linf.omega_worst=measured.omega_worst;
            }
            if (measured.smooth_l2 < baseline_l2*(1.-1e-9*alpha)) {
                ++out.smooth_improving;
                if (!out.smooth_accepted || measured.smooth_l2<out.selected_l2) {
                    out.smooth_accepted=true;
                    out.smooth_fields=std::move(measured.fields);
                    out.smooth_alpha=alpha;
                    out.selected_linf=measured.merit;
                    out.selected_l2=measured.smooth_l2;
                }
            }
        } catch (const std::exception&) {
            // Failed wall refresh / nonfinite trial cannot promote a state.
        }
    }
    out.linf.reason=out.linf.accepted?"strict physical merit decreased":
                        "no projected physical step passed strict decrease";
    return out;
}

} // namespace cfdx::physics
