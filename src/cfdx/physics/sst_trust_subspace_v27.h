#pragma once

#include "cfdx/physics/sst_globalization_v26.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

// Experimental V27: sampled *projected physical* Cauchy directions in a
// matrix-free 2--4 dimensional subspace, not a PETSc trust-region solver.
// Basis: spatial JFNK, admissible Picard log increment, and one active cell
// from each FV equation. The smooth objective is Phi=0.5*M2^2, with dynamic
// row normalization from the production FV assembler. All finite differences
// pass through the same wall projection, nu_t/G refresh and physical evaluator
// as the FINAL Linf gate. Neither a linear model nor M2 alone accepts a step.
struct SSTV27TrustSearch {
    bool accepted=false;
    std::vector<double> fields;
    double merit=std::numeric_limits<double>::infinity();
    double smooth_l2=std::numeric_limits<double>::infinity();
    double initial_linf=std::numeric_limits<double>::infinity();
    double initial_l2=std::numeric_limits<double>::infinity();
    double chosen_radius=0.;
    double ratio=0.;
    double directional_l2_slope=0.;
    std::size_t bases=0;
    std::size_t evaluations=0;
    std::size_t valid_trials=0;
    std::size_t rejected_linf=0;
    std::size_t rejected_l2=0;
    std::size_t rejected_trust=0;
    std::size_t k_active=0,omega_active=0;
    std::string direction="none";
    std::string reason;
};

inline SSTV27TrustSearch sst_cauchy_subspace_v27(
    const std::vector<double>& initial,
    const std::vector<double>& spatial_log_direction,
    const std::vector<double>& picard_state,
    std::size_t n_cells,double k_floor,double omega_floor,
    double initial_linf,double initial_l2,
    std::size_t k_active,std::size_t omega_active,
    const SSTV25PhysicalEvaluator& evaluate,
    double initial_radius=0.5,std::size_t maximum_halvings=7,
    double probe=1e-3)
{
    if (!n_cells || initial.size()!=2*n_cells ||
        spatial_log_direction.size()!=initial.size() ||
        (!picard_state.empty() && picard_state.size()!=initial.size()) ||
        !evaluate || k_active>=n_cells || omega_active>=n_cells ||
        !(k_floor>0.) || !(omega_floor>0.) ||
        !(initial_linf>0.) || !(initial_l2>0.) ||
        !std::isfinite(initial_linf) || !std::isfinite(initial_l2) ||
        !(initial_radius>0.) || initial_radius>1.5 ||
        maximum_halvings>12 || !(probe>0.) || probe>1e-2)
        throw std::invalid_argument("V27 SST subspace: invalid contract");
    SSTV27TrustSearch result;
    result.initial_linf=initial_linf;
    result.initial_l2=initial_l2;
    result.k_active=k_active;
    result.omega_active=omega_active;
    for (std::size_t i=0;i<initial.size();++i) {
        const double floor=i<n_cells?k_floor:omega_floor;
        if (!std::isfinite(initial[i]) || !(initial[i]>floor) ||
            !std::isfinite(spatial_log_direction[i]))
            throw std::invalid_argument("V27 SST subspace: non-interior starting state");
    }
    // Reject a callback whose baseline differs from the current FV certificate.
    // No caller-owned state is ever changed on failure.
    auto checked=[&](const std::vector<double>& candidate) {
        auto trial=evaluate(candidate);
        ++result.evaluations;
        if (trial.fields.size()!=initial.size() ||
            !std::isfinite(trial.merit) || trial.merit<0. ||
            !std::isfinite(trial.smooth_l2) || trial.smooth_l2<0. ||
            !std::isfinite(trial.k_local) || trial.k_local<0. ||
            !std::isfinite(trial.omega_local) || trial.omega_local<0. ||
            trial.k_worst>=n_cells || trial.omega_worst>=n_cells ||
            std::abs(trial.merit-std::max(trial.k_local,trial.omega_local))>
                1e-12*std::max(1.,trial.merit))
            throw std::domain_error("V27 SST subspace: invalid FV certificate");
        for (std::size_t i=0;i<trial.fields.size();++i) {
            const double floor=i<n_cells?k_floor:omega_floor;
            if (!std::isfinite(trial.fields[i]) || !(trial.fields[i]>=floor))
                throw std::domain_error("V27 SST subspace: invalid projected state");
        }
        return trial;
    };
    try {
        const auto base=checked(initial);
        const double tol=1e-7;
        if (std::abs(base.merit-initial_linf)>tol*std::max(1.,initial_linf) ||
            std::abs(base.smooth_l2-initial_l2)>tol*std::max(1.,initial_l2)) {
            result.reason="projected baseline differs from original FV certificate";
            return result;
        }
    } catch (const std::exception& e) {
        result.reason=std::string("invalid baseline: ")+e.what();
        return result;
    }
    const std::size_t size=initial.size();
    std::vector<std::vector<double>> basis;
    auto append_normalized=[&](std::vector<double> direction) {
        double norm=0.;
        for (double d:direction) {
            if (!std::isfinite(d)) return;
            norm=std::max(norm,std::abs(d));
        }
        if (!(norm>1e-14)) return;
        for (double& d:direction) d/=norm;
        basis.push_back(std::move(direction));
    };
    append_normalized(spatial_log_direction);
    if (!picard_state.empty()) {
        std::vector<double> pd(size);
        bool valid=true;
        for (std::size_t i=0;i<size;++i) {
            const double floor=i<n_cells?k_floor:omega_floor;
            if (!std::isfinite(picard_state[i]) || picard_state[i]<floor) {
                valid=false;break;
            }
            pd[i]=std::log(picard_state[i]/initial[i]);
            if (!std::isfinite(pd[i])) {valid=false;break;}
        }
        if (valid) append_normalized(std::move(pd));
    }
    // The active coordinates are selected by the largest *measured* FV row
    // defect (k and omega separately), not by a hard-coded mesh cell id.
    auto ek=std::vector<double>(size,0.);
    ek[k_active]=1.;
    append_normalized(std::move(ek));
    auto ew=std::vector<double>(size,0.);
    ew[n_cells+omega_active]=1.;
    append_normalized(std::move(ew));
    result.bases=basis.size();

    // Sample genuine PROJECTED nonlinear merits in each subspace direction.
    // This differentiates the wall constraints and dynamic row scaling, unlike
    // the frozen-scale GMRES Jacobian. Near an active floor, skip that basis.
    std::vector<double> grad_l2(basis.size(),0.),grad_linf(basis.size(),0.);
    const auto perturb=[&](const std::vector<double>& b,double h) {
        auto candidate=initial;
        for (std::size_t i=0;i<size;++i) {
            candidate[i]*=std::exp(h*b[i]);
            const double floor=i<n_cells?k_floor:omega_floor;
            if (!std::isfinite(candidate[i]) || !(candidate[i]>floor))
                throw std::domain_error("V27 SST subspace: positivity during probing");
        }
        return candidate;
    };
    for (std::size_t j=0;j<basis.size();++j) {
        try {
            const auto plus=checked(perturb(basis[j],probe));
            const auto minus=checked(perturb(basis[j],-probe));
            grad_l2[j]=(0.5*plus.smooth_l2*plus.smooth_l2-
                         0.5*minus.smooth_l2*minus.smooth_l2)/(2.*probe);
            grad_linf[j]=(plus.merit-minus.merit)/(2.*probe);
            if (!std::isfinite(grad_l2[j]) || !std::isfinite(grad_linf[j]))
                throw std::domain_error("V27 SST subspace: invalid gradient");
        } catch (const std::exception&) {
            grad_l2[j]=0.;grad_linf[j]=0.;
        }
    }
    // The Cauchy direction in this measured subspace is the negative
    // combination of coordinate derivatives. The true directional slope is
    // -sum_j grad_j^2 / (normalization), even for a nonorthogonal basis.
    auto make_direction=[&](const std::vector<double>& gradient) {
        std::vector<double> d(size,0.);
        double largest=0.;
        for (double g:gradient) largest=std::max(largest,std::abs(g));
        if (!(largest>1e-18)) return d;
        for (std::size_t j=0;j<basis.size();++j)
            for (std::size_t i=0;i<size;++i)
                d[i]-=(gradient[j]/largest)*basis[j][i];
        double dnorm=0.;
        for (double v:d) dnorm=std::max(dnorm,std::abs(v));
        if (dnorm>0.) for (double& v:d) v/=dnorm;
        return d;
    };
    auto l2_cauchy=make_direction(grad_l2);
    auto linf_restoration=make_direction(grad_linf);
    std::vector<std::pair<std::string,std::vector<double>>> directions;
    directions.emplace_back("projected_l2_cauchy",std::move(l2_cauchy));
    directions.emplace_back("active_linf_restoration",std::move(linf_restoration));
    const double phi0=0.5*initial_l2*initial_l2;
    for (const auto& entry:directions) {
        const auto& d=entry.second;
        double maximum=0.;
        for (double v:d) maximum=std::max(maximum,std::abs(v));
        if (!(maximum>1e-14)) continue;
        // Differentiate Phi along d using actual projected FV evaluations;
        // protects against basis dependence, source switches and rounding.
        double directional_slope=0.;
        try {
            const auto plus=checked(perturb(d,probe));
            const auto minus=checked(perturb(d,-probe));
            directional_slope=(0.5*plus.smooth_l2*plus.smooth_l2-
                0.5*minus.smooth_l2*minus.smooth_l2)/(2.*probe);
        } catch (const std::exception&) {continue;}
        if (!std::isfinite(directional_slope) || !(directional_slope<-1e-15))
            continue;
        for (std::size_t step=0;step<=maximum_halvings;++step) {
            const double radius=std::ldexp(initial_radius,-static_cast<int>(step));
            try {
                const auto trial=checked(perturb(d,radius));
                ++result.valid_trials;
                // Enforce trust radius on the post-projection physical state,
                // not merely on the trial before the owner-wall projection.
                double projected_log_size=0.;
                for (std::size_t i=0;i<size;++i)
                    projected_log_size=std::max(projected_log_size,
                        std::abs(std::log(trial.fields[i]/initial[i])));
                if (!std::isfinite(projected_log_size) ||
                    projected_log_size>radius*(1.+1e-9)) {
                    ++result.rejected_trust;continue;
                }
                if (!(trial.merit<initial_linf*(1.-1e-9*radius))) {
                    ++result.rejected_linf;continue;
                }
                if (!(trial.smooth_l2<initial_l2*(1.-1e-9*radius))) {
                    ++result.rejected_l2;continue;
                }
                // First-order Cauchy prediction: positive only if locally
                // descending for the actual dynamic-L2 projected evaluator.
                const double predicted=-directional_slope*radius;
                if (!(predicted>0.)) continue;
                const double actual=phi0-0.5*trial.smooth_l2*trial.smooth_l2;
                const double ratio=actual/predicted;
                if (!std::isfinite(ratio) || !(ratio>1e-3)) continue;
                if (!result.accepted || trial.merit<result.merit ||
                    (trial.merit==result.merit && trial.smooth_l2<result.smooth_l2)) {
                    result.accepted=true;
                    result.fields=trial.fields;
                    result.merit=trial.merit;
                    result.smooth_l2=trial.smooth_l2;
                    result.chosen_radius=radius;
                    result.ratio=ratio;
                    result.directional_l2_slope=directional_slope;
                    result.direction=entry.first;
                }
            } catch (const std::exception&) {
                // Trial is rejected; the caller's state remains intact.
            }
        }
    }
    result.reason=result.accepted?
        "strict projected Linf and L2 reduction with positive trust ratio":
        "no Cauchy or active-row step passed strict physical gates";
    return result;
}

} // namespace cfdx::physics
