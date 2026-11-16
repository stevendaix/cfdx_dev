#pragma once

// CFDX V82 -- bounded, positive, fully reassembled SST trial selection.
// The caller builds a frozen [k,omega] block/spatial preconditioner from the
// production FV equations. Its direction is NOT a Newton certificate. Each
// candidate is projected through the CURRENT wall operator and judged by a
// separate production FV reassembly through `evaluator`.
//
// This is a steady nonlinear search: no pseudo-time term enters F(k,omega).
// U, pressure and Rhie-Chow faces are frozen during this INNER search; the
// caller's global U/p/flux/SST transaction remains the final authority.
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

struct SSTV82PhysicalTrial {
    std::vector<double> state;  // post-wall projection, [k0..kn,w0..wn]
    double k_local = std::numeric_limits<double>::infinity();
    double omega_local = std::numeric_limits<double>::infinity();
    double k_row = std::numeric_limits<double>::infinity();
    double omega_row = std::numeric_limits<double>::infinity();
};
using SSTV82Evaluator =
    std::function<SSTV82PhysicalTrial(const std::vector<double>&)>;

struct SSTV82SearchResult {
    bool accepted = false;
    // Preserves legacy Picard as a PRIVATE continuation if a strict FV
    // descent does not exist. This is NOT a nonlinear acceptance certificate.
    bool legacy_fallback_available = false;
    std::vector<double> legacy_fallback_state;
    double legacy_fallback_merit = std::numeric_limits<double>::infinity();
    std::vector<double> state;
    std::size_t evaluations = 0;
    double alpha = 0.0;
    double merit_before = std::numeric_limits<double>::infinity();
    double merit_after = std::numeric_limits<double>::infinity();
    double k_before = std::numeric_limits<double>::infinity();
    double omega_before = std::numeric_limits<double>::infinity();
    double k_after = std::numeric_limits<double>::infinity();
    double omega_after = std::numeric_limits<double>::infinity();
    std::string selected = "none";
    std::string reason = "no_reassembled_merit_descent";
};

inline bool sst_v82_compatible_trial(const SSTV82PhysicalTrial& t,
                                     std::size_t n, double k_floor,
                                     double omega_floor)
{
    if (t.state.size()!=2*n || !(t.k_local>=0.) || !(t.omega_local>=0.) ||
        !(t.k_row>=0.) || !(t.omega_row>=0.) ||
        !std::isfinite(t.k_local) || !std::isfinite(t.omega_local) ||
        !std::isfinite(t.k_row) || !std::isfinite(t.omega_row)) return false;
    for (std::size_t i=0; i<n; ++i)
        if (!std::isfinite(t.state[i]) || !std::isfinite(t.state[n+i]) ||
            t.state[i]<k_floor || t.state[n+i]<omega_floor) return false;
    return true;
}

// The input directions are in LOG variables. A single global scaling preserves
// their off-cell and k/omega couplings (individual clipping would not). Try
// both the implicit spatial block direction and the original Picard direction
// on the SAME accepted-state baseline, choosing only a real FV dual descent.
// Except for `baseline`, no state is implicitly accepted or mutated.
inline SSTV82SearchResult sst_implicit_spatial_search_v82(
    const SSTV82PhysicalTrial& baseline,
    const std::vector<double>& spatial_log_direction,
    const std::vector<double>& picard_log_direction,
    std::size_t n, double k_floor, double omega_floor,
    const SSTV82Evaluator& evaluator,
    const std::vector<double>& raw_picard_state = {},
    double max_log_increment = 0.5)
{
    if (!n || !evaluator || !(k_floor>0.) || !(omega_floor>0.) ||
        !std::isfinite(k_floor) || !std::isfinite(omega_floor) ||
        !(max_log_increment>0.) || !(max_log_increment<=2.) ||
        !std::isfinite(max_log_increment) ||
        !sst_v82_compatible_trial(baseline,n,k_floor,omega_floor))
        throw std::invalid_argument("V82 SST: invalid physical baseline");
    SSTV82SearchResult out;
    out.k_before = baseline.k_local;
    out.omega_before = baseline.omega_local;
    out.merit_before = std::max(baseline.k_local,baseline.omega_local);
    out.merit_after = out.merit_before;
    out.k_after = baseline.k_local;
    out.omega_after = baseline.omega_local;
    // This is a norm globalization, NOT a Pareto requirement on each
    // equation. Initially omega can dominate by 1e3-1e6; forbidding any
    // transient increase of Rk needlessly blocks a massive reduction of
    // max(Rk,Romega). Both component defects are still logged and the FINAL
    // global U,p,continuity,k,omega guards are unchanged in the carrier.
    const auto consider=[&](const SSTV82PhysicalTrial& t,
                            double alpha,const char* name) {
        if (!sst_v82_compatible_trial(t,n,k_floor,omega_floor)) return;
        const double merit=std::max(t.k_local,t.omega_local);
        if (!(merit<out.merit_before*(1.-1e-9*alpha)) ||
            !(merit<out.merit_after)) return;
        out.accepted=true;
        out.state=t.state;
        out.merit_after=merit;
        out.k_after=t.k_local;
        out.omega_after=t.omega_local;
        out.alpha=alpha;
        out.selected=name;
        out.reason="reassembled_global_max_descent";
    };
    // Critical for near-wall stiffness: full physical Picard can change
    // omega by orders of magnitude while remaining positive. Truncating its
    // logarithm to 0.5 (as for approximate-Newton directions) can create a
    // vastly inferior candidate. Evaluate full Picard BEFORE any clipping.
    if (raw_picard_state.size()==2*n) {
        bool admissible=true;
        for (std::size_t i=0;i<2*n;++i) {
            const double floor=i<n?k_floor:omega_floor;
            if (!std::isfinite(raw_picard_state[i]) || raw_picard_state[i]<floor) {
                admissible=false;break;
            }
        }
        if (admissible) {
            try {
                const auto t=evaluator(raw_picard_state);
                ++out.evaluations;
                if (sst_v82_compatible_trial(t,n,k_floor,omega_floor)) {
                    out.legacy_fallback_available=true;
                    out.legacy_fallback_state=t.state;
                    out.legacy_fallback_merit=std::max(t.k_local,t.omega_local);
                }
                consider(t,1.,"picard_full");
            } catch (const std::exception&) {/*candidate rejected*/}
        }
    }
    const auto attempt=[&](const std::vector<double>& direction,const char* name) {
        if (direction.size()!=2*n) return;
        double largest=0.;
        for (double v:direction) {
            if (!std::isfinite(v)) return;
            largest=std::max(largest,std::abs(v));
        }
        if (!(largest>0.)) return;
        const double factor=std::min(1.,max_log_increment/largest);
        for (double alpha: {1.,0.25,0.0625}) {
            auto proposed=baseline.state;
            bool viable=true;
            for (std::size_t i=0;i<2*n;++i) {
                const double z=alpha*factor*direction[i];
                proposed[i]*=std::exp(z);
                const double floor=i<n?k_floor:omega_floor;
                if (!std::isfinite(proposed[i]) || proposed[i]<floor) {
                    viable=false;break;
                }
            }
            if (!viable) continue;
            try {
                const auto t=evaluator(proposed);
                ++out.evaluations;
                consider(t,alpha*factor,name);
            } catch (const std::exception&) {/*candidate rejected*/}
        }
    };
    attempt(spatial_log_direction,"fv_spatial_2x2");
    if (!out.accepted) attempt(picard_log_direction,"picard_log");
    return out;
}

} // namespace cfdx::physics