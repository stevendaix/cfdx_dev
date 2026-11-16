#pragma once
// V85: strictly opt-in *progression* ledger for the same independent FV U,p,phi,k,omega
// residual as V8. Does not certify convergence, does not modify SST equations.
// Larger startup windows must be paid back: each temporary excursion consumes a
// bounded debt, after which recovery is required. All real flux/positivity guards
// remain with V8 and the production transaction.
#include "cfdx/physics/global_nonlinear_residual_v8.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <stdexcept>
#include <string>

namespace cfdx::physics {
struct GlobalOuterTrialV85 {
    bool allowed=false;
    bool baseline=false;
    bool nonmonotone=false;
    std::string reason;
    std::size_t k_debt=0, omega_debt=0;
    double merit_ratio=0., momentum_ratio=0., k_ratio=0., omega_ratio=0.;
};

class GlobalOuterProgressLedgerV85 {
public:
    static constexpr std::size_t window=6;
    // Progress policy ONLY; final convergence values remain those of V8.
    static constexpr double sst_best_multiplier=8.;
    static constexpr double sst_previous_multiplier=3.5;
    GlobalOuterTrialV85 evaluate(const GlobalNonlinearResidualV8& s,
                                const ConvergenceCriteria& tol,
                                bool same_operator=true,
                                double merit_window_growth_cap=1.01) const {
        const auto bad=[](const std::string& reason){
            GlobalOuterTrialV85 r;r.reason=reason;return r;
        };
        if (!s.all_finite || !s.sst || !s.bounded ||
            s.algorithm!=PressureVelocityAlgorithm::SIMPLE ||
            (s.converged && !s.sst_transport_committed) ||
            !finite_nonnegative(s.merit) || !finite_nonnegative(s.continuity) ||
            !finite_nonnegative(s.k) || !finite_nonnegative(s.omega) ||
            !finite_nonnegative(s.eddy_viscosity_closure) ||
            !finite_nonnegative(s.momentum.global_relative) ||
            !finite_nonnegative(s.momentum.max_active_component_relative))
            throw std::invalid_argument("V85 invalid/incomplete physical state");
        if (!std::isfinite(merit_window_growth_cap) ||
            merit_window_growth_cap < 1.0 || merit_window_growth_cap > 1.05)
            throw std::invalid_argument("V85/V89 invalid bounded merit growth cap");
        if (!std::isfinite(tol.relative_tolerance) || tol.relative_tolerance<=0. ||
            !std::isfinite(tol.continuity_tolerance) || tol.continuity_tolerance<=0.)
            throw std::invalid_argument("V85 invalid final tolerances");
        if (!recent_.empty() && !same_operator)
            throw std::invalid_argument("V85 operator changed; refuse ledger comparison");
        GlobalOuterTrialV85 t;
        t.baseline=recent_.empty();
        if (t.baseline) {
            // First provisional state must still be finite and conservative.
            if (s.continuity>5.*tol.continuity_tolerance)
                return bad("initial_continuity_guard");
            t.allowed=true;t.reason="provisional_baseline";return t;
        }
        const auto& prev=recent_.back();
        if (s.algorithm!=prev.algorithm || s.schur!=prev.schur)
            throw std::invalid_argument("V85 algorithm/schur changed");
        const auto ratio=[](double a,double b){return a/std::max(b,1.e-300);};
        t.merit_ratio=ratio(s.merit,prev.merit);
        t.momentum_ratio=ratio(s.momentum.global_relative,prev.momentum.global_relative);
        t.k_ratio=ratio(s.k,prev.k);
        t.omega_ratio=ratio(s.omega,prev.omega);
        if (s.converged && s.sst_transport_committed) {
            // Do not trust asserted convergence; V8 and its independent actual
            // component residuals still govern the final certificate.
            if (s.momentum.max_active_component_relative>tol.relative_tolerance ||
                s.momentum.global_relative>tol.relative_tolerance ||
                s.continuity>tol.continuity_tolerance ||
                s.k>tol.relative_tolerance ||s.omega>tol.relative_tolerance ||
                s.eddy_viscosity_closure>tol.relative_tolerance)
                return bad("false_final_certificate");
            t.allowed=true;t.reason="v8_final_residual_verified";return t;
        }
        if (s.continuity>5.*tol.continuity_tolerance ||
            s.continuity>std::max(5.*tol.continuity_tolerance,10.*prev.continuity))
            return bad("continuity_guard");
        if (s.momentum.max_active_component_relative>1.02*prev.momentum.max_active_component_relative ||
            s.momentum.global_relative>1.02*prev.momentum.global_relative)
            return bad("momentum_guard");
        if (s.eddy_viscosity_closure>std::max(10.*tol.relative_tolerance,
                                                1.2*prev.eddy_viscosity_closure))
            return bad("nut_closure_guard");
        double best_merit=prev.merit,best_k=prev.k,best_omega=prev.omega;
        for (const auto& x:recent_) {
            best_merit=std::min(best_merit,x.merit);
            best_k=std::min(best_k,x.k);
            best_omega=std::min(best_omega,x.omega);
        }
        if (s.merit>merit_window_growth_cap*best_merit)
            return bad("merit_window_guard");
        const double k_scale=std::max(best_k,10.*tol.relative_tolerance);
        const double omega_scale=std::max(best_omega,10.*tol.relative_tolerance);
        if (s.k>sst_best_multiplier*k_scale ||
            s.k>sst_previous_multiplier*std::max(prev.k,10.*tol.relative_tolerance))
            return bad("k_excursion_guard");
        if (s.omega>sst_best_multiplier*omega_scale ||
            s.omega>sst_previous_multiplier*std::max(prev.omega,10.*tol.relative_tolerance))
            return bad("omega_excursion_guard");
        t.k_debt=s.k>1.2*k_scale?k_debt_+1:0;
        t.omega_debt=s.omega>1.2*omega_scale?omega_debt_+1:0;
        if (t.k_debt>=window || t.omega_debt>=window)
            return bad("sst_exhausted_recovery_window");
        const bool excursion=t.k_debt>0 || t.omega_debt>0;
        if (excursion && t.momentum_ratio>0.98)
            return bad("excursion_without_momentum_progress");
        // At least one independent residual family must have contracted;
        // never accept a trial solely because a monitor flag says 'progress'.
        if (!(s.momentum.global_relative<0.995*prev.momentum.global_relative ||
              s.k<0.98*prev.k ||s.omega<0.98*prev.omega ||
              s.merit<0.9999*best_merit))
            return bad("no_physical_progress");
        if (recent_.size()>=window &&
            s.merit>recent_.front().merit*(1.-1.e-4))
            return bad("no_net_window_merit_progress");
        t.allowed=true;
        t.nonmonotone=s.merit>prev.merit;
        t.reason=excursion?"bounded_provisional_sst_excursion":"bounded_progress";
        return t;
    }
    void commit(const GlobalNonlinearResidualV8& state,const ConvergenceCriteria& gates,
                bool same_operator=true,double merit_window_growth_cap=1.01) {
        const auto t=evaluate(state,gates,same_operator,merit_window_growth_cap);
        if (!t.allowed) throw std::logic_error("V85 rejected state cannot enter ledger");
        k_debt_=t.k_debt;omega_debt_=t.omega_debt;
        recent_.push_back(state);
        while (recent_.size()>window)recent_.pop_front();
        ++steps_;
    }
    std::size_t committed_steps()const noexcept{return steps_;}
    std::size_t k_debt()const noexcept{return k_debt_;}
    std::size_t omega_debt()const noexcept{return omega_debt_;}
private:
    static bool finite_nonnegative(double x){return std::isfinite(x)&&x>=0.;}
    std::deque<GlobalNonlinearResidualV8> recent_;
    std::size_t steps_=0,k_debt_=0,omega_debt_=0;
};
} // namespace cfdx::physics