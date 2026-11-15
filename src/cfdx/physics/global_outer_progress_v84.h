#pragma once
// V84: opt-in, bounded, finite-window progress policy for the FULL returned
// physical U,p,phi,k,omega,nu_t FV state. It NEVER creates a final certificate.
// Compare only identical operators, after reassembling momentum and both SST
// equations. Evaluate has NO side effects; commit only after the FV transaction.
#include "cfdx/physics/global_nonlinear_residual_v8.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

struct GlobalOuterTrialV84 {
    bool allowed = false;
    bool baseline = false;
    bool nonmonotone = false;
    double merit_ratio = 0.0;
    double momentum_global_ratio = 0.0;
    double k_ratio = 0.0;
    double omega_ratio = 0.0;
    std::size_t k_debt = 0;
    std::size_t omega_debt = 0;
    std::string reason;
};

class GlobalOuterProgressLedgerV84 {
public:
    // These are acceptance safeguards, NOT physical convergence tolerances.
    // The k/omega debt allows at most two accepted temporary SST excursions
    // before a trial must recover close to its best recent value.
    explicit GlobalOuterProgressLedgerV84(std::size_t window = 3,
                                          double sst_budget = 3.5,
                                          double merit_envelope = 1.01)
        : window_(window), sst_budget_(sst_budget),
          merit_envelope_(merit_envelope) {
        if (window_ < 2 || window_ > 8 ||
            !std::isfinite(sst_budget_) || sst_budget_ < 1.2 || sst_budget_ > 4. ||
            !std::isfinite(merit_envelope_) ||
            merit_envelope_ < 1. || merit_envelope_ > 1.03)
            throw std::invalid_argument("V84 invalid bounded progress policy");
    }

    GlobalOuterTrialV84 evaluate(const GlobalNonlinearResidualV8& trial,
                                 const ConvergenceCriteria& final_gates,
                                 bool same_physical_operator = true) const {
        validate(trial);
        if (!std::isfinite(final_gates.relative_tolerance) ||
            !std::isfinite(final_gates.continuity_tolerance) ||
            final_gates.relative_tolerance <= 0 ||
            final_gates.continuity_tolerance <= 0)
            throw std::invalid_argument("V84 invalid original convergence criteria");
        GlobalOuterTrialV84 verdict;
        verdict.baseline = recent_.empty();
        verdict.allowed = verdict.baseline;
        verdict.reason = verdict.baseline ? "first_provisional_baseline" : "ok";
        if (verdict.baseline) return verdict;
        if (!same_physical_operator)
            throw std::invalid_argument("V84 cannot compare different physical operators");
        const auto& prev = recent_.back();
        if (prev.algorithm != trial.algorithm || prev.schur != trial.schur)
            throw std::invalid_argument("V84 algorithm or Schur changed");
        const auto ratio=[](double a,double b){return a/std::max(b,1e-300);};
        verdict.merit_ratio=ratio(trial.merit,prev.merit);
        verdict.momentum_global_ratio=ratio(trial.momentum.global_relative,
                                            prev.momentum.global_relative);
        verdict.k_ratio=ratio(trial.k,prev.k);
        verdict.omega_ratio=ratio(trial.omega,prev.omega);
        const auto reject=[&](const char* reason){
            verdict.allowed=false; verdict.reason=reason; return verdict;
        };
        if (trial.converged && trial.sst_transport_committed) {
            // Never rely on an asserted convergence flag alone. V8 is the
            // authority, and all represented physical residuals must also
            // satisfy the ORIGINAL criteria on this returned state.
            if (trial.momentum.global_relative > final_gates.relative_tolerance ||
                trial.momentum.max_active_component_relative >
                    final_gates.relative_tolerance ||
                trial.continuity > final_gates.continuity_tolerance ||
                trial.k > final_gates.relative_tolerance ||
                trial.omega > final_gates.relative_tolerance ||
                trial.eddy_viscosity_closure > final_gates.relative_tolerance)
                return reject("inconsistent_final_certificate");
            verdict.allowed=true; verdict.reason="independent_final_certificate";
            return verdict;
        }
        // Pressure must remain conservative. A transient increase of SST may
        // never excuse a lost Rhie-Chow balance or momentum blow-up.
        if (trial.continuity > std::max(5.*final_gates.continuity_tolerance,
                                        10.*prev.continuity))
            return reject("pressure_continuity_guard");
        if (trial.momentum.max_active_component_relative >
                1.02*prev.momentum.max_active_component_relative ||
            trial.momentum.global_relative > 1.05*prev.momentum.global_relative)
            return reject("momentum_growth_guard");
        if (trial.eddy_viscosity_closure >
            std::max(10.*final_gates.relative_tolerance,
                     1.2*prev.eddy_viscosity_closure))
            return reject("eddy_viscosity_closure_guard");

        double best_merit=prev.merit;
        double best_k=prev.k;
        double best_omega=prev.omega;
        for (const auto& old:recent_) {
            best_merit=std::min(best_merit,old.merit);
            best_k=std::min(best_k,old.k);
            best_omega=std::min(best_omega,old.omega);
        }
        if (trial.merit > merit_envelope_*best_merit)
            return reject("nonmonotone_merit_envelope");
        // The best values come from ACTUALLY COMMITTED states, not line-search
        // trials. Do not rebase on a rejected state to manufacture progress.
        const double k_floor=10.*final_gates.relative_tolerance;
        const double omega_floor=10.*final_gates.relative_tolerance;
        if (trial.k > sst_budget_*std::max(best_k,k_floor) ||
            trial.k > sst_budget_*std::max(prev.k,k_floor))
            return reject("k_growth_budget");
        if (trial.omega > sst_budget_*std::max(best_omega,omega_floor) ||
            trial.omega > sst_budget_*std::max(prev.omega,omega_floor))
            return reject("omega_growth_budget");
        verdict.k_debt = trial.k > 1.20*std::max(best_k,k_floor)
                            ? k_debt_ + 1 : 0;
        verdict.omega_debt = trial.omega > 1.20*std::max(best_omega,omega_floor)
                                ? omega_debt_ + 1 : 0;
        if (verdict.k_debt >= window_ || verdict.omega_debt >= window_)
            return reject("unrecovered_sst_growth_debt");
        // A temporarily worse SST component is possible ONLY with a clearly
        // improving independently-reassembled global momentum residual.
        const bool sst_excursion=trial.k > 1.20*std::max(best_k,k_floor) ||
                                 trial.omega > 1.20*std::max(best_omega,omega_floor);
        if (sst_excursion && verdict.momentum_global_ratio > 0.98)
            return reject("sst_excursion_without_carrier_progress");
        // No permanent nonmonotone drift. Once a 3-state window is full,
        // the new merit must be lower than its oldest committed anchor.
        if (recent_.size() >= window_ &&
            trial.merit > recent_.front().merit*(1.-1e-4))
            return reject("no_net_merit_improvement_over_window");
        const bool merit_progress = trial.merit < best_merit*(1.-1e-4);
        const bool carrier_progress = verdict.momentum_global_ratio <= 0.995;
        const bool sst_progress = trial.k < prev.k*0.98 ||
                                  trial.omega < prev.omega*0.98;
        if (!merit_progress && !carrier_progress && !sst_progress)
            return reject("no_independent_physical_progress");
        verdict.allowed=true;
        verdict.nonmonotone=trial.merit>prev.merit;
        verdict.reason=verdict.nonmonotone ? "bounded_nonmonotone_progress" :
                                            "bounded_progress";
        return verdict;
    }

    void commit(const GlobalNonlinearResidualV8& state,
                const ConvergenceCriteria& gates,
                bool same_physical_operator = true) {
        const auto decision=evaluate(state,gates,same_physical_operator);
        if (!decision.allowed)
            throw std::logic_error("V84 rejected state cannot enter committed ledger");
        if (!recent_.empty() && (state.converged || !same_physical_operator)) {
            // In the current experimental driver operator changes are refused;
            // never silently carry a nonmonotone debt to another problem.
            k_debt_=omega_debt_=0;
        } else {
            k_debt_=decision.k_debt;
            omega_debt_=decision.omega_debt;
        }
        recent_.push_back(state);
        while (recent_.size()>window_) recent_.pop_front();
        ++committed_steps_;
    }
    std::size_t committed_steps() const noexcept {return committed_steps_;}
    std::size_t k_debt() const noexcept {return k_debt_;}
    std::size_t omega_debt() const noexcept {return omega_debt_;}
private:
    static void validate(const GlobalNonlinearResidualV8& s) {
        if (!s.all_finite || !s.sst || !s.bounded ||
            (s.converged && !s.sst_transport_committed) ||
            s.algorithm!=PressureVelocityAlgorithm::SIMPLE ||
            !std::isfinite(s.merit) || s.merit<0. ||
            !std::isfinite(s.continuity) || s.continuity<0. ||
            !std::isfinite(s.k) || s.k<0. ||
            !std::isfinite(s.omega) || s.omega<0. ||
            !std::isfinite(s.eddy_viscosity_closure) || s.eddy_viscosity_closure<0. ||
            !std::isfinite(s.momentum.global_relative) ||
            s.momentum.global_relative<0. ||
            !std::isfinite(s.momentum.max_active_component_relative) ||
            s.momentum.max_active_component_relative<0.)
            throw std::invalid_argument("V84 incomplete/unbounded/fictitious physical state");
    }
    std::size_t window_;
    double sst_budget_;
    double merit_envelope_;
    std::deque<GlobalNonlinearResidualV8> recent_;
    std::size_t committed_steps_=0;
    std::size_t k_debt_=0;
    std::size_t omega_debt_=0;
};
} // namespace cfdx::physics
