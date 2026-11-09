#pragma once
// V14 strict, opt-in contraction of the reassembled PHYSICAL U,p,k,omega,
// nu_t merit. Commit only after the complete carrier/SST state is accepted;
// retry candidates must not advance this ledger. No convergence tolerance is
// changed and a successful commit does not imply physical convergence.
#include "cfdx/physics/global_nonlinear_residual_v8.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct GlobalContractionTrialV14 {
    bool allowed=false;
    bool baseline=false;
    double old_merit=0.;
    double new_merit=0.;
    double ratio=0.;
};

class GlobalNonlinearContractionLedgerV14 {
public:
    explicit GlobalNonlinearContractionLedgerV14(double max_ratio=0.995)
        : max_ratio_(max_ratio) {
        if (!(max_ratio_>0. && max_ratio_<1.) || !std::isfinite(max_ratio_))
            throw std::invalid_argument("V14 merit contraction must be strictly below one");
    }

    GlobalContractionTrialV14 evaluate(const GlobalNonlinearResidualV8& trial,
                                      bool same_physical_operator=true) const {
        validate(trial);
        GlobalContractionTrialV14 result;
        result.baseline=!has_last_ || !same_physical_operator;
        result.new_merit=trial.merit;
        result.old_merit=result.baseline ? 0. : last_.merit;
        if (!result.baseline) {
            if (last_.algorithm != trial.algorithm || last_.schur != trial.schur ||
                last_.sst != trial.sst)
                throw std::invalid_argument("V14 compared different nonlinear operators");
            result.ratio=trial.merit/std::max(last_.merit,1e-300);
            // When all physical tests already pass, the physical convergence
            // result takes precedence over a meaningless machine-noise ratio.
            result.allowed=trial.converged ||
                (std::isfinite(result.ratio) && result.ratio<=max_ratio_);
        } else {
            result.ratio=0.;
            result.allowed=true;
        }
        return result;
    }

    void commit(const GlobalNonlinearResidualV8& accepted,
                bool same_physical_operator=true) {
        const auto trial=evaluate(accepted,same_physical_operator);
        if (!trial.allowed)
            throw std::logic_error("V14 rejected state cannot enter committed ledger");
        if (trial.baseline) contraction_steps_=0;
        else ++contraction_steps_;
        last_=accepted;
        has_last_=true;
        ++committed_steps_;
    }
    std::size_t committed_steps() const noexcept { return committed_steps_; }
    std::size_t contraction_steps() const noexcept { return contraction_steps_; }
    bool has_last() const noexcept { return has_last_; }

private:
    static void validate(const GlobalNonlinearResidualV8& state) {
        if (!state.all_finite || !state.sst || !state.bounded ||
            !state.sst_transport_committed || !std::isfinite(state.merit) ||
            state.merit<0.)
            throw std::invalid_argument("V14 expected a complete finite SST physical certificate");
    }
    double max_ratio_;
    GlobalNonlinearResidualV8 last_{};
    std::size_t committed_steps_=0;
    std::size_t contraction_steps_=0;
    bool has_last_=false;
};
} // namespace cfdx::physics
