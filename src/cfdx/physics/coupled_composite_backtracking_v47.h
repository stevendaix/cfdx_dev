#pragma once
// V47 optional reduced-space backtracking for a certified COUPLED carrier
// direction. Each beta is followed by a NEW nonlinear SST transport and a
// fully reassembled U,p,k,omega,nu_t merit (V7/V8/V13/V14) BEFORE commit.
// This is NOT a Newton six-field Jacobian or an interpolation of SST states.
// No field, flux or convergence certificate is published before acceptance.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

inline void validate_coupled_composite_backtracking_v47(
    std::size_t retries, double shrink) {
    if (retries < 1 || retries > 12 ||
        !std::isfinite(shrink) || shrink <= 0.0 || shrink >= 1.0 ||
        std::pow(shrink, static_cast<double>(retries)) < 1.0e-8)
        throw std::invalid_argument(
            "V47 backtracking requires 1..12 retries and a representable step >=1e-8");
}

inline double coupled_composite_beta_v47(std::size_t rejected_trials,
                                         double shrink) {
    if (!std::isfinite(shrink) || shrink <= 0.0 || shrink >= 1.0 ||
        rejected_trials > 12)
        throw std::invalid_argument("V47 invalid composite line-search controls");
    const double beta = std::pow(shrink, static_cast<double>(rejected_trials));
    if (!std::isfinite(beta) || beta < 1.0e-8 || beta > 1.0)
        throw std::invalid_argument("V47 composite line-search step below safety floor");
    return beta;
}

// Safeguarded nonlinear *trial* selection. After three rejected full-state
// physical evaluations, fit a quadratic to their V13 physical growth scores.
// The minimum is only a proposal: the solver must re-solve SST and pass V13
// AND V14 before committing any field or Rhie-Chow flux. If the data do not
// bracket a strict local minimum, continue ordinary geometric backtracking.
class CoupledCompositeSearchV47 {
public:
    double trial_beta(std::size_t retry, double shrink) const {
        const double geometric=coupled_composite_beta_v47(retry,shrink);
        if (retry!=3 || samples_!=3) return geometric;
        const double x0=beta_[0],x1=beta_[1],x2=beta_[2];
        const double y0=score_[0],y1=score_[1],y2=score_[2];
        if (!(x0>x1 && x1>x2)) return geometric;
        const double slope01=(y0-y1)/(x0-x1);
        const double slope12=(y1-y2)/(x1-x2);
        const double a=(slope01-slope12)/(x0-x2);
        const double b=slope01-a*(x0+x1);
        if (!std::isfinite(a) || !std::isfinite(b) || !(a>0.0))
            return geometric;
        const double proposal=-b/(2.0*a);
        const double margin=0.10*(x1-x2);
        if (!std::isfinite(proposal) ||
            proposal<=x2+margin || proposal>=x1-margin)
            return geometric;
        return proposal;
    }
    void observe_rejection(double beta, double true_physical_growth) {
        if (samples_>=3) return;
        if (!std::isfinite(beta) || beta<=0.0 || beta>1.0 ||
            !std::isfinite(true_physical_growth) || true_physical_growth<0.0)
            throw std::invalid_argument("V47 invalid physical growth sample");
        beta_[samples_]=beta;
        score_[samples_]=true_physical_growth;
        ++samples_;
    }
    std::size_t samples() const noexcept { return samples_; }
private:
    std::array<double,3> beta_{};
    std::array<double,3> score_{};
    std::size_t samples_ = 0;
};

} // namespace cfdx::physics
