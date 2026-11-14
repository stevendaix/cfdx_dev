#pragma once
// V32 experimental: relative, symmetric, bounded pseudo-time amplitudes
// for momentum and (optionally) the coupled k-omega JFNK direction. The
// external CFL is shared, while SST may adapt a transaction-local inner CFL;
// it is NOT dimensionally meaningful to equate two raw spectral
// masses with different physical units. Instead bound both by the SAME
// fraction of their own *stationary* FV diagonal. No physical equation,
// wall constraint, Rhie-Chow flux, or certification threshold is changed.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {
struct SSTV32CouplingBudget {
    // Conservative initial screening magnitude measured on the reduced BFS.
    // Not universal: must earn expansion through actual outer commits.
    double fraction=1.0e-6;
    std::size_t committed_outer=0;
    static constexpr double minimum_fraction=1.0e-8;
    static constexpr double maximum_fraction=4.0;
    void validate() const {
        if (!std::isfinite(fraction) || fraction<minimum_fraction ||
            fraction>maximum_fraction)
            throw std::invalid_argument("V32 outer PTC: invalid relative diagonal budget");
    }
    // Failed outer retries may use a smaller *trial* budget, but can never
    // publish it. The controller itself is unmodified until outer commit.
    double trial_fraction(std::size_t retries) const {
        validate();
        const int power=-static_cast<int>(std::min<std::size_t>(retries,32));
        return std::max(minimum_fraction,std::ldexp(fraction,power));
    }
    void committed(double old_merit, double new_merit,
                   bool globally_converged=false) {
        validate();
        if (!std::isfinite(old_merit) || !std::isfinite(new_merit) ||
            old_merit<0.0 || new_merit<0.0 ||
            (!globally_converged && !(new_merit<old_merit)))
            throw std::invalid_argument("V32 outer PTC: unverified physical contraction");
        if (!globally_converged && old_merit>0.0) {
            const double relative=new_merit/old_merit;
            const double gain=relative<0.8 ? 1.2 : relative>0.98 ? 0.75 : 1.0;
            fraction=std::clamp(fraction*gain,minimum_fraction,maximum_fraction);
        }
        ++committed_outer;
    }
};
} // namespace cfdx::physics
