#pragma once
#include "krylov_reductions.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::core {

struct KrylovControls {
    int restart_min = 10;
    int restart_max = 512;
    bool adaptive_restart = true;
    std::size_t residual_replacement_period = 50;
    bool residual_replacement = true;
    // Re-evaluate the UPDATED Krylov candidate periodically. Restart only if
    // the true original-system residual exceeds the Givens projection by this
    // factor; blindly restarting can degrade very ill-conditioned pressure
    // systems even when their Arnoldi relation is accurate.
    double residual_replacement_gap_ratio = 10.0;
    double target_reduction_per_cycle = 0.1;
    // When enabled, reductions are performed over the caller-owned local range
    // and then combined across ranks. This keeps solver math independent of
    // the distributed execution policy while providing deterministic mode.
    KrylovReductionPolicy reduction{};
};

inline int choose_gmres_restart(
    int current,
    double cycle_reduction,
    const KrylovControls& c = {}) {
    if (c.restart_min <= 0 || c.restart_min > c.restart_max ||
        !std::isfinite(c.target_reduction_per_cycle) ||
        c.target_reduction_per_cycle <= 0.0 || c.target_reduction_per_cycle >= 1.0 ||
        !std::isfinite(cycle_reduction) || cycle_reduction < 0.0 ||
        current <= 0) {
        throw std::invalid_argument("invalid Krylov controls");
    }
    if (!c.adaptive_restart) return std::clamp(current, c.restart_min, c.restart_max);
    // A cycle ratio close to zero is good convergence: a smaller restart is
    // sufficient. A ratio near one is poor convergence: retain more Krylov
    // directions. The previous implementation did the opposite and made
    // difficult coupled systems less robust.
    if (cycle_reduction > c.target_reduction_per_cycle * 2.0) {
        // Avoid signed overflow when current is close to INT_MAX.
        return current > c.restart_max - 5 ? c.restart_max : current + 5;
    }
    if (cycle_reduction < c.target_reduction_per_cycle * 0.5)
        return std::max(c.restart_min, current - 5);
    return std::clamp(current, c.restart_min, c.restart_max);
}

} // namespace cfdx::core
