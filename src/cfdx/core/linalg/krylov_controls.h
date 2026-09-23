#pragma once
#include "krylov_reductions.h"
#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace cfdx::core {

struct KrylovControls {
    int restart_min = 10;
    int restart_max = 40;
    bool adaptive_restart = true;
    std::size_t residual_replacement_period = 50;
    bool residual_replacement = true;
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
        c.target_reduction_per_cycle <= 0.0 || c.target_reduction_per_cycle >= 1.0 ||
        current <= 0) {
        throw std::invalid_argument("invalid Krylov controls");
    }
    if (!c.adaptive_restart) return std::clamp(current, c.restart_min, c.restart_max);
    if (cycle_reduction > c.target_reduction_per_cycle * 2.0)
        return std::max(c.restart_min, current - 5);
    if (cycle_reduction < c.target_reduction_per_cycle * 0.5)
        return std::min(c.restart_max, current + 5);
    return std::clamp(current, c.restart_min, c.restart_max);
}

} // namespace cfdx::core
