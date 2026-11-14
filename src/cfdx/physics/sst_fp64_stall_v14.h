#pragma once
// An identical returned FP64 (k,omega) pair and an unmet independent true
// equation gate cannot improve by repeating the same deterministic assembly.
// Stop after three repeats, preserving UNCONVERGED status and ALL tolerances.
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct SSTFp64StallV14 {
    std::size_t identical_unconverged=0;
    bool observe(double dk, double domega, bool all_gates_pass) {
        if (!std::isfinite(dk) || !std::isfinite(domega) || dk < 0 || domega < 0)
            throw std::invalid_argument("SST FP64 stall: invalid update");
        if (all_gates_pass || dk!=0. || domega!=0.) identical_unconverged=0;
        else ++identical_unconverged;
        return identical_unconverged>=3;
    }
};

} // namespace cfdx::physics
