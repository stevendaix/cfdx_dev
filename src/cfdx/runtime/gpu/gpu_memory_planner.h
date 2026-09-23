#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::runtime::gpu {

struct MemoryBudget {
    std::size_t device_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    std::size_t runtime_reserve_bytes = 0;
};

inline std::size_t checked_add(std::size_t a, std::size_t b)
{
    if (b > std::numeric_limits<std::size_t>::max() - a)
        throw std::overflow_error("GPU memory budget overflow");
    return a + b;
}

inline bool fits_memory_budget(
    std::size_t required_bytes,
    const MemoryBudget& budget)
{
    const auto reserve = checked_add(
        budget.safety_margin_bytes, budget.runtime_reserve_bytes);
    if (reserve > budget.device_bytes)
        return false;
    return required_bytes <= budget.device_bytes - reserve;
}

inline std::size_t available_working_set_bytes(const MemoryBudget& budget)
{
    const auto reserve = checked_add(
        budget.safety_margin_bytes, budget.runtime_reserve_bytes);
    return reserve >= budget.device_bytes ? 0 : budget.device_bytes - reserve;
}

inline std::size_t plan_working_set(
    std::size_t requested_bytes,
    const MemoryBudget& budget)
{
    if (!fits_memory_budget(requested_bytes, budget))
        throw std::runtime_error("GPU working set exceeds memory budget");
    return requested_bytes;
}

} // namespace cfdx::runtime::gpu
