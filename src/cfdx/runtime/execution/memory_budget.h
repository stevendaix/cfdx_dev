#pragma once

#include "cfdx/runtime/execution/execution_policy.h"

#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::runtime {

struct MemoryBudgetPolicy {
    // 0 means: use the hardware-reported available device memory.
    std::size_t configured_budget_bytes = 0;
    std::size_t runtime_reserve_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    double safety_margin_fraction = 0.0;
};

struct MemoryBudgetHardware {
    std::size_t total_device_bytes = 0;
    std::size_t free_device_bytes = 0;
    bool free_memory_known = false;
};

struct MemoryBudgetResult {
    std::size_t hardware_budget_bytes = 0;
    std::size_t budget_bytes = 0;
    std::size_t estimated_bytes = 0;
    std::size_t runtime_reserve_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    std::size_t required_bytes = 0;
    std::size_t usable_bytes = 0;
    bool overflow = false;
    bool fits = false;
};

inline MemoryBudgetResult make_memory_budget(
    const MemoryBudgetHardware& hardware,
    const MemoryBudgetPolicy& policy,
    std::size_t estimated_bytes)
{
    if (!(policy.safety_margin_fraction >= 0.0 &&
          policy.safety_margin_fraction <= 1.0)) {
        throw std::invalid_argument(
            "memory safety margin fraction must be in [0,1]");
    }

    MemoryBudgetResult result;
    result.estimated_bytes = estimated_bytes;
    result.runtime_reserve_bytes = policy.runtime_reserve_bytes;

    // A reported free-memory value is the immediately usable hardware budget.
    // If it is unavailable, total device memory is the conservative fallback
    // for the policy calculation; actual allocation must still be checked by
    // the backend.
    result.hardware_budget_bytes =
        hardware.free_memory_known
            ? hardware.free_device_bytes
            : hardware.total_device_bytes;

    result.budget_bytes = result.hardware_budget_bytes;
    if (policy.configured_budget_bytes != 0)
        result.budget_bytes =
            result.budget_bytes == 0
                ? policy.configured_budget_bytes
                : (policy.configured_budget_bytes < result.budget_bytes
                       ? policy.configured_budget_bytes
                       : result.budget_bytes);

    const auto max_size = std::numeric_limits<std::size_t>::max();
    const std::size_t proportional = static_cast<std::size_t>(
        static_cast<long double>(estimated_bytes) *
        policy.safety_margin_fraction);
    result.safety_margin_bytes =
        proportional > policy.safety_margin_bytes
            ? proportional
            : policy.safety_margin_bytes;

    if (estimated_bytes > max_size - result.runtime_reserve_bytes ||
        estimated_bytes + result.runtime_reserve_bytes >
            max_size - result.safety_margin_bytes) {
        result.overflow = true;
        result.required_bytes = max_size;
        result.usable_bytes = 0;
        return result;
    }

    result.required_bytes = estimated_bytes +
                            result.runtime_reserve_bytes +
                            result.safety_margin_bytes;

    if (result.budget_bytes >= result.runtime_reserve_bytes +
                                  result.safety_margin_bytes) {
        result.usable_bytes = result.budget_bytes -
                             result.runtime_reserve_bytes -
                             result.safety_margin_bytes;
    } else {
        result.usable_bytes = 0;
    }

    result.fits = !result.overflow &&
                  result.required_bytes <= result.budget_bytes;
    return result;
}

inline MemoryBudgetResult make_memory_budget(
    const MemoryBudgetHardware& hardware,
    const MemoryBudgetPolicy& policy,
    const RuntimeWorkload& workload)
{
    // RuntimeWorkload already carries its explicit reserve and margin. The
    // policy is used only to impose/override the hardware/configured budget.
    MemoryBudgetResult result;
    result.estimated_bytes = workload.estimated_bytes;
    result.runtime_reserve_bytes = workload.runtime_reserve_bytes;
    result.safety_margin_bytes = workload.safety_margin_bytes;
    result.hardware_budget_bytes =
        hardware.free_memory_known
            ? hardware.free_device_bytes
            : hardware.total_device_bytes;
    result.budget_bytes = result.hardware_budget_bytes;
    if (policy.configured_budget_bytes != 0)
        result.budget_bytes =
            result.budget_bytes == 0
                ? policy.configured_budget_bytes
                : (policy.configured_budget_bytes < result.budget_bytes
                       ? policy.configured_budget_bytes
                       : result.budget_bytes);

    const auto max_size = std::numeric_limits<std::size_t>::max();
    if (result.estimated_bytes > max_size - result.runtime_reserve_bytes ||
        result.estimated_bytes + result.runtime_reserve_bytes >
            max_size - result.safety_margin_bytes) {
        result.overflow = true;
        result.required_bytes = max_size;
        return result;
    }

    result.required_bytes = result.estimated_bytes +
                            result.runtime_reserve_bytes +
                            result.safety_margin_bytes;
    result.usable_bytes =
        result.budget_bytes >=
                result.runtime_reserve_bytes + result.safety_margin_bytes
            ? result.budget_bytes - result.runtime_reserve_bytes -
                  result.safety_margin_bytes
            : 0;
    result.fits = !result.overflow &&
                  result.required_bytes <= result.budget_bytes;
    return result;
}

} // namespace cfdx::runtime
