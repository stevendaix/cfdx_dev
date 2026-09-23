#include "cfdx/runtime/execution/memory_budget.h"
#include "common/test_harness.h"

#include <limits>
#include <stdexcept>

using namespace cfdx::testing;
using namespace cfdx::runtime;

int main() {
    run_case("hardware_free_memory_is_used_as_budget", [] {
        MemoryBudgetHardware hw{16000, 9000, true};
        MemoryBudgetPolicy policy;
        const auto b = make_memory_budget(hw, policy, 8000);
        EXPECT_TRUE(b.hardware_budget_bytes == 9000);
        EXPECT_TRUE(b.budget_bytes == 9000);
        EXPECT_TRUE(b.fits);
        EXPECT_TRUE(b.usable_bytes == 9000);
    });

    run_case("configured_budget_caps_hardware_budget", [] {
        MemoryBudgetHardware hw{16000, 12000, true};
        MemoryBudgetPolicy policy;
        policy.configured_budget_bytes = 10000;
        policy.runtime_reserve_bytes = 500;
        policy.safety_margin_bytes = 500;
        const auto b = make_memory_budget(hw, policy, 9000);
        EXPECT_TRUE(b.budget_bytes == 10000);
        EXPECT_TRUE(b.required_bytes == 10000);
        EXPECT_TRUE(b.usable_bytes == 9000);
        EXPECT_TRUE(b.fits);
    });

    run_case("fractional_margin_is_explicit_and_bounded", [] {
        MemoryBudgetHardware hw{10000, 10000, true};
        MemoryBudgetPolicy policy;
        policy.safety_margin_fraction = 0.10;
        const auto b = make_memory_budget(hw, policy, 4000);
        EXPECT_TRUE(b.safety_margin_bytes == 400);
        EXPECT_TRUE(b.required_bytes == 4400);
        EXPECT_TRUE(b.usable_bytes == 9600);
        EXPECT_TRUE(b.fits);
    });

    run_case("budget_exhaustion_is_fail_fast", [] {
        MemoryBudgetHardware hw{10000, 5000, true};
        MemoryBudgetPolicy policy;
        policy.runtime_reserve_bytes = 500;
        const auto b = make_memory_budget(hw, policy, 5000);
        EXPECT_TRUE(!b.fits);
        EXPECT_TRUE(b.usable_bytes == 4500);
    });

    run_case("unknown_free_memory_uses_total_capacity", [] {
        MemoryBudgetHardware hw{16000, 0, false};
        const auto b = make_memory_budget(hw, {}, 15000);
        EXPECT_TRUE(b.hardware_budget_bytes == 16000);
        EXPECT_TRUE(b.fits);
    });

    run_case("invalid_fraction_is_rejected", [] {
        MemoryBudgetPolicy policy;
        policy.safety_margin_fraction = 1.01;
        EXPECT_THROW(
            make_memory_budget(MemoryBudgetHardware{100, 100, true},
                               policy, 1),
            std::invalid_argument);
    });

    run_case("required_bytes_overflow_is_rejected", [] {
        MemoryBudgetPolicy policy;
        policy.runtime_reserve_bytes = 1;
        const auto b = make_memory_budget(
            MemoryBudgetHardware{std::numeric_limits<std::size_t>::max(),
                                 std::numeric_limits<std::size_t>::max(), true},
            policy, std::numeric_limits<std::size_t>::max());
        EXPECT_TRUE(b.overflow);
        EXPECT_TRUE(!b.fits);
    });

    run_case("workload_margins_are_not_double_counted", [] {
        MemoryBudgetHardware hw{10000, 9000, true};
        MemoryBudgetPolicy policy;
        policy.safety_margin_fraction = 0.50;
        policy.runtime_reserve_bytes = 5000;
        RuntimeWorkload workload{3000, 500, 500};
        const auto b = make_memory_budget(hw, policy, workload);
        EXPECT_TRUE(b.required_bytes == 4000);
        EXPECT_TRUE(b.fits);
        EXPECT_TRUE(b.safety_margin_bytes == 500);
        EXPECT_TRUE(b.runtime_reserve_bytes == 500);
    });

    return run_all();
}
