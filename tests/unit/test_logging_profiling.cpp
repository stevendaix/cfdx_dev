#include "cfdx/utils/logging_profiling.h"
#include "common/test_harness.h"

#include <cmath>
#include <stdexcept>

using namespace cfdx::testing;

int main()
{
    run_case("logging_residual_accepts_finite_values", [] {
        EXPECT_NO_THROW(cfdx::utils::log_residual("1e-8"));
        EXPECT_NO_THROW(cfdx::utils::log_residual("1e-3", 1e-2));
    });

    run_case("logging_residual_rejects_nonfinite_values", [] {
        EXPECT_THROW(cfdx::utils::log_residual("nan"), std::invalid_argument);
        EXPECT_THROW(cfdx::utils::log_residual("inf"), std::invalid_argument);
    });

    run_case("profiling_hook_measures_void_callable", [] {
        const auto elapsed = cfdx::utils::ProfilingHook::profile([] {});
        EXPECT_TRUE(elapsed.count() >= 0);
    });

    return run_all();
}
