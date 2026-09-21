#include "cfdx/physics/sst_solver.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::testing;
using namespace cfdx::physics;

int main()
{
    run_case("sst_blending_is_bounded", [] {
        const auto [f1, f2] = compute_sst_blending(0.5, 10.0, 0.01, 1.5e-5, 0.09);
        EXPECT_TRUE(f1 >= 0.0 && f1 <= 1.0);
        EXPECT_TRUE(f2 >= 0.0 && f2 <= 1.0);
    });

    run_case("sst_blending_responds_to_wall_distance", [] {
        const auto near_wall = compute_sst_blending(0.5, 10.0, 1.0e-3, 1.5e-5, 0.09);
        const auto far_wall = compute_sst_blending(0.5, 10.0, 10.0, 1.5e-5, 0.09);
        EXPECT_TRUE(near_wall.first > far_wall.first);
        EXPECT_TRUE(near_wall.second > far_wall.second);
    });

    run_case("sst_blending_rejects_invalid_distance", [] {
        EXPECT_THROW(compute_sst_blending(0.5, 10.0, 0.0, 1.5e-5, 0.09), std::invalid_argument);
    });

    return run_all();
}
