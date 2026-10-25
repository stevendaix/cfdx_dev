#include "cfdx/physics/vof.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("planar_interface_has_zero_curvature", [] {
        EXPECT_NEAR(interface_curvature_from_normal_divergence(0.0, 0.0, 0.0),
                    0.0, 1e-14);
    });

    run_case("circle_curvature_uses_transverse_derivative", [] {
        // Outward radial normal at (R,0): d(nx)/dx=0, d(ny)/dy=1/R.
        constexpr double radius = 2.0;
        EXPECT_NEAR(interface_curvature_from_normal_divergence(
                        0.0, 1.0 / radius),
                    -1.0 / radius, 1e-14);
    });

    run_case("sphere_curvature_sums_all_three_axis_derivatives", [] {
        constexpr double radius = 4.0;
        EXPECT_NEAR(interface_curvature_from_normal_divergence(
                        0.0, 1.0 / radius, 1.0 / radius),
                    -2.0 / radius, 1e-14);
    });

    run_case("curvature_rejects_non_finite_derivatives", [] {
        EXPECT_THROW(interface_curvature_from_normal_divergence(
                         std::numeric_limits<double>::quiet_NaN(), 0.0),
                     std::invalid_argument);
    });

    return run_all();
}
