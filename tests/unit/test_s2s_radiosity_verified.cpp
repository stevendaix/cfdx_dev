#include "cfdx/physics/radiation_s2s.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("s2s_rejects_nonrepresentable_blackbody_emission", [] {
        const std::vector<double> area{1.0};
        const std::vector<double> emissivity{1.0};
        const std::vector<double> temperature{1.0e100};
        const std::vector<double> F{0.0};
        EXPECT_THROW_WITH(solve_s2s_radiosity(area, emissivity, temperature, F),
                          std::runtime_error,
                          "S2S: non-representable blackbody emission");
    });
    run_case("s2s_rejects_nonrepresentable_integrated_power", [] {
        const std::vector<double> area{1.0e308};
        const std::vector<double> emissivity{1.0};
        const std::vector<double> temperature{1000.0};
        const std::vector<double> F{0.0};
        EXPECT_THROW_WITH(solve_s2s_radiosity(area, emissivity, temperature, F),
                          std::runtime_error,
                          "S2S: non-representable irradiation or surface power");
    });
    run_case("s2s_reports_true_radiosity_backward_error", [] {
        const std::vector<double> area{1.0, 1.0};
        const std::vector<double> emissivity{0.75, 0.45};
        const std::vector<double> temperature{650.0, 470.0};
        const std::vector<double> F{0.0, 1.0, 1.0, 0.0};
        const auto r = solve_s2s_radiosity(area, emissivity, temperature, F);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(std::isfinite(r.radiosity_backward_error));
        EXPECT_TRUE(r.radiosity_backward_error <= 1e-10);
        const double e0 = STEFAN_BOLTZMANN * std::pow(temperature[0], 4);
        const double e1 = STEFAN_BOLTZMANN * std::pow(temperature[1], 4);
        const double lhs0 = r.radiosity[0] - (1.0 - emissivity[0]) * r.radiosity[1];
        const double lhs1 = r.radiosity[1] - (1.0 - emissivity[1]) * r.radiosity[0];
        EXPECT_NEAR(lhs0, emissivity[0] * e0, 1e-10 * e0);
        EXPECT_NEAR(lhs1, emissivity[1] * e1, 1e-10 * e1);
    });
    run_case("s2s_user_tolerance_blocks_false_direct_solver_pass", [] {
        S2SControls strict;
        strict.tolerance = 1e-20; // more demanding than roundoff in this case
        const auto r = solve_s2s_radiosity(
            {1.0, 1.0}, {0.75, 0.45}, {650.0, 470.0},
            {0.0, 1.0, 1.0, 0.0}, strict);
        // A direct solve may be usable physically, but its true residual
        // must not be silently promoted to an unachievable tighter gate.
        EXPECT_FALSE(r.converged);
        EXPECT_TRUE(r.radiosity_backward_error > strict.tolerance);
    });
    run_case("s2s_open_enclosure_flux_can_be_nonzero", [] {
        const auto r=solve_s2s_radiosity({1.0}, {1.0}, {500.0}, {0.0});
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(r.total_power>0.0);
    });
    return run_all();
}
