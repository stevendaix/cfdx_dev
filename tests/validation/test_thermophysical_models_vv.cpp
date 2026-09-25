#include "cfdx/physics/thermophysical_models.h"
#include "cfdx/physics/turbulence_models.h"
#include "common/test_harness.h"
#include <cmath>
#include <iostream>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("thermophysical_vv_analytic_linear", [] {
        ScalarPropertyControls k;
        k.model = ScalarPropertyModel::LINEAR;
        k.reference_temperature = 300;
        k.reference_value = 10;
        k.coefficients = {0.25};
        EXPECT_NEAR(evaluate_scalar_property(k, 340), 20, 1e-12);
        EXPECT_NEAR(
            (evaluate_scalar_property(k, 340) - evaluate_scalar_property(k, 300)) / 40,
            0.25, 1e-12);
    });

    run_case("thermophysical_vv_table_exact_linear", [] {
        ScalarPropertyControls c;
        c.model = ScalarPropertyModel::TABLE;
        c.table.temperature = {300.0, 400.0, 500.0};
        c.table.value = {1000.0, 1200.0, 1400.0};
        EXPECT_NEAR(evaluate_scalar_property(c, 350), 1100, 1e-12);

        ThermophysicalProperties p;
        p.heat_capacity = c;
        EXPECT_NEAR(p.enthalpy(500, 300), 240000, 1e-10);

        c.table.extrapolation = ExtrapolationPolicy::REJECT;
        bool threw = false;
        try {
            (void)evaluate_scalar_property(c, 250);
        } catch (...) {
            threw = true;
        }
        EXPECT_TRUE(threw);
    });

    run_case("boussinesq_validity_vv", [] {
        DensityControls d;
        d.rho0 = 1000;
        d.model = DensityModel::BOUSSINESQ;
        d.beta = 1e-3;
        EXPECT_NEAR(evaluate_density(d, 330), 970, 1e-12);

        bool threw = false;
        try {
            (void)evaluate_density(d, 450);
        } catch (...) {
            threw = true;
        }
        EXPECT_TRUE(threw);
    });

    run_case("hybrid_shielding_vv", [] {
        EXPECT_NEAR(ddes_shielding(0), 1, 1e-12);
        EXPECT_TRUE(ddes_shielding(1) < 1e-10);
        EXPECT_NEAR(ddes_length_scale_from_rd(0.5, 0.2, 0.65, 0), 0.13, 1e-12);
        EXPECT_NEAR(ddes_length_scale_from_rd(0.5, 0.2, 0.65, 1), 0.5, 1e-12);
        EXPECT_NEAR(iddes_length_scale_from_rd(0.5, 0.2, 0.65, 0, 1), 0.13, 1e-12);
        EXPECT_NEAR(iddes_length_scale_from_rd(0.5, 0.2, 0.65, 1, 1), 0.5, 1e-12);
    });

    std::cout << "THERMAL_RESIDUAL: suite=thermophysical_models PASS\n";
    return run_all();
}
