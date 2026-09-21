#include "cfdx/thermodynamics/thermo_state.h"
#include "cfdx/thermodynamics/thermo_cache.h"
#include "cfdx/transport/conductivity/conductivity.h"
#include "cfdx/transport/diffusivity/diffusivity.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace cfdx::testing;

int main()
{
    using namespace cfdx::thermodynamics;

    run_case("ideal_gas_default_is_thermodynamically_consistent", [] {
        IdealGasThermoModel model;
        const auto s = model.state(101325.0, 300.0);
        EXPECT_TRUE(s.is_valid());
        EXPECT_NEAR(s.cp, model.gamma * model.R / (model.gamma - 1.0), 1e-12);
        EXPECT_TRUE(s.speed_of_sound > 0.0);
    });

    run_case("ideal_gas_rejects_invalid_state", [] {
        IdealGasThermoModel model;
        EXPECT_THROW(model.state(0.0, 300.0), std::invalid_argument);
        EXPECT_THROW(model.state(101325.0, 0.0), std::invalid_argument);
        EXPECT_THROW(model.state(std::numeric_limits<double>::quiet_NaN(), 300.0),
                     std::invalid_argument);
    });

    run_case("ideal_gas_rejects_inconsistent_cp", [] {
        IdealGasThermoModel model;
        model.cp += 1.0;
        EXPECT_THROW(model.state(101325.0, 300.0), std::invalid_argument);
    });

    run_case("thermo_cache_update_is_transactional", [] {
        IdealGasThermoModel model;
        ThermoCache cache;
        cache.update(model, {101325.0}, {300.0});
        EXPECT_TRUE(cache.valid);
        const double old_rho = cache.state[0].rho;
        EXPECT_THROW(cache.update(model, {101325.0, 101325.0}, {300.0, 0.0}),
                     std::invalid_argument);
        EXPECT_TRUE(cache.valid);
        EXPECT_TRUE(cache.state.size() == 1);
        EXPECT_NEAR(cache.state[0].rho, old_rho, 1e-14);
    });

    run_case("conductivity_contracts", [] {
        EXPECT_NEAR(cfdx::transport::conductivity::constant(2.5), 2.5, 1e-14);
        EXPECT_NEAR(cfdx::transport::conductivity::prandtl(1.8e-5, 1004.5, 0.71),
                    1.8e-5 * 1004.5 / 0.71, 1e-14);
        EXPECT_THROW(cfdx::transport::conductivity::constant(-1.0), std::invalid_argument);
        EXPECT_THROW(cfdx::transport::conductivity::prandtl(0.0, 1004.5, 0.71),
                     std::invalid_argument);
        EXPECT_THROW(cfdx::transport::conductivity::prandtl(
                         1.8e-5, 1004.5, std::numeric_limits<double>::infinity()),
                     std::invalid_argument);
    });

    run_case("diffusivity_contracts", [] {
        EXPECT_NEAR(cfdx::transport::diffusivity::constant(1e-5), 1e-5, 1e-14);
        EXPECT_NEAR(cfdx::transport::diffusivity::schmidt(1.8e-5, 1.2, 0.7),
                    1.8e-5 / (1.2 * 0.7), 1e-14);
        EXPECT_THROW(cfdx::transport::diffusivity::constant(-1.0), std::invalid_argument);
        EXPECT_THROW(cfdx::transport::diffusivity::schmidt(1.8e-5, 0.0, 0.7),
                     std::invalid_argument);
        EXPECT_THROW(cfdx::transport::diffusivity::schmidt(
                         1.8e-5, 1.2, std::numeric_limits<double>::quiet_NaN()),
                     std::invalid_argument);
    });

    return run_all();
}
