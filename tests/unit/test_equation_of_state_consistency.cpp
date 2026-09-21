#include "cfdx/physics/equation_of_state.h"
#include "common/test_harness.h"

#include <cmath>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("ideal_gas_thermo_identities_are_consistent", [] {
        IdealGasEOS eos;
        const double p = 101325.0;
        const double T = 450.0;
        const double rho = eos.density(p, T);
        const double h = eos.enthalpy(p, T);
        const double recovered_T = eos.temperature_from_enthalpy(p, h);

        EXPECT_NEAR(eos.pressure_from_density_temp(rho, T), p, 1e-10);
        EXPECT_NEAR(recovered_T, T, 1e-12);
        EXPECT_NEAR(eos.cp(p, T) - eos.cv(p, T), 8.314462618 / 0.02896546, 1e-12);
        EXPECT_NEAR(eos.dp_drho_s(p, T), 1.4 * p / rho, 1e-12);
        EXPECT_NEAR(eos.dp_dT_rho(p, T), rho * (8.314462618 / 0.02896546), 1e-12);
    });

    run_case("ideal_gas_rejects_inconsistent_cp", [] {
        IdealGasParams params;
        params.Cp *= 1.01;
        EXPECT_THROW((IdealGasEOS(params)), std::invalid_argument);
    });

    run_case("ideal_gas_set_params_preserves_consistency_contract", [] {
        IdealGasEOS eos;
        IdealGasParams params;
        params.T_ref = 500.0;
        params.p_ref = 90000.0;
        eos.set_params(params);
        EXPECT_NEAR(eos.enthalpy(90000.0, 550.0), params.Cp * 50.0, 1e-12);
        EXPECT_NEAR(eos.temperature_from_enthalpy(90000.0, params.Cp * 50.0), 550.0, 1e-12);
    });

    return run_all();
}
