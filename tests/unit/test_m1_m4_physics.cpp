#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/turbulence.h"
#include "cfdx/physics/thermal.h"
#include "cfdx/physics/radiation.h"
#include "common/test_harness.h"
#include <cmath>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("pressure_velocity_relaxation_and_rhie_chow", [] {
        CouplingControls c;
        validate_coupling_controls(c);
        EXPECT_NEAR(relaxed_value(1.0, 3.0, 0.5), 2.0, 1e-12);
        EXPECT_NEAR(rhie_chow_face_flux(2.0, 5.0, 3.0, 1.0, 2.0, 2.0, 1.0),
                    1.5, 1e-12);
        EXPECT_TRUE(piso_correction_gain(2.0, 0.5) > 0.0);
    });

    run_case("rans_and_les_models_are_positive", [] {
        EXPECT_NEAR(turbulent_kinematic_viscosity_kepsilon(4.0, 2.0), 0.72, 1e-12);
        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(1.0, 2.0, 0.5) > 0.0);
        EXPECT_NEAR(smagorinsky_eddy_viscosity(0.1, 20.0),
                    0.1156, 1e-12);
        EXPECT_TRUE(des_eddy_viscosity(0.1, 0.02, 20.0) >= 0.0);
    });

    run_case("thermal_and_cht_fluxes", [] {
        EXPECT_NEAR(heat_flux_conduction(10.0, 300.0, 290.0, 0.5),
                    200.0, 1e-12);
        const double g = cht_interface_conductance(10.0, 2.0, 0.01, 0.02, 1.0);
        EXPECT_NEAR(g, 83.33333333333333, 1e-10);
        EXPECT_NEAR(cht_interface_heat_flux(g, 320.0, 300.0),
                    1666.6666666666665, 1e-9);
    });

    run_case("surface_radiation_and_view_factors", [] {
        EXPECT_NEAR(blackbody_emissive_power(0.0), 0.0, 1e-12);
        EXPECT_NEAR(M_PI * blackbody_intensity(300.0), blackbody_emissive_power(300.0), 1e-10);
        const double q = gray_surface_emissivity_flux(1.0, 300.0, 0.0);
        EXPECT_NEAR(q, STEFAN_BOLTZMANN * std::pow(300.0, 4), 1e-10);
        EXPECT_TRUE(two_surface_net_exchange(1.0, 1.0, 400.0, 300.0, 1.0) > 0.0);
        EXPECT_NEAR(two_surface_net_exchange(1.0, 1.0, 400.0, 300.0, 0.0), 0.0, 1e-12);

        std::vector<double> F{0.0, 1.0, 1.0, 0.0};
        validate_view_factor_matrix(F, 2);
        validate_discrete_directions({
            {1.0, 0.0, 0.0, 2.0 * M_PI},
            {-1.0, 0.0, 0.0, 2.0 * M_PI}
        });
    });

    return run_all();
}
