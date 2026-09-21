#include "cfdx/physics/m1_m4_models.h"
#include "common/test_harness.h"
#include <cmath>
#include <vector>

using namespace cfdx::physics::m1m4;
using namespace cfdx::testing;

int main()
{
    run_case("M1_coupling_controls_and_relaxation", [] {
        CouplingOptions c;
        validate_coupling(c);
        EXPECT_NEAR(under_relax(2.0, 6.0, 0.25), 3.0, 1e-14);
        EXPECT_NEAR(pressure_velocity_coefficient(2.0, 4.0), 0.5, 1e-14);
        EXPECT_THROW(validate_coupling(CouplingOptions{1, 1, -1, 0.7, 0.3, false}),
                     std::invalid_argument);
    });

    run_case("M1_Rhie_Chow_zero_pressure_correction", [] {
        const double phi = rhie_chow_mass_flux(
            1.0, 2.5, 10.0, 10.0, 0.0, 0.0, 0.5, 1.0);
        EXPECT_NEAR(phi, 2.5, 1e-14);
    });

    run_case("M2_k_epsilon_and_Smagorinsky", [] {
        EXPECT_NEAR(k_epsilon_nut(4.0, 2.0), 0.72, 1e-14);
        EXPECT_TRUE(k_epsilon_production(0.5, 2.0) > 0.0);
        EXPECT_TRUE(k_epsilon_dissipation_source(1.0, 1.92, 0.5, 2.0) > 0.0);
        EXPECT_NEAR(smagorinsky_nut(0.1, 20.0), 0.00578, 1e-14);
        EXPECT_NEAR(des_length_scale(0.1, 0.02), 0.02, 1e-14);
    });

    run_case("M2_SST_limiter_and_wall_law", [] {
        const double nut = sst_nut(1.0, 2.0, 5.0, 0.5);
        EXPECT_TRUE(nut > 0.0);
        EXPECT_NEAR(van_driest_damping(0.0), 0.0, 1e-14);
        EXPECT_NEAR(wall_function_u_plus(5.0), 5.0, 1e-14);
        EXPECT_TRUE(wall_function_u_plus(100.0) > wall_function_u_plus(20.0));
    });

    run_case("M3_conduction_and_CHT", [] {
        EXPECT_NEAR(thermal_diffusivity(10.0, 2.0, 5.0), 1.0, 1e-14);
        EXPECT_NEAR(conductive_flux(10.0, 300.0, 290.0, 0.5), 200.0, 1e-14);
        EXPECT_NEAR(interface_conductance(10.0, 2.0, 0.01, 0.02, 1.0),
                    90.9090909090909, 1e-12);
        EXPECT_NEAR(interface_heat_flux(10.0, 320.0, 300.0), 200.0, 1e-14);
        double Su=0.0, Sp=0.0;
        energy_source_linearization(100.0, 5.0, 20.0, Su, Sp);
        EXPECT_NEAR(Sp, 0.0, 1e-14);
        EXPECT_NEAR(Su, 100.0, 1e-14);
    });

    run_case("M4_blackbody_gray_and_exchange", [] {
        EXPECT_NEAR(blackbody(0.0), 0.0, 1e-14);
        EXPECT_NEAR(gray_emission(1.0, 300.0), blackbody(300.0), 1e-12);
        EXPECT_TRUE(two_surface_exchange(1.0, 1.0, 400.0, 300.0, 1.0) > 0.0);
        EXPECT_NEAR(two_surface_exchange(1.0, 1.0, 300.0, 300.0, 1.0), 0.0, 1e-14);
    });

    run_case("M4_view_factor_reciprocity_and_DOM", [] {
        const std::vector<double> F{0.0, 1.0, 1.0, 0.0};
        validate_view_factors(F, 2);
        validate_reciprocity(F, {1.0, 1.0}, 2);
        const double pi = std::acos(-1.0);
        validate_dom({
            {1.0, 0.0, 0.0, 2.0*pi},
            {-1.0, 0.0, 0.0, 2.0*pi}
        });
        EXPECT_THROW(validate_view_factors({0.0, 0.2, 0.5, 0.0}, 2),
                     std::invalid_argument);
    });

    run_case("M1_poiseuille_exact_oracle", [] {
        const double H=1.0, mu=2.0, dpdx=-4.0;
        EXPECT_NEAR(plane_poiseuille_velocity(0.0,H,dpdx,mu),0.0,1e-14);
        EXPECT_NEAR(plane_poiseuille_velocity(H,H,dpdx,mu),0.0,1e-14);
        EXPECT_NEAR(plane_poiseuille_velocity(0.5,H,dpdx,mu),0.25,1e-14);
        EXPECT_NEAR(plane_poiseuille_flow_rate_per_width(H,dpdx,mu),
                    1.0/6.0,1e-14);
    });

    run_case("M3_energy_balance_oracle", [] {
        const double T=fully_developed_temperature(
            2.0,300.0,100.0,1.0,1000.0,2.0,1.0,4.0);
        EXPECT_NEAR(T,300.4,1e-14);
    });

    return run_all();
}
