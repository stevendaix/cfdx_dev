#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/turbulence.h"
#include "cfdx/physics/thermal.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_advanced.h"
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
                    1.0, 1e-12);
        EXPECT_TRUE(piso_correction_gain(2.0, 0.5) > 0.0);
    });

    run_case("rans_and_les_models_are_positive", [] {
        EXPECT_NEAR(turbulent_kinematic_viscosity_kepsilon(4.0, 2.0), 0.72, 1e-12);
        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(1.0, 2.0, 0.5) > 0.0);
        EXPECT_NEAR(smagorinsky_eddy_viscosity(0.1, 20.0),
                    0.00578, 1e-12);
        EXPECT_TRUE(des_eddy_viscosity(0.1, 0.02, 20.0) >= 0.0);
    });

    run_case("thermal_and_cht_fluxes", [] {
        EXPECT_NEAR(heat_flux_conduction(10.0, 300.0, 290.0, 0.5),
                    200.0, 1e-12);
        const double g = cht_interface_conductance(10.0, 2.0, 0.01, 0.02, 1.0);
        EXPECT_NEAR(g, 90.9090909090909, 1e-10);
        EXPECT_NEAR(cht_interface_heat_flux(g, 320.0, 300.0),
                    1818.181818181818, 1e-9);
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
        validate_view_factor_matrix(F, 2, {2.0, 2.0});

        const double q_equal=two_surface_net_exchange(0.5,0.5,900,500,1.0);
        const double q_area=two_surface_net_exchange(0.5,0.5,900,500,1.0,2.0,4.0);
        EXPECT_TRUE(q_area > q_equal);
        EXPECT_NEAR(gray_diffuse_wall_intensity(1.0,300.0,0.0),
                    blackbody_intensity(300.0),1e-12);

        validate_discrete_directions({
            {1.0, 0.0, 0.0, 2.0 * M_PI / 3.0},
            {-1.0, 0.0, 0.0, 2.0 * M_PI / 3.0},
            {0.0, 1.0, 0.0, 2.0 * M_PI / 3.0},
            {0.0, -1.0, 0.0, 2.0 * M_PI / 3.0},
            {0.0, 0.0, 1.0, 2.0 * M_PI / 3.0},
            {0.0, 0.0, -1.0, 2.0 * M_PI / 3.0}
        });
    });

    run_case("radiation_advanced_physical_invariants", [] {
        RadiationOpticalProperties p;
        p.absorption=2.0; p.scattering=1.0; p.emissivity=0.8;
        p.validate();
        EXPECT_NEAR(p.extinction(),3.0,1e-14);
        EXPECT_NEAR(p.optical_thickness(2.0),6.0,1e-14);

        std::vector<RadiationBand> bands(2);
        bands[0].wavelength_min=1e-6; bands[0].wavelength_max=2e-6;
        bands[0].weight=0.4; bands[0].properties=p;
        bands[1].wavelength_min=2e-6; bands[1].wavelength_max=3e-6;
        bands[1].weight=0.6; bands[1].properties=p;
        EXPECT_NEAR(weighted_band_absorption(bands,1000.0),2.0,1e-14);

        std::vector<DiscreteDirection> dirs{
            {1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},
            {0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},
            {0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};
        std::vector<double> incoming(6,blackbody_intensity(300.0));
        std::vector<double> wall(6,0.0);
        apply_diffuse_gray_wall(dirs,incoming,{1,0,0},1.0,300.0,wall);
        EXPECT_TRUE(wall[0] > 0.0);
        EXPECT_NEAR(wall[2],0.0,1e-14);

        std::vector<ViewFactorPatch> patches{
            {{0,0,0},{0,0,1},1.0},
            {{0,0,1},{0,0,-1},1.0}};
        auto F=estimate_view_factor_matrix(patches);
        EXPECT_TRUE(F[1] >= 0.0 && F[1] <= 1.0);
        EXPECT_TRUE(F[0] >= 0.0 && F[0] <= 1.0);
        EXPECT_TRUE(F[1] >= 0.0 && F[1] <= 1.0);

        std::vector<RadiationTriangle> s1{{
            {0,0,0},{1,0,0},{0,1,0}}};
        std::vector<RadiationTriangle> s2{{
            {0,0,1},{0,1,1},{1,0,1}}};
        const double F12=estimate_view_factor_ray_traced(s1,s2,{},2000);
        const double F21=estimate_view_factor_ray_traced(s2,s1,{},2000);
        EXPECT_TRUE(F12>0.0 && F12<1.0);
        EXPECT_NEAR(F12,F21,0.05);

        auto b=radiation_balance(100.0,100.0);
        EXPECT_NEAR(b.net,0.0,1e-14);
        EXPECT_NEAR(b.relative_error,0.0,1e-14);
    });

    return run_all();
}
