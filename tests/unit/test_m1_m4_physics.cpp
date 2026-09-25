#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/turbulence.h"
#include "cfdx/physics/thermal.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_advanced.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
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

        {
            // Explicit area-weighted reciprocity oracle: A1 F12 = A2 F21.
            const std::vector<double> areas{1.0,2.0};
            const std::vector<double> reciprocal_F{
                0.5,0.5,
                0.25,0.75};
            validate_view_factor_matrix(reciprocal_F,2,areas,1e-12);
            auto nonreciprocal_F=reciprocal_F;
            nonreciprocal_F[2]=0.30;
            EXPECT_THROW(validate_view_factor_matrix(nonreciprocal_F,2,areas,1e-12),
                         std::invalid_argument);
        }

        auto b=radiation_balance(100.0,100.0);
        EXPECT_NEAR(b.net,0.0,1e-14);
        EXPECT_NEAR(b.relative_error,0.0,1e-14);
    });

    run_case("s2s_bvh_matches_bruteforce", [] {
        // More than the leaf size so the regression exercises an actual BVH split.
        std::vector<RadiationTriangle> target;
        target.reserve(16);
        target.push_back({{0,0,1},{0,1,1},{1,0,1}});
        for (int i=1; i<16; ++i) {
            const double x=10.0 + static_cast<double>(i);
            target.push_back({{x,0,1},{x,1,1},{x+1,0,1}});
        }

        RadiationBvh bvh(target);
        const auto brute_force = [&target](
            const std::array<double,3>& origin,
            const std::array<double,3>& direction,
            double& distance) {
            distance=std::numeric_limits<double>::infinity();
            bool hit=false;
            for(const auto& tri:target) {
                double d=0.0;
                if(radiation_ray_triangle_hit(origin,direction,tri,d) &&
                   d<distance) {
                    distance=d;
                    hit=true;
                }
            }
            return hit;
        };
        const std::vector<std::pair<std::array<double,3>,std::array<double,3>>> rays{
            {{{0.1,0.1,0.0}},{{0.0,0.0,1.0}}},
            {{{0.1,0.1,2.0}},{{0.0,0.0,-1.0}}},
            {{{-0.5,0.1,0.0}},{{0.6,0.0,1.0}}},
            {{{0.1,0.1,0.0}},{{0.0,1.0,0.0}}}
        };
        for(const auto& ray:rays) {
            double direct=0.0, accelerated=0.0;
            const bool direct_hit=brute_force(ray.first,ray.second,direct);
            const bool accelerated_hit=bvh.nearest_hit(
                ray.first,ray.second,accelerated);
            EXPECT_TRUE(accelerated_hit==direct_hit);
            if(direct_hit) EXPECT_NEAR(accelerated,direct,1e-14);
            else EXPECT_TRUE(std::isinf(accelerated));
        }

        std::vector<RadiationTriangle> blockers{
            {{0,0,0.5},{0,1,0.5},{1,0,0.5}}};
        RadiationBvh blocker_bvh(blockers);
        double blocker_distance=0.0;
        EXPECT_TRUE(blocker_bvh.nearest_hit(
            {0.1,0.1,0.0},{0.0,0.0,1.0},blocker_distance));
        EXPECT_NEAR(blocker_distance,0.5,1e-14);

        RadiationBvh empty_bvh;
        double empty_distance=123.0;
        EXPECT_TRUE(!empty_bvh.nearest_hit(
            {0.1,0.1,0.0},{0.0,0.0,1.0},empty_distance));
        EXPECT_TRUE(std::isinf(empty_distance));
        EXPECT_THROW(
            bvh.nearest_hit({0.0,0.0,0.0},{0.0,0.0,0.0},empty_distance),
            std::invalid_argument);

        auto invalid=target;
        invalid[0].a[0]=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(bvh.build(invalid),std::invalid_argument);
        double preserved_distance=0.0;
        EXPECT_TRUE(bvh.nearest_hit(
            {0.1,0.1,0.0},{0.0,0.0,1.0},preserved_distance));
        EXPECT_NEAR(preserved_distance,1.0,1e-14);
    });

    return run_all();
}
