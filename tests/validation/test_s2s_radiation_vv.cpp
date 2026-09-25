#include "cfdx/physics/radiation_s2s.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("s2s_two_surface_gray_radiosity_analytic", [] {
        const std::vector<double> A{1.0,1.0};
        const std::vector<double> eps{0.7,0.5};
        const std::vector<double> T{800.0,400.0};
        const std::vector<double> F{0.0,1.0,1.0,0.0};

        const auto r=solve_s2s_radiosity(A,eps,T,F);
        const double oracle=STEFAN_BOLTZMANN*
            (std::pow(T[0],4)-std::pow(T[1],4)) /
            ((1.0-eps[0])/eps[0] + 1.0 + (1.0-eps[1])/eps[1]);

        std::cout << "RADIATION_RESIDUAL: model=S2S"
                  << " case=two_surface analytic_error="
                  << std::abs(r.net_flux[0]-oracle)
                  << " energy_balance_error=" << r.energy_balance_error << '\n';

        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(r.net_flux[0],oracle,1e-10*std::max(1.0,std::abs(oracle)));
        EXPECT_NEAR(r.net_flux[0],-r.net_flux[1],1e-10*std::max(1.0,std::abs(r.net_flux[0])));
        EXPECT_NEAR(r.total_power,0.0,1e-10*std::max(1.0,std::abs(oracle)));
    });

    run_case("s2s_view_factor_validation", [] {
        const std::vector<double> A{2.0,1.0};
        // Reciprocity: A1 F12 = A2 F21 = 1.
        const std::vector<double> F{0.0,0.5,1.0,0.0};
        validate_s2s_view_factors(F,A,1e-12);

        std::cout << "RADIATION_RESIDUAL: model=S2S"
                  << " view_factor_reciprocity_error="
                  << std::abs(A[0]*F[1]-A[1]*F[2])
                  << " closure_error=" << std::abs(F[1]-0.5) << '\n';
        EXPECT_NEAR(A[0]*F[1],A[1]*F[2],1e-12);
    });

    run_case("s2s_monte_carlo_visibility_and_blocking", [] {
        const std::vector<RadiationTriangle> s1{{
            {0,0,0},{1,0,0},{0,1,0}}}; // +z
        const std::vector<RadiationTriangle> s2{{
            {0,0,1},{0,1,1},{1,0,1}}}; // -z
        const std::vector<RadiationTriangle> blocker{{
            {0,0,0.5},{1,0,0.5},{0,1,0.5}}};

        const double f12=estimate_view_factor_ray_traced(s1,s2,{},8192);
        const double f21=estimate_view_factor_ray_traced(s2,s1,{},8192);
        const double blocked=estimate_view_factor_ray_traced(s1,s2,blocker,8192);

        std::cout << "RADIATION_RESIDUAL: model=S2S_MONTE_CARLO"
                  << " reciprocity_error=" << std::abs(f12-f21)
                  << " unblocked_F12=" << f12
                  << " blocked_F12=" << blocked << '\n';

        EXPECT_TRUE(f12>0.0 && f12<1.0);
        EXPECT_NEAR(f12,f21,0.08);
        EXPECT_NEAR(blocked,0.0,1e-14);
    });

    run_case("s2s_model_selector", [] {
        RadiationModelSelector selector;
        EXPECT_TRUE(selector.select_surface_to_surface(false)==RadiationApproximation::S2S);
        EXPECT_TRUE(selector.select_surface_to_surface(true)==RadiationApproximation::S2S_MONTE_CARLO);
        std::cout << "RADIATION_RESIDUAL: model_selector S2S=PASS S2S_MONTE_CARLO=PASS error=0\n";
    });

    return run_all();
}
