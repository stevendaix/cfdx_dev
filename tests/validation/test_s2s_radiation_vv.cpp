#include "cfdx/physics/radiation_s2s.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <limits>
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

        const auto q_into_fluid = s2s_heat_flux_into_fluid(r);
        EXPECT_NEAR(q_into_fluid[0], -oracle,
                    1e-10*std::max(1.0,std::abs(oracle)));
        EXPECT_NEAR(q_into_fluid[1], oracle,
                    1e-10*std::max(1.0,std::abs(oracle)));
    });

    run_case("s2s_ambient_irradiation_uses_open_view_fraction", [] {
        const double ambient = 250.0;
        const std::vector<double> areas{1.0, 1.0};
        const std::vector<double> eps{1.0, 1.0};
        const std::vector<double> temperatures{700.0, 400.0};
        S2SControls controls;
        controls.ambient_irradiation = ambient;

        // A closed enclosure has no view to ambient: ambient irradiation
        // must not be added on top of the complete surface-to-surface view.
        const std::vector<double> closed_F{0.0, 1.0, 1.0, 0.0};
        const auto closed = solve_s2s_radiosity(
            areas, eps, temperatures, closed_F, controls);
        EXPECT_NEAR(closed.irradiation[0], closed.radiosity[1],
                    1e-12 * std::max(1.0, closed.radiosity[1]));

        // With only half of each hemisphere covered by modeled surfaces,
        // the remaining half sees the ambient reservoir.
        const std::vector<double> open_F{0.0, 0.5, 0.5, 0.0};
        const auto open = solve_s2s_radiosity(
            areas, eps, temperatures, open_F, controls);
        const double expected = 0.5 * open.radiosity[1] + 0.5 * ambient;
        EXPECT_NEAR(open.irradiation[0], expected,
                    1e-12 * std::max(1.0, expected));

        // An explicit per-surface external irradiation vector keeps the
        // established API semantics: it replaces ambient irradiation.
        const std::vector<double> explicit_external{1000.0, 1000.0};
        const auto explicit_result = solve_s2s_radiosity(
            areas, eps, temperatures, closed_F, controls, explicit_external);
        EXPECT_NEAR(explicit_result.irradiation[0],
                    explicit_result.radiosity[1] + explicit_external[0],
                    1e-12 * std::max(1.0, explicit_result.irradiation[0]));
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

    run_case("s2s_centroid_geometry_rejects_nonfinite_inputs", [] {
        const ViewFactorPatch source{{0, 0, 0}, {0, 0, 1}, 2.0};
        const ViewFactorPatch target{{0, 0, 2}, {0, 0, -1}, 1.0};
        const double invalid[] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
        for (double value : invalid) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                auto bad = source;
                bad.center[axis] = value;
                EXPECT_THROW(patch_pair_view_factor(bad, target),
                             std::invalid_argument);
                EXPECT_THROW(patch_pair_view_factor(target, bad),
                             std::invalid_argument);
                EXPECT_THROW(estimate_view_factor_matrix({bad}),
                             std::invalid_argument);
                bad = source;
                bad.normal[axis] = value;
                EXPECT_THROW(patch_pair_view_factor(bad, target),
                             std::invalid_argument);
                EXPECT_THROW(patch_pair_view_factor(target, bad),
                             std::invalid_argument);
                EXPECT_THROW(estimate_view_factor_matrix({bad}),
                             std::invalid_argument);
            }
            auto bad = source;
            bad.area = value;
            EXPECT_THROW(patch_pair_view_factor(bad, target),
                         std::invalid_argument);
            EXPECT_THROW(patch_pair_view_factor(target, bad),
                         std::invalid_argument);
            EXPECT_THROW(estimate_view_factor_matrix({bad}),
                         std::invalid_argument);
        }
        auto bad = source;
        bad.normal = {0, 0, 0};
        EXPECT_THROW(estimate_view_factor_matrix({bad}),
                     std::invalid_argument);
        // Reject geometry before even the coincident-centroid early return.
        EXPECT_THROW(patch_pair_view_factor(bad, source),
                     std::invalid_argument);
        for (double area : {0.0, -1.0}) {
            bad = source;
            bad.area = area;
            EXPECT_THROW(patch_pair_view_factor(bad, source),
                         std::invalid_argument);
            EXPECT_THROW(estimate_view_factor_matrix({bad}),
                         std::invalid_argument);
        }
        auto distant = target;
        distant.center[2] = std::numeric_limits<double>::max();
        auto opposite = source;
        opposite.center[2] = -std::numeric_limits<double>::max();
        EXPECT_THROW(patch_pair_view_factor(opposite, distant),
                     std::overflow_error);
        // Finite geometry may still exceed the representable exchange range.
        // Reject it before a zero closure scale could multiply infinity.
        auto tiny = source;
        tiny.area = 1e-300;
        auto huge = target;
        huge.area = 1e308;
        huge.center[2] = 0.1;
        EXPECT_THROW(estimate_view_factor_matrix({tiny, huge}),
                     std::overflow_error);
    });

    run_case("s2s_centroid_geometry_preserves_reciprocity_and_limits", [] {
        const ViewFactorPatch source{{0, 0, 0}, {0, 0, 2}, 2.0};
        const ViewFactorPatch target{{0, 0, 2}, {0, 0, -3}, 1.0};
        const auto F = estimate_view_factor_matrix({source, target});
        const std::vector<double> areas{source.area, target.area};
        validate_s2s_view_factors(F, areas, 1e-12);
        EXPECT_NEAR(F[1], 1.0 / (4.0 * M_PI), 1e-14);
        EXPECT_NEAR(F[2], 2.0 / (4.0 * M_PI), 1e-14);
        EXPECT_NEAR(areas[0] * F[1], areas[1] * F[2], 1e-14);
        EXPECT_NEAR(F[0], 0.0, 0.0);
        EXPECT_NEAR(F[3], 0.0, 0.0);
        EXPECT_NEAR(patch_pair_view_factor(source, source), 0.0, 0.0);
        auto facing_away = target;
        facing_away.normal = {0, 0, 1};
        EXPECT_NEAR(patch_pair_view_factor(source, facing_away), 0.0, 0.0);

        // Scaling normals cannot change the projected solid angle.
        auto large_normal = source;
        large_normal.normal[2] = 1e200;
        EXPECT_NEAR(patch_pair_view_factor(large_normal, target), F[1], 1e-14);
        const auto enclosure = solve_s2s_radiosity(
            {1.0, 1.0}, {0.7, 0.5}, {800.0, 400.0},
            {0.0, 1.0, 1.0, 0.0});
        EXPECT_TRUE(enclosure.converged);
        EXPECT_NEAR(enclosure.total_power, 0.0, 1e-9);
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

        const auto matrix = estimate_s2s_view_factors_monte_carlo(
            {s1, s2}, {}, S2SMonteCarloControls{8192});
        const std::vector<double> areas{
            radiation_triangle_area(s1.front()),
            radiation_triangle_area(s2.front())};
        validate_s2s_view_factors(matrix, areas, 1e-12);
        EXPECT_NEAR(areas[0] * matrix[1], areas[1] * matrix[2], 1e-12);
    });

    run_case("s2s_triangle_kernel_rejects_nonfinite_geometry", [] {
        const RadiationTriangle triangle{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
        const double invalid[] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
        for (double value : invalid) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                for (std::size_t vertex = 0; vertex < 3; ++vertex) {
                    auto bad = triangle;
                    auto& coordinate = vertex == 0 ? bad.a : vertex == 1 ? bad.b : bad.c;
                    coordinate[axis] = value;
                    double distance = 0.0;
                    EXPECT_THROW(radiation_ray_triangle_hit(
                        {0.1, 0.1, 0}, {0, 0, 1}, bad, distance), std::invalid_argument);
                    EXPECT_THROW(radiation_triangle_area(bad), std::invalid_argument);
                    EXPECT_THROW(radiation_triangle_normal(bad), std::invalid_argument);
                    EXPECT_THROW(estimate_view_factor_ray_traced(
                        {bad}, {triangle}, {}, 1), std::invalid_argument);
                }
                std::array<double, 3> bad_ray{0, 0, 1};
                bad_ray[axis] = value;
                double distance = 0.0;
                EXPECT_THROW(radiation_ray_triangle_hit(
                    bad_ray, {0, 0, 1}, triangle, distance), std::invalid_argument);
                EXPECT_THROW(radiation_ray_triangle_hit(
                    {0.1, 0.1, 0}, bad_ray, triangle, distance), std::invalid_argument);
            }
        }
        double distance = 0.0;
        EXPECT_TRUE(radiation_ray_triangle_hit(
            {0.1, 0.1, 0}, {0, 0, 1}, triangle, distance));
        EXPECT_NEAR(distance, 1.0, 1e-14);
        EXPECT_NEAR(radiation_triangle_area(triangle), 0.5, 1e-14);
    });

    run_case("s2s_model_selector", [] {
        RadiationModelSelector selector;
        EXPECT_TRUE(selector.select_surface_to_surface(false)==RadiationApproximation::S2S);
        EXPECT_TRUE(selector.select_surface_to_surface(true)==RadiationApproximation::S2S_MONTE_CARLO);
        std::cout << "RADIATION_RESIDUAL: model_selector S2S=PASS S2S_MONTE_CARLO=PASS error=0\n";
    });

    return run_all();
}
