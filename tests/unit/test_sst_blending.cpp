#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/sst_strain_rate.h"
#include <array>
#include "common/test_harness.h"
#include <cmath>
#include <algorithm>
#include <limits>

using namespace cfdx::testing;
using namespace cfdx::physics;

int main()
{
    run_case("sst_blending_is_bounded", [] {
        const auto [f1, f2] = compute_sst_blending(0.5, 10.0, 0.01, 1.5e-5, 0.09);
        EXPECT_TRUE(f1 >= 0.0 && f1 <= 1.0);
        EXPECT_TRUE(f2 >= 0.0 && f2 <= 1.0);
    });

    run_case("sst_blending_responds_to_wall_distance", [] {
        const auto near_wall = compute_sst_blending(0.5, 10.0, 1.0e-3, 1.5e-5, 0.09);
        const auto far_wall = compute_sst_blending(0.5, 10.0, 10.0, 1.5e-5, 0.09);
        EXPECT_TRUE(near_wall.first > far_wall.first);
        EXPECT_TRUE(near_wall.second > far_wall.second);
    });

    run_case("sst_production_viscosity_uses_cellwise_sst_closure", [] {
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(2, "k", "m2/s2", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> omega(2, "omega", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> strain(2, "strain", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> f2(2, "F2", "1", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> nu_t(2, "nu_t", "m2/s", 1);
        k(0) = 0.5; k(1) = 1.0;
        omega(0) = 10.0; omega(1) = 20.0;
        strain(0) = 2.0; strain(1) = 8.0;
        f2(0) = 0.5; f2(1) = 1.0;
        nu_t.fill(7.0);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.a1 = 0.31;

        refresh_sst_production_eddy_viscosity(k, omega, strain, f2, controls, nu_t);
        EXPECT_NEAR(nu_t(0), 0.31 * 0.5 / std::max(0.31 * 10.0, 2.0 * 0.5), 1e-14);
        EXPECT_NEAR(nu_t(1), 0.31 * 1.0 / std::max(0.31 * 20.0, 8.0), 1e-14);
        EXPECT_TRUE(nu_t(0) != nu_t(1));
    });

    run_case("sst_production_viscosity_rejects_invalid_state_atomically", [] {
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(2), omega(2), strain(2), f2(2), nu_t(2);
        k(0) = 0.5; k(1) = 1.0;
        omega(0) = 10.0; omega(1) = 20.0;
        strain(0) = 2.0; strain(1) = 4.0;
        f2(0) = 0.5; f2(1) = std::numeric_limits<double>::quiet_NaN();
        nu_t(0) = 3.0; nu_t(1) = 4.0;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(refresh_sst_production_eddy_viscosity(k, omega, strain, f2, controls, nu_t),
                     std::invalid_argument);
        EXPECT_NEAR(nu_t(0), 3.0, 0.0);
        EXPECT_NEAR(nu_t(1), 4.0, 0.0);
    });

    run_case("sst_transaction_preflight_rejects_field_mismatch_without_mutation", [] {
        cfdx::core::Mesh mesh;
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(1), omega(1), strain(1), wall(1), nu_t(1);
        cfdx::core::Field<double, cfdx::core::Location::FACE> flux(0);
        k(0) = 0.4; omega(0) = 8.0; strain(0) = 1.0; wall(0) = 0.01; nu_t(0) = 0.02;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, FvGeometry{}, flux, k, omega, strain, wall, controls, nu_t),
            std::invalid_argument);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });


    run_case("sst_strain_rate_magnitude_matches_simple_shear", [] {
        // Ux = gamma*y, Uy = Uz = 0 gives |S| = |gamma|.
        const double gamma = 3.25;
        const std::array<double, 9> grad{
            0.0, gamma, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0};
        EXPECT_NEAR(sst_strain_rate_magnitude(grad), gamma, 1e-14);
    });

    run_case("sst_strain_rate_magnitude_matches_extension", [] {
        // Ux = a*x and Uy = -a*y gives |S| = 2*|a|.
        const double a = 1.75;
        const std::array<double, 9> grad{
            a, 0.0, 0.0,
            0.0, -a, 0.0,
            0.0, 0.0, 0.0};
        EXPECT_NEAR(sst_strain_rate_magnitude(grad), 2.0 * a, 1e-14);
    });

    run_case("sst_strain_rate_magnitude_rejects_nonfinite_gradient", [] {
        const std::array<double, 9> grad{
            0.0, 0.0, 0.0,
            0.0, std::numeric_limits<double>::infinity(), 0.0,
            0.0, 0.0, 0.0};
        EXPECT_THROW(sst_strain_rate_magnitude(grad), std::invalid_argument);
    });

    run_case("sst_blending_rejects_invalid_distance", [] {
        EXPECT_THROW(compute_sst_blending(0.5, 10.0, 0.0, 1.5e-5, 0.09), std::invalid_argument);
    });

    return run_all();
}
