#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_single_cell_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        m.points().set(i, p[i][0], p[i][1], p[i][2]);

    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});

    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

static void make_zero_velocity_bcs(
    VelocityBoundaryConditions& ubc,
    ScalarBoundaryConditions& pbc)
{
    ubc["wall"] = {
        VelocityBoundaryCondition::Type::FIXED_VALUE,
        {0.0, 0.0, 0.0}};
    pbc["wall"] = {
        ScalarBoundaryType::ZERO_GRADIENT,
        0.0, 0.0};
}

static IncompressibleSolverControls base_controls(PressureVelocityAlgorithm algorithm)
{
    IncompressibleSolverControls c;
    c.algorithm = algorithm;
    c.convergence.max_iterations = 5;
    c.convergence.continuity_tolerance = 1e-12;
    c.linear_tolerance = 1e-12;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    return c;
}

int main()
{
    run_case("n9_algorithm_controls_cover_all_pressure_velocity_paths", [] {
        const std::array<PressureVelocityAlgorithm, 6> algorithms{
            PressureVelocityAlgorithm::SIMPLE,
            PressureVelocityAlgorithm::SIMPLEC,
            PressureVelocityAlgorithm::PISO,
            PressureVelocityAlgorithm::PIMPLE,
            PressureVelocityAlgorithm::FRACTIONAL_STEP,
            PressureVelocityAlgorithm::COUPLED};

        for (const auto algorithm : algorithms) {
            IncompressibleSolverControls c = base_controls(algorithm);
            validate_incompressible_controls(c, 1);
        }

        EXPECT_TRUE(!is_fractional_step_algorithm(PressureVelocityAlgorithm::PISO));
        EXPECT_TRUE(is_fractional_step_algorithm(PressureVelocityAlgorithm::FRACTIONAL_STEP));
        EXPECT_TRUE(!is_coupled_algorithm(PressureVelocityAlgorithm::SIMPLE));
        EXPECT_TRUE(is_coupled_algorithm(PressureVelocityAlgorithm::COUPLED));
    });

    run_case("n9_all_segregated_algorithms_preserve_zero_state", [] {
        const Mesh m = make_single_cell_cube();
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        make_zero_velocity_bcs(ubc, pbc);

        const std::array<PressureVelocityAlgorithm, 5> algorithms{
            PressureVelocityAlgorithm::SIMPLE,
            PressureVelocityAlgorithm::SIMPLEC,
            PressureVelocityAlgorithm::PISO,
            PressureVelocityAlgorithm::PIMPLE,
            PressureVelocityAlgorithm::FRACTIONAL_STEP};

        for (const auto algorithm : algorithms) {
            Field<double, Location::CELL> U(1, "U", "m/s", 3);
            Field<double, Location::CELL> p(1, "p", "Pa", 1);
            U.fill(0.0);
            p.fill(0.0);

            auto c = base_controls(algorithm);
            c.coupling.n_outer_correctors = 2;
            c.coupling.n_pressure_correctors = 2;
            c.coupling.n_fractional_steps = 2;

            const auto r = solve_steady_incompressible(m, U, p, ubc, pbc, c);
            EXPECT_TRUE(r.converged);
            EXPECT_TRUE(!r.history.empty());
            EXPECT_NEAR(U(0,0), 0.0, 1e-12);
            EXPECT_NEAR(U(0,1), 0.0, 1e-12);
            EXPECT_NEAR(U(0,2), 0.0, 1e-12);
            EXPECT_NEAR(p(0), 0.0, 1e-12);
            EXPECT_NEAR(r.history.back().mass_local_linf, 0.0, 1e-12);
            EXPECT_NEAR(r.history.back().corrected_flux_continuity_linf, 0.0, 1e-12);
        }
    });

    run_case("n9_pressure_reference_removes_only_constant_pressure_mode", [] {
        const Mesh m = make_single_cell_cube();
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        make_zero_velocity_bcs(ubc, pbc);

        Field<double, Location::CELL> U0(1, "U", "m/s", 3);
        Field<double, Location::CELL> p0(1, "p", "Pa", 1);
        Field<double, Location::CELL> U1(1, "U", "m/s", 3);
        Field<double, Location::CELL> p1(1, "p", "Pa", 1);
        U0.fill(0.0);
        U1.fill(0.0);
        p0(0) = 0.0;
        p1(0) = 137.25;

        auto c = base_controls(PressureVelocityAlgorithm::SIMPLE);
        const auto r0 = solve_steady_incompressible(m, U0, p0, ubc, pbc, c);
        const auto r1 = solve_steady_incompressible(m, U1, p1, ubc, pbc, c);

        EXPECT_TRUE(r0.converged);
        EXPECT_TRUE(r1.converged);
        EXPECT_NEAR(U0(0,0), U1(0,0), 1e-12);
        EXPECT_NEAR(U0(0,1), U1(0,1), 1e-12);
        EXPECT_NEAR(U0(0,2), U1(0,2), 1e-12);
        EXPECT_NEAR(p0(0), c.pressure_reference_value, 1e-12);
        EXPECT_NEAR(p1(0), c.pressure_reference_value, 1e-12);
    });

    run_case("n9_rhie_chow_constant_pressure_invariance", [] {
        const Mesh m = make_single_cell_cube();
        const auto geometry = build_fv_geometry(m);

        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p0(1, "p", "Pa", 1);
        Field<double, Location::CELL> p1(1, "p", "Pa", 1);
        U.fill(0.0);
        p0(0) = 0.0;
        p1(0) = 1.0e8;

        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        make_zero_velocity_bcs(ubc, pbc);

        const std::vector<double> rAU{0.25};
        const auto phi0 = make_rhie_chow_mass_flux(
            m, geometry, U, p0, rAU, 1.0, ubc);
        const auto phi1 = make_rhie_chow_mass_flux(
            m, geometry, U, p1, rAU, 1.0, ubc);

        EXPECT_TRUE(phi0.size() == phi1.size());
        for (std::size_t f = 0; f < phi0.size(); ++f) {
            EXPECT_TRUE(std::isfinite(phi0(f)));
            EXPECT_TRUE(std::isfinite(phi1(f)));
            EXPECT_NEAR(phi0(f), phi1(f), 1e-12);
            EXPECT_NEAR(phi0(f), 0.0, 1e-12);
        }
    });

    run_case("n9_rhie_chow_rejects_invalid_inverse_momentum_diagonal", [] {
        const Mesh m = make_single_cell_cube();
        const auto geometry = build_fv_geometry(m);
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        make_zero_velocity_bcs(ubc, pbc);

        const std::array<std::vector<double>, 3> nan_rAU{
            std::vector<double>{std::numeric_limits<double>::quiet_NaN()},
            std::vector<double>{0.25},
            std::vector<double>{0.25}};
        EXPECT_THROW(
            make_rhie_chow_mass_flux(
                m, geometry, U, p, nan_rAU,
                1.0, ubc, pbc),
            std::invalid_argument);

        const std::array<std::vector<double>, 3> zero_rAU{
            std::vector<double>{0.0},
            std::vector<double>{0.25},
            std::vector<double>{0.25}};
        EXPECT_THROW(
            make_rhie_chow_mass_flux(
                m, geometry, U, p, zero_rAU,
                1.0, ubc, pbc),
            std::invalid_argument);
    });

    run_case("n9_piso_corrector_count_is_executed_not_decorative", [] {
        const Mesh m = make_single_cell_cube();
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        make_zero_velocity_bcs(ubc, pbc);

        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        auto c = base_controls(PressureVelocityAlgorithm::PISO);
        c.coupling.n_pressure_correctors = 3;

        const auto r = solve_steady_incompressible(m, U, p, ubc, pbc, c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
        EXPECT_TRUE(r.pressure_linear_context.solves ==
                    3 * r.iterations);
    });

    run_case("n9_invalid_coupling_parameters_are_rejected", [] {
        EXPECT_THROW(
            validate_coupling_controls(
                CouplingControls{0.7, 0.3, 0, 2, 2, 100, 1e-10}),
            std::invalid_argument);
        EXPECT_THROW(
            validate_coupling_controls(
                CouplingControls{0.7, 0.3, 1, 0, 2, 100, 1e-10}),
            std::invalid_argument);
        EXPECT_THROW(
            validate_coupling_controls(
                CouplingControls{0.7, 0.3, 1, 2, 0, 100, 1e-10}),
            std::invalid_argument);
        EXPECT_THROW(
            validate_coupling_controls(
                CouplingControls{0.0, 0.3, 1, 2, 2, 100, 1e-10}),
            std::invalid_argument);
    });

    return run_all();
}
