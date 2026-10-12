#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});
    m.faces().push_face({1,2,6,5});
    m.faces().push_face({8,9,11,10});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static void configure_case(IncompressibleSolverControls& c, PressureGaugePolicy policy)
{
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.pressure_gauge_policy = policy;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 3.0;
    c.convergence.max_iterations = 2;
    c.convergence.continuity_tolerance = 1e-12;
    c.convergence.relative_tolerance = 1e-12;
    c.linear_max_iterations = 100;
    c.linear_tolerance = 1e-12;
}

int main()
{
    run_case("n9_s10_reference_cell_policy", [] {
        const Mesh mesh = make_two_cell_channel();
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(7.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        IncompressibleSolverControls c;
        configure_case(c, PressureGaugePolicy::REFERENCE_CELL);
        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(p(0), 3.0, 1e-12);
        EXPECT_NEAR(p(1), 3.0, 1e-12);
    });

    run_case("n9_s10_zero_mean_policy", [] {
        const Mesh mesh = make_two_cell_channel();
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 2.0;
        p(1) = 6.0;
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        IncompressibleSolverControls c;
        configure_case(c, PressureGaugePolicy::ZERO_MEAN);
        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(0.5 * (p(0) + p(1)), 0.0, 1e-12);
        EXPECT_NEAR(p(1) - p(0), 4.0, 1e-12);
    });

    run_case("n9_s10_none_policy_preserves_pressure_level", [] {
        const Mesh mesh = make_two_cell_channel();
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(7.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        IncompressibleSolverControls c;
        configure_case(c, PressureGaugePolicy::NONE);
        c.pressure_linear_solver.null_space = NullSpaceModel::Constant;
        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(p(0), 7.0, 1e-12);
        EXPECT_NEAR(p(1), 7.0, 1e-12);
    });
    return 0;
}
