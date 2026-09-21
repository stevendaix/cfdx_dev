#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0; i<8; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("steady_incompressible_zero_state_is_fixed_point", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
        EXPECT_NEAR(U(0,0),0.0,1e-14);
        EXPECT_NEAR(U(0,1),0.0,1e-14);
        EXPECT_NEAR(U(0,2),0.0,1e-14);
        EXPECT_NEAR(p(0),0.0,1e-14);
        EXPECT_NEAR(r.history.back().continuity_linf,0.0,1e-14);
    });


    run_case("simplec_zero_state_uses_consistent_momentum_diagonal", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLEC;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
    });

    run_case("piso_zero_state_runs_multiple_pressure_corrections", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::PISO;
        c.coupling.n_pressure_correctors = 3;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
    });

    run_case("fixed_pressure_boundary_is_accepted_by_coupled_solver", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(100.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::FIXED_VALUE,100.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 100.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(p(0),100.0,1e-12);
    });

    return run_all();
}
