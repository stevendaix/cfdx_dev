// Regression: a full SIMPLE solve must carry the TVD memory across accepted
// iterations exactly once. Previously a duplicated commit emptied the saved
// face arrays and the second iteration threw an incompatible-size exception.
#include "cfdx/physics/steady_incompressible_solver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;

#define REQUIRE(condition) do { \
    if (!(condition)) throw std::runtime_error("REQUIRE failed: " #condition); \
} while (false)

static Mesh make_two_cell_mesh() {
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}};
    for (std::size_t i=0;i<12;++i)
        m.points().set(i,xyz[i][0],xyz[i][1],xyz[i][2]);
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
    for (std::size_t f=0;f<11;++f) {
        m.ownership().set_owner(f, (f==2 || f>=7) ? 1u : 0u);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1,0);
    m.ownership().set_neighbour(1,1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL;
    wall.face_ids={0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

int main() {
    const Mesh mesh=make_two_cell_mesh();
    for (bool recursive : {false,true}) {
        Field<double,Location::CELL> velocity(2,"U","m/s",3);
        Field<double,Location::CELL> pressure(2,"p","Pa",1);
        velocity.fill(0.0); pressure.fill(0.0);
        velocity.set(0,0.30,0.04,0.0);
        velocity.set(1,0.22,-0.03,0.0);
        VelocityBoundaryConditions velocity_bc;
        velocity_bc["wall"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ScalarBoundaryConditions pressure_bc;
        pressure_bc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        IncompressibleSolverControls controls;
        controls.algorithm=PressureVelocityAlgorithm::SIMPLE;
        controls.use_bounded_convection=true;
        controls.convection_scheme=ConvectionScheme::TVD;
        controls.tvd_deferred_memory.enabled=true;
        controls.tvd_deferred_memory.previous_weight=0.25;
        controls.tvd_deferred_memory.recursive=recursive;
        controls.convergence.max_iterations=4;
        controls.convergence.relative_tolerance=1e-14;
        controls.convergence.continuity_tolerance=1e-14;
        controls.pressure_reference_cell=0;
        controls.pressure_reference_value=0;
        controls.body_force={0.06,0.015,0.0};
        const auto result=solve_steady_incompressible(
            mesh,velocity,pressure,velocity_bc,pressure_bc,controls);
        REQUIRE(result.history.size()>=2);
        for (const auto& h:result.history) {
            REQUIRE(std::isfinite(h.tvd_memory_source_relative));
            REQUIRE(std::isfinite(h.continuity_normalized));
            REQUIRE(std::isfinite(h.momentum_equation_residual_rhs_relative));
        }
        for (std::size_t c=0;c<2;++c) {
            REQUIRE(std::isfinite(velocity(c,0)));
            REQUIRE(std::isfinite(velocity(c,1)));
            REQUIRE(std::isfinite(pressure(c)));
        }
        std::cout << "TVD_MEMORY_SOLVER_INTEGRATION: PASS recursive="
                  << recursive << " accepted=" << result.history.size() << '\n';
    }
}
