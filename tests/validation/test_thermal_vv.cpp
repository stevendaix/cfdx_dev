#include "cfdx/physics/energy_solver.h"
#include "common/test_harness.h"

#include <iostream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh two_cell_1d()
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{0,1,0},{0,1,1},{0,0,1},
        {0.5,0,0},{0.5,1,0},{0.5,1,1},{0.5,0,1},
        {1,0,0},{1,1,0},{1,1,1},{1,0,1}
    };
    for (std::size_t i=0; i<12; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.faces().push_face({8,9,10,11});
    m.faces().push_face({4,8,9,5});
    m.faces().push_face({7,6,10,11});
    m.faces().push_face({4,7,11,8});
    m.faces().push_face({5,9,10,6});

    m.ownership().resize(11);
    m.ownership().set_owner(0,0); m.ownership().set_neighbour(0,FaceOwnership::BOUNDARY);
    m.ownership().set_owner(1,0); m.ownership().set_neighbour(1,1);
    for (std::size_t f=2; f<=5; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(6,1); m.ownership().set_neighbour(6,FaceOwnership::BOUNDARY);
    for (std::size_t f=7; f<=10; ++f) {
        m.ownership().set_owner(f,1);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }

    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({1,6,7,8,9,10});

    Patch left; left.name="left"; left.type=PatchType::WALL; left.face_ids={0};
    Patch right; right.name="right"; right.type=PatchType::WALL; right.face_ids={6};
    Patch walls; walls.name="walls"; walls.type=PatchType::WALL;
    walls.face_ids={2,3,4,5,7,8,9,10};
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(walls);
    return m;
}

static Mesh one_cell_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for (std::size_t i=0;i<8;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0;f<6;++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall; wall.name="walls"; wall.type=PatchType::WALL; wall.face_ids={0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

static void print_energy_history(const std::string& name, const EnergySolveResult& result)
{
    std::cout << "THERMAL_VV: " << name
              << " converged=" << (result.converged ? "YES" : "NO")
              << " iterations=" << result.iterations << '\n';
    for (const auto& h : result.history) {
        std::cout << "THERMAL_RESIDUAL: case=" << name
                  << " iter=" << h.iteration
                  << " residual=" << h.residual
                  << " energy_imbalance=" << h.energy_imbalance << '\n';
    }
}

int main()
{
    run_case("thermal_1d_conduction_linear_profile", [] {
        Mesh m=two_cell_1d();
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> T(2,"T","K",1), source(2,"source","W/m3",1);
        T(0)=350.0; T(1)=350.0; source.fill(0.0);

        ScalarBoundaryConditions bc;
        bc["left"]={ScalarBoundaryType::FIXED_VALUE,400.0,0.0};
        bc["right"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        bc["walls"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        EnergySolverControls c;
        c.conductivity=1.0; c.relaxation=1.0; c.max_iterations=20; c.tolerance=1e-12;
        const auto result=solve_energy(m,g,phi,T,source,c,bc);
        print_energy_history("1d_conduction",result);

        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),375.0,1e-10);
        EXPECT_NEAR(T(1),325.0,1e-10);
        EXPECT_TRUE(result.history.back().energy_imbalance < 1e-12);
        std::cout << "THERMAL_RESIDUAL: source=volumetric_generation PASS"
                  << " value=" << source(0)
                  << " final_T=" << T(0)
                  << " energy_imbalance=" << result.history.back().energy_imbalance << '\n';
    });

    run_case("thermal_uniform_generation_conduction_balance", [] {
        Mesh m=one_cell_cube();
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=100.0;

        ScalarBoundaryFaceValues fv;
        fv.values["walls"].assign(m.n_faces(),310.0);

        EnergySolverControls c;
        c.conductivity=1.0; c.relaxation=1.0; c.max_iterations=20; c.tolerance=1e-12;
        const auto result=solve_energy(m,g,phi,T,source,c,{},&fv);
        print_energy_history("uniform_generation",result);

        // Six unit-area faces, each at d=0.5 from the cell centre:
        // 6*k*(T_cell-310)/0.5 = q'''*V, hence T_cell=310+100/12.
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),310.0 + 100.0/12.0,1e-10);
        EXPECT_TRUE(result.history.back().energy_imbalance < 1e-12);
    });

    run_case("thermal_convection_transports_enthalpy", [] {
        Mesh m=one_cell_cube();
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);
        phi.fill(0.0);
        phi(0)=-2.0;
        phi(1)=2.0;
        Field<double,Location::CELL> source(1,"source","W/m3",1);
        Field<double,Location::CELL> old_temperature(1,"T_old","K",1);
        source.fill(0.0);
        old_temperature(0)=300.0;

        ScalarBoundaryConditions bc;
        bc["walls"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        EnergySolverControls c;
        c.cp=1250.0;
        c.conductivity=0.0;

        const auto equation=assemble_energy_equation(
            m,g,phi,source,old_temperature,c,bc);
        EXPECT_NEAR(equation.matrix(0,0),2.0*c.cp,1e-12);
        EXPECT_NEAR(equation.rhs(0),2.0*c.cp*300.0,1e-9);
    });

    return run_all();
}
