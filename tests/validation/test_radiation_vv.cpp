#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

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

int main()
{
    run_case("radiation_p1_isothermal_cavity_exact_equilibrium", [] {
        Mesh m=one_cell_cube();
        auto g=build_fv_geometry(m);
        const double T=600.0;
        const double expected_G=4.0*blackbody_emissive_power(T);

        Field<double,Location::CELL> temperature(1,"T","K",1);
        Field<double,Location::CELL> G(1,"G","W/m2",1);
        Field<double,Location::CELL> qrad(1,"qrad","W/m3",1);
        temperature(0)=T; G(0)=0.0; qrad(0)=0.0;

        ScalarBoundaryConditions bc;
        bc["walls"]={ScalarBoundaryType::FIXED_VALUE,expected_G,0.0};

        const auto result=solve_p1_radiation(
            m,g,temperature,G,qrad,1.0,0.0,bc,100,1e-12);

        std::cout << "RADIATION_RESIDUAL: case=p1_isothermal_cavity"
                  << " iterations=" << result.iterations
                  << " residual=" << result.residual << '\n';

        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(G(0),expected_G,1e-10*expected_G);
        EXPECT_NEAR(qrad(0),0.0,1e-8);
    });

    run_case("radiation_gray_surface_exchange_analytic", [] {
        const double e1=0.7, e2=0.5, T1=800.0, T2=400.0;
        const double q=two_surface_net_exchange(e1,e2,T1,T2,1.0);
        const double oracle=STEFAN_BOLTZMANN*(std::pow(T1,4)-std::pow(T2,4)) /
            ((1.0-e1)/e1 + 1.0 + (1.0-e2)/e2);
        std::cout << "RADIATION_RESIDUAL: case=gray_surface_exchange"
                  << " analytic_error=" << std::abs(q-oracle)
                  << " flux=" << q << '\n';
        EXPECT_NEAR(q,oracle,1e-10*std::max(1.0,std::abs(oracle)));
    });

    run_case("radiation_diffuse_gray_reflection_identity", [] {
        const double eps=0.35, T=700.0, G=2500.0;
        const double I=gray_diffuse_wall_intensity(eps,T,G);
        const double oracle=eps*blackbody_intensity(T)+(1.0-eps)*G/M_PI;
        const double error=std::abs(I-oracle);
        std::cout << "RADIATION_RESIDUAL: case=diffuse_gray_reflection"
                  << " error=" << error << '\n';
        EXPECT_NEAR(I,oracle,1e-12*std::max(1.0,std::abs(oracle)));
    });

    return run_all();
}
