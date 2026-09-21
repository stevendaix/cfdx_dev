#include "cfdx/core/numerics/convection.h"
#include "cfdx/physics/incompressible.h"
#include "cfdx/physics/pressure_velocity.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

Mesh make_cube() {
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

Mesh make_two_tetra() {
    Mesh m;
    const double p[5][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
    m.points().resize(5);
    for(std::size_t i=0;i<5;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,2,1}); m.faces().push_face({0,1,3}); m.faces().push_face({1,2,3});
    m.faces().push_face({2,0,3}); m.faces().push_face({1,4,2}); m.faces().push_face({1,3,4}); m.faces().push_face({3,2,4});
    m.ownership().resize(7);
    for(std::size_t f=0;f<7;++f){m.ownership().set_owner(f,f==2?0:(f<4?0:1));m.ownership().set_neighbour(f, f==2?1:FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3}); m.cells().push_cell({2,4,5,6});
    return m;
}

int main(){
    run_case("constant_scalar_convection_is_zero",[](){
        Mesh m=make_cube();
        Field<double,Location::CELL> phi(1,"phi","1",1); phi(0)=2.0;
        Field<double,Location::FACE> flux(6,"flux","m3/s",1);
        for(std::size_t f=0;f<6;++f) flux(f)=0.0;
        auto c=compute_convection(phi,flux,m,InterpScheme::UPWIND);
        EXPECT_NEAR(c(0),0.0,1e-12);
    });
    run_case("uniform_velocity_continuity_is_zero",[](){
        Mesh m=make_cube();
        Field<double,Location::CELL> u(1,"U","m/s",3);
        u.component_data(0)[0]=1.0; u.component_data(1)[0]=0.0; u.component_data(2)[0]=0.0;
        auto r=compute_continuity_residual(u,m);
        EXPECT_NEAR(r(0),0.0,1e-12);
    });
    run_case("momentum_terms_zero_for_constant_state",[](){
        Mesh m=make_cube();
        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        u.component_data(0)[0]=1.0; u.component_data(1)[0]=2.0; u.component_data(2)[0]=3.0; p(0)=5.0;
        auto terms=evaluate_momentum_terms(u,p,m,1e-3);
        for(std::size_t d=0;d<3;++d){
            EXPECT_NEAR(terms.convection.component_data(d)[0],0.0,1e-12);
            EXPECT_NEAR(terms.diffusion.component_data(d)[0],0.0,1e-12);
            EXPECT_NEAR(terms.pressure_gradient.component_data(d)[0],0.0,1e-12);
        }
    });
    run_case("pressure_correction_has_conservative_row_sum",[](){
        Mesh m=make_two_tetra();
        Field<double,Location::CELL> ap(2,"aP","1",1);
        Field<double,Location::CELL> r(2,"r","m3/s",1);
        ap(0)=ap(1)=1.0;
        const double v0 = 1.0/6.0, v1 = 1.0/3.0;
        r(0)=1.0; r(1)=-(v0/v1);
        auto pc=assemble_pressure_correction(m,ap,r);
        EXPECT_NEAR(pc.matrix(0,0)+pc.matrix(0,1),0.0,1e-12);
        EXPECT_NEAR(pc.matrix(1,0)+pc.matrix(1,1),0.0,1e-12);
        EXPECT_NEAR(pc.rhs(0)*v0+pc.rhs(1)*v1,0.0,1e-12);
        EXPECT_TRUE(pc.matrix.n_rows()==2);
    });
    return run_all();
}
