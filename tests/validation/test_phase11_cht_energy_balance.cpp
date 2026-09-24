#include "cfdx/physics/energy_solver.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::physics;
using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh make_unit_cube() {
    Mesh m; m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall; wall.name="interface"; wall.type=PatchType::WALL; wall.face_ids={0,1,2,3,4,5}; m.boundary().add_patch(wall);
    return m;
}

int main() {
    run_case("M3_face_override_enters_independent_energy_balance", [] {
        Mesh m=make_unit_cube();
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1), Told(1,"Told","K",1), source(1,"S","W/m3",1);
        T(0)=300.0; Told(0)=300.0; source(0)=0.0;
        EnergySolverControls c; c.conductivity=1.0;
        ScalarBoundaryFaceValues fv;
        fv.values["interface"].assign(m.n_faces(),400.0);
        const double balance=energy_balance_relative(m,g,flux,T,Told,source,c,{},&fv);
        EXPECT_TRUE(balance>0.0);
        EXPECT_NEAR(balance,1.0,1e-14);
    });
    return 0;
}
