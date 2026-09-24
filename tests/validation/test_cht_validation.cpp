#include "cfdx/physics/cht_solver.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh cube(const char* patch_name)
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5});
    Patch patch;patch.name=patch_name;patch.type=PatchType::WALL;patch.face_ids={0,1,2,3,4,5};
    m.boundary().add_patch(patch);
    return m;
}

int main()
{
    run_case("two_region_cht_matches_interface_and_balances_flux",[] {
        Mesh m1=cube("interface1"),m2=cube("interface2");
        auto g1=build_fv_geometry(m1),g2=build_fv_geometry(m2);
        Field<double,Location::FACE> phi1(m1.n_faces(),"phi1","kg/s",1);
        Field<double,Location::FACE> phi2(m2.n_faces(),"phi2","kg/s",1);
        phi1.fill(0);phi2.fill(0);
        Field<double,Location::CELL> T1(1,"T1","K",1),T2(1,"T2","K",1);
        Field<double,Location::CELL> q1(1,"q1","W/m3",1),q2(1,"q2","W/m3",1);
        T1(0)=400;T2(0)=300;q1(0)=0;q2(0)=0;
        EnergySolverControls e1;e1.conductivity=1;e1.max_iterations=10;e1.tolerance=1e-10;e1.relaxation=1.0;
        EnergySolverControls e2=e1;
        ChtInterfaceControls c;c.region1_patch="interface1";c.region2_patch="interface2";
        c.conductivity1=1;c.conductivity2=1;c.matching_tolerance=1e-12;
        auto r=solve_two_region_cht(m1,g1,m2,g2,phi1,phi2,T1,T2,q1,q2,e1,e2,c);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(T1(0),350.0,1e-8);
        EXPECT_NEAR(T2(0),350.0,1e-8);
        EXPECT_NEAR(r.interface_imbalance,0.0,1e-10);
    });
    run_case("cht_face_override_is_used_by_independent_energy_balance",[] {
        Mesh m=cube("interface");
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1),oldT(1,"oldT","K",1),source(1,"S","W/m3",1);
        T(0)=300.0; oldT(0)=300.0; source(0)=0.0;
        EnergySolverControls e; e.conductivity=1.0;
        ScalarBoundaryFaceValues fv; fv.values["interface"].assign(m.n_faces(),400.0);
        const double b=energy_balance_relative(m,g,phi,T,oldT,source,e,{},&fv);
        EXPECT_NEAR(b,1.0,1e-14);
    });
    return run_all();
}
