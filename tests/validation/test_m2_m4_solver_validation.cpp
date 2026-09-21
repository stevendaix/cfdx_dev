#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/radiation_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL; wall.face_ids={0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("kepsilon_transport_bounds_and_wall_closure",[] {
        Mesh m=make_unit_cube(); auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),e(1,"epsilon","m2/s3",1),S(1,"S","1/s",1);
        k(0)=1.0;e(0)=1.0;S(0)=1.0;
        TurbulenceTransportControls c;c.model=TurbulenceModel::KEPSILON;c.density=1.0;c.molecular_viscosity=1e-3;
        ScalarBoundaryConditions kb,eb; kb["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-3,0};
        eb["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-3,0};
        auto r=solve_kepsilon_transport(m,g,phi,k,e,S,c,kb,eb,10,1e-6);
        EXPECT_TRUE(r.iterations>0); EXPECT_TRUE(k(0)>=c.k_min); EXPECT_TRUE(e(0)>=c.epsilon_min);
        EXPECT_NEAR(wall_epsilon_from_k(1.0,0.1),std::pow(0.09,0.75)/0.1,1e-12);
    });

    run_case("sst_transport_preserves_positive_turbulence",[] {
        Mesh m=make_unit_cube(); auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),w(1,"omega","1/s",1),S(1,"S","1/s",1);
        Field<double,Location::CELL> F1(1,"F1","1",1),F2(1,"F2","1",1);
        k(0)=0.1;w(0)=1.0;S(0)=1.0;F1(0)=1.0;F2(0)=1.0;
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;c.molecular_viscosity=1e-5;
        ScalarBoundaryConditions bc;bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0};
        auto r=solve_sst_transport(m,g,phi,k,w,S,F1,F2,c,bc,bc,10,1e-6);
        EXPECT_TRUE(r.iterations>0);EXPECT_TRUE(k(0)>=c.k_min);EXPECT_TRUE(w(0)>=c.omega_min);
    });

    run_case("transient_energy_relaxes_to_constant_source_balance",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);phi.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1),old(1,"old","K",1),source(1,"Q","W/m3",1);
        T(0)=300;source(0)=1000;old(0)=300;
        ScalarBoundaryConditions bc;bc["wall"]={ScalarBoundaryType::FIXED_VALUE,300,0};
        EnergySolverControls c;c.density=1;c.cp=1000;c.conductivity=1;c.dt=1;c.max_iterations=10;c.tolerance=1e-10;c.relaxation=1.0;
        auto r=solve_energy(m,g,phi,T,source,c,bc);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(T(0),300.0 + 1000.0/(1000.0+12.0),1e-8);
    });

    run_case("dom_isotropic_blackbody_equilibrium",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::CELL> T(1,"T","K",1),G(1,"G","W/m2",1),q(1,"qrad","W/m3",1);
        T(0)=1000;G(0)=0;q(0)=0;
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);phi.fill(0);
        std::vector<DiscreteDirection> dirs={
            {1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},
            {0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},
            {0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};
        ScalarBoundaryConditions bc;bc["wall"]={ScalarBoundaryType::FIXED_VALUE,blackbody_intensity(1000),0};
        RadiationTransportControls c;c.absorption=1;c.scattering=0;c.max_iterations=100;c.tolerance=1e-10;
        auto r=solve_participating_radiation(m,g,T,G,q,dirs,c,bc);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(q(0),0.0,1e-5*blackbody_emissive_power(1000));
    });

    return run_all();
}
