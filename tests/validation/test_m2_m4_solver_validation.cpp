#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/radiation_solver.h"
#include "cfdx/physics/radiation_advanced.h"
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

    run_case("rosseland_isothermal_equilibrium",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);phi.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1),S(1,"S","W/m3",1),a(1,"a","1/m",1);
        T(0)=1000.0;S(0)=0.0;a(0)=1.0;
        ScalarBoundaryConditions bc;bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1000.0,0.0};
        EnergySolverControls ec;ec.density=1.0;ec.cp=1000.0;ec.conductivity=0.0;
        ec.max_iterations=20;ec.tolerance=1e-10;ec.relaxation=1.0;
        RosselandSolveControls rc;rc.max_iterations=20;rc.tolerance=1e-10;rc.relaxation=1.0;
        auto r=solve_rosseland_energy(m,g,phi,T,S,a,ec,rc,bc);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(T(0),1000.0,1e-8);
    });

    run_case("advanced_radiation_variable_properties_and_spectral",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::CELL> T(1,"T","K",1),G(1,"G","W/m2",1),Q(1,"Q","W/m3",1);
        T(0)=800.0;G(0)=0.0;Q(0)=0.0;
        std::vector<DiscreteDirection> dirs{
            {1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},
            {0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},
            {0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};
        RadiationOpticalPropertyField props;
        props.cells.resize(1);
        props[0].absorption=0.5;props[0].scattering=0.1;props[0].emissivity=1.0;
        RadiationTransportControls rc;rc.max_iterations=30;rc.linear_max_iterations=200;rc.tolerance=1e-7;
        auto vr=solve_participating_radiation_variable_properties(
            m,g,T,G,Q,dirs,props,rc);
        EXPECT_TRUE(vr.converged);
        EXPECT_TRUE(G(0)>=0.0);

        std::vector<RadiationBand> bands(2);
        bands[0].wavelength_min=1e-6;bands[0].wavelength_max=2e-6;bands[0].weight=0.5;
        bands[0].properties=props[0];
        bands[1].wavelength_min=2e-6;bands[1].wavelength_max=4e-6;bands[1].weight=0.5;
        bands[1].properties=props[0];
        G(0)=0.0;Q(0)=0.0;
        auto sr=solve_spectral_dom(m,g,T,G,Q,dirs,bands,rc);
        EXPECT_TRUE(sr.converged);
        EXPECT_TRUE(std::isfinite(G(0)) && G(0)>=0.0);
        EXPECT_TRUE(std::isfinite(Q(0)));
    });

    run_case("dom_diffuse_gray_wall_boundary_is_energy_consistent",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::CELL> T(1,"T","K",1),G(1,"G","W/m2",1),q(1,"qrad","W/m3",1);
        T(0)=1000.0;G(0)=0.0;q(0)=0.0;
        const double wa=M_PI*(0.5-std::sqrt(3.0)/6.0);
        const double wc=M_PI*(1.0+std::sqrt(3.0))/8.0;
        const double mu=1.0/std::sqrt(3.0);
        std::vector<DiscreteDirection> dirs{
            {1,0,0,wa},{-1,0,0,wa},
            {0,1,0,wa},{0,-1,0,wa},
            {0,0,1,wa},{0,0,-1,wa},
            {mu,mu,mu,wc},{mu,mu,-mu,wc},{mu,-mu,mu,wc},{mu,-mu,-mu,wc},
            {-mu,mu,mu,wc},{-mu,mu,-mu,wc},{-mu,-mu,mu,wc},{-mu,-mu,-mu,wc}};
        DomWallBoundaryConditions walls;
        for (std::size_t p=0;p<m.boundary().n_patches();++p)
            walls[m.boundary().patch(p).name]={0.5,1000.0};
        RadiationTransportControls rc;
        rc.absorption=1.0;rc.scattering=0.0;rc.max_iterations=100;
        rc.tolerance=1e-9;rc.linear_tolerance=1e-11;
        const auto r=solve_participating_radiation_diffuse_gray_walls(
            m,g,T,G,q,dirs,rc,walls);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(q(0),0.0,1e-6*blackbody_emissive_power(1000.0));
        EXPECT_NEAR(G(0),4.0*blackbody_emissive_power(1000.0),
                    1e-5*blackbody_emissive_power(1000.0));
    });

    run_case("rosseland_scattering_uses_transport_opacity",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);phi.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1),S(1,"S","W/m3",1);
        Field<double,Location::CELL> a(1,"a","1/m",1),s(1,"s","1/m",1);
        T(0)=1000.0;S(0)=0.0;a(0)=2.0;s(0)=3.0;
        EXPECT_NEAR(rosseland_conductivity(1000.0,2.0,3.0),
                    rosseland_conductivity(1000.0,5.0),1e-12);
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1000.0,0.0};
        EnergySolverControls ec;ec.density=1.0;ec.cp=1000.0;ec.conductivity=0.0;
        ec.max_iterations=20;ec.tolerance=1e-10;
        RosselandSolveControls rc;rc.max_iterations=20;rc.tolerance=1e-10;
        const auto r=solve_rosseland_energy_with_scattering(
            m,g,phi,T,S,a,s,0.0,ec,rc,bc);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(T(0),1000.0,1e-8);
    });

    run_case("p1_isothermal_blackbody_equilibrium",[] {
        Mesh m=make_unit_cube();auto g=build_fv_geometry(m);
        Field<double,Location::CELL> T(1,"T","K",1),G(1,"G","W/m2",1),q(1,"qrad","W/m3",1);
        T(0)=1000.0; G(0)=4.0*M_PI*blackbody_intensity(1000.0); q(0)=0.0;
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,4.0*M_PI*blackbody_intensity(1000.0),0.0};
        const auto r=solve_p1_radiation(m,g,T,G,q,1.0,0.0,bc,2000,1e-12);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(q(0),0.0,1e-8*blackbody_emissive_power(1000.0));
        EXPECT_NEAR(G(0),4.0*M_PI*blackbody_intensity(1000.0),1e-10*G(0));
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
