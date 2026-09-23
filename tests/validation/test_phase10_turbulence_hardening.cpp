#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/spalart_allmaras.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

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
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL; wall.face_ids={0,1,2,3,4,5}; m.boundary().add_patch(wall);
    return m;
}

int main() {
    run_case("M2_SA_canonical_closures", [] {
        SpalartAllmarasModel sa;
        EXPECT_NEAR(sa.fv1(1.0), 1.0/(1.0+7.1*7.1*7.1), 1e-14);
        EXPECT_TRUE(sa.fv2(1.0) >= 0.0 && sa.fv2(1.0) <= 1.0);
        EXPECT_TRUE(sa.ft2(1.0) > 0.0);
        EXPECT_NEAR(sa.cw1, sa.cb1/(sa.kappa*sa.kappa)+(1.0+sa.cb2)/sa.sigma,1e-14);
        EXPECT_NEAR(sa.destruction_coefficient(0.0), 0.0, 1e-14);
        EXPECT_TRUE(sa.turbulent_viscosity(1.0,1e-4,1e-5)>0.0);
    });

    run_case("M2_RNG_reference_constants_and_correction", [] {
        TurbulenceTransportControls c; c.model=TurbulenceModel::RNG_KEPSILON;
        EXPECT_NEAR(c.rng_C_mu,0.0845,1e-14);
        const double eta=2.0;
        const double c1star=c.rng_C1-eta*(1.0-eta/c.rng_eta0)/(1.0+c.rng_beta*eta*eta*eta);
        EXPECT_TRUE(c1star < c.rng_C1);
        EXPECT_TRUE(std::isfinite(c1star));
    });

    run_case("M2_SST_blending_cross_diffusion_bounds", [] {
        const auto near_wall=compute_sst_blending(1.0,10.0,1e-4,1e-5,0.09,1.0e-3);
        const auto far_wall=compute_sst_blending(1.0,1.0,10.0,1e-5,0.09,1.0e-6);
        EXPECT_TRUE(near_wall.first>=0.0 && near_wall.first<=1.0);
        EXPECT_TRUE(near_wall.second>=0.0 && near_wall.second<=1.0);
        EXPECT_TRUE(far_wall.first>=0.0 && far_wall.first<=1.0);
    });

    run_case("M2_komega_transport", [] {
        Mesh m=make_unit_cube(); auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),w(1,"omega","1/s",1),S(1,"S","1/s",1);
        k(0)=0.1; w(0)=2.0; S(0)=1.0;
        TurbulenceTransportControls c; c.model=TurbulenceModel::KOMEGA; c.molecular_viscosity=1e-5;
        ScalarBoundaryConditions bc; bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-4,0.0};
        const auto r=solve_komega_transport(m,g,flux,k,w,S,c,bc,bc,10,1e-6);
        EXPECT_TRUE(r.iterations>0); EXPECT_TRUE(k(0)>=c.k_min); EXPECT_TRUE(w(0)>=c.omega_min);
    });

    run_case("M2_rng_kepsilon_transport", [] {
        Mesh m=make_unit_cube(); auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),e(1,"epsilon","m2/s3",1),S(1,"S","1/s",1);
        k(0)=0.1; e(0)=0.05; S(0)=2.0;
        TurbulenceTransportControls c; c.model=TurbulenceModel::RNG_KEPSILON; c.molecular_viscosity=1e-5;
        ScalarBoundaryConditions bc; bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-5,0.0};
        const auto r=solve_rng_kepsilon_transport(m,g,flux,k,e,S,c,bc,bc,10,1e-6);
        EXPECT_TRUE(r.iterations>0); EXPECT_TRUE(k(0)>=c.k_min); EXPECT_TRUE(e(0)>=c.epsilon_min);
    });

    run_case("M2_spalart_allmaras_transport", [] {
        Mesh m=make_unit_cube(); auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> nut(1,"nu_tilde","m2/s",1),vort(1,"vorticity","1/s",1),wall(1,"wall_distance","m",1);
        nut(0)=1e-4; vort(0)=10.0; wall(0)=0.05;
        TurbulenceTransportControls c; c.model=TurbulenceModel::SPALART_ALLMARAS; c.molecular_viscosity=1e-5;
        ScalarBoundaryConditions bc; bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};
        const auto r=solve_spalart_allmaras_transport(m,g,flux,nut,vort,wall,c,bc,10,1e-6);
        EXPECT_TRUE(r.iterations>0); EXPECT_TRUE(std::isfinite(nut(0))); EXPECT_TRUE(nut(0)>=0.0);
    });

    run_case("M2_controls_reject_nonfinite_values", [] {
        TurbulenceTransportControls c; c.k_min=std::numeric_limits<double>::quiet_NaN();
        bool rejected=false; try { validate_turbulence_controls(c); } catch(const std::invalid_argument&) { rejected=true; }
        EXPECT_TRUE(rejected);
    });
    return run_all();
}
