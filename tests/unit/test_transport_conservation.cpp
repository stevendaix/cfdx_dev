#include "cfdx/physics/conservation_boundedness.h"
#include "cfdx/physics/energy_solver.h"
#include "common/test_harness.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include <limits>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh one_cell_mesh() {
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    for (std::size_t f = 0; f < 6; ++f) wall.face_ids.push_back(f);
    m.boundary().add_patch(wall);
    return m;
}

int main() {
    run_case("source_balance_closes_boundary_flux", [] {
        const Mesh m=one_cell_mesh();
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> source(m.n_cells(),"source","kg/s",1);
        phi.fill(0.0); source.fill(0.0);
        phi(0)=2.0; source(0)=-2.0;
        const auto r=audit_transport_balance(m,phi,source);
        EXPECT_NEAR(r.boundary_flux,2.0,1e-15);
        EXPECT_NEAR(r.source_sum,-2.0,1e-15);
        EXPECT_NEAR(r.residual,0.0,1e-15);
        EXPECT_NEAR(r.normalized_residual,0.0,1e-15);
    });
    run_case("transient_accumulation_is_explicit", [] {
        const Mesh m=one_cell_mesh();
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> source(m.n_cells(),"source","kg/s",1);
        Field<double,Location::CELL> accumulation(m.n_cells(),"dMdt","kg/s",1);
        phi.fill(0.0); source.fill(0.0); accumulation.fill(0.0);
        phi(0)=3.0; accumulation(0)=3.0;
        const auto r=audit_transport_balance(m,phi,source,&accumulation);
        EXPECT_NEAR(r.residual,0.0,1e-15);
    });
    run_case("scalar_bounds_report_worst_cell_and_nonfinite", [] {
        Field<double,Location::CELL> k(4,"k","m2/s2",1);
        k(0)=0.1; k(1)=-0.02; k(2)=1.5; k(3)=std::numeric_limits<double>::quiet_NaN();
        const auto r=audit_scalar_bounds(k,0.0,1.0);
        EXPECT_TRUE(!r.bounded());
        EXPECT_TRUE(r.below==1);
        EXPECT_TRUE(r.above==1);
        EXPECT_TRUE(r.nonfinite==1);
        EXPECT_TRUE(r.worst_violation>0.0);
    });
    run_case("strict_turbulence_floor_is_model_selectable", [] {
        Field<double,Location::CELL> omega(2,"omega","1/s",1);
        omega(0)=1e-8; omega(1)=2e-3;
        EXPECT_TRUE(audit_positive_scalar(omega,1e-10).bounded());
        omega(0)=0.0;
        EXPECT_TRUE(!audit_positive_scalar(omega,1e-10).bounded());
    });

    run_case("independent_scalar_flux_reconstruction_is_face_conservative", [] {
        Mesh m = one_cell_mesh();
        auto g = build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> q(m.n_cells(),"q","1",1);
        phi.fill(0.0); q(0)=2.0;
        phi(0)=3.0;
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,1.0,0.0};
        const auto reconstructed =
            reconstruct_scalar_transport_flux(
                m,g,phi,q,0.0,bc,true,ConvectionScheme::UPWIND);
        EXPECT_NEAR(reconstructed(0),3.0,1e-15);
        const auto audit=audit_face_flux_conservation(m,reconstructed);
        EXPECT_NEAR(audit.global_boundary_flux,3.0,1e-15);
        EXPECT_NEAR(audit.global_cell_balance,3.0,1e-15);
    });


    run_case("independent_energy_balance_is_zero_for_closed_constant_state", [] {
        Mesh m = one_cell_mesh();
        auto g = build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(m.n_cells(),"T","K",1);
        Field<double,Location::CELL> Told(m.n_cells(),"Told","K",1);
        Field<double,Location::CELL> source(m.n_cells(),"source","W/m3",1);
        phi.fill(0.0); T(0)=300.0; Told(0)=300.0; source.fill(0.0);
        EnergySolverControls controls;
        controls.density=1.0; controls.cp=1000.0; controls.conductivity=1.0; controls.dt=0.0;
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        const auto r=reconstruct_energy_balance(m,g,phi,T,Told,source,controls,bc);
        EXPECT_NEAR(r.residual,0.0,1e-15);
        EXPECT_NEAR(r.normalized_residual,0.0,1e-15);
        EXPECT_TRUE(r.nonfinite_faces==0);
        EXPECT_TRUE(r.nonfinite_cells==0);
    });

    return run_all();
}
