#include "cfdx/physics/sst_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;

static cfdx::core::Mesh one_cell_fv_cube() {
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double pts[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                            {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i) mesh.points().set(i,pts[i][0],pts[i][1],pts[i][2]);
    const std::vector<std::vector<FaceIndex>> faces={
        {0,3,2,1},{4,5,6,7},{0,1,5,4},
        {3,7,6,2},{0,4,7,3},{1,2,6,5}};
    for(const auto& face:faces) mesh.faces().push_face(face);
    mesh.ownership().resize(faces.size());
    for(std::size_t f=0;f<faces.size();++f) {
        mesh.ownership().set_owner(f,0);
        mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    Patch wall;wall.name="wall";wall.type=PatchType::WALL;
    wall.face_ids={0,1,2,3,4,5};mesh.boundary().add_patch(wall);
    return mesh;
}

int main() {
    run_case("v22_newton_log_unknowns_couple_k_omega_and_reduce_true_residual", [] {
        const SSTPhysicalResidualV22 F=[](const std::vector<double>& v) {
            return std::vector<double>{v[0]+0.1*v[1]-2.0,v[1]-3.0*v[0]};
        };
        CoupledNewtonKrylovControls ctl;
        ctl.krylov_dimension=2;ctl.linear_relative_tolerance=0.5;
        ctl.finite_difference_step=1e-7;
        const auto result=sst_jfnk_propose_v22({1.,4.},{1.,4.},1,1e-12,1e-12,F,ctl);
        EXPECT_TRUE(result.available);
        EXPECT_TRUE(result.function_evaluations>=5);
        EXPECT_TRUE(result.linear_true_relative<=0.5);
        EXPECT_TRUE(result.frozen_merit_after<result.frozen_merit_before);
        EXPECT_TRUE(result.fields[0]>0. && result.fields[1]>0.);
        // A Newton step in log(k),log(omega) is not the exact root after
        // one iteration even when F is linear in physical units.
        auto x=result.fields;
        for (int it=0;it<6;++it) {
            const auto f=F(x);
            if (std::max(std::abs(f[0]),std::abs(f[1]))<1e-12) break;
            const auto next=sst_jfnk_propose_v22(x,{1.,4.},1,1e-12,1e-12,F,ctl);
            EXPECT_TRUE(next.available);
            x=next.fields;
        }
        EXPECT_NEAR(x[0],2./1.3,1e-10);
        EXPECT_NEAR(x[1],6./1.3,1e-10);
    });
    run_case("v22_jacobi_left_preconditioner_preserves_coupled_root", [] {
        const SSTPhysicalResidualV22 F=[](const std::vector<double>& v) {
            return std::vector<double>{0.01*(v[0]+0.1*v[1]-2.),
                                       100.*(v[1]-3.*v[0])};
        };
        CoupledNewtonKrylovControls ctl;
        ctl.krylov_dimension=2;ctl.linear_relative_tolerance=0.5;
        auto x=std::vector<double>{1.,4.};
        for (int i=0;i<7;++i) {
            const auto next=sst_jfnk_propose_v22(x,{1.,1.},1,1e-12,1e-12,
                                                  F,ctl,{0.01,100.});
            EXPECT_TRUE(next.available);
            x=next.fields;
            const auto f=F(x);
            if (std::max(std::abs(f[0]),std::abs(f[1]))<1e-10) break;
        }
        EXPECT_NEAR(x[0],2./1.3,1e-9);
        EXPECT_NEAR(x[1],6./1.3,1e-9);
        EXPECT_THROW(sst_jfnk_propose_v22(x,{1.,1.},1,1e-12,1e-12,
            F,ctl,{1.,0.}),std::invalid_argument);
    });
    run_case("v22_near_floor_rejects_without_modifying_physical_state", [] {
        const auto initial=std::vector<double>{1e-12,8.};
        const auto r=sst_jfnk_propose_v22(initial,{1.,1.},1,1e-12,1e-12,
            [](const auto& v){return std::vector<double>{v[0]-1.,v[1]-2.};});
        EXPECT_TRUE(!r.available);
        EXPECT_TRUE(r.fields.empty());
        EXPECT_TRUE(initial[0]==1e-12 && initial[1]==8.);
    });
    run_case("v22_rejects_incomplete_vector_operator", [] {
        CoupledNewtonKrylovControls ctl;
        ctl.krylov_dimension=2;ctl.linear_relative_tolerance=0.4;
        const auto r=sst_jfnk_propose_v22({1.,4.},{1.,1.},1,1e-12,1e-12,
            [](const auto&){return std::vector<double>{0.};},ctl);
        EXPECT_TRUE(!r.available);
        EXPECT_TRUE(r.fields.empty());
        EXPECT_TRUE(!r.reason.empty());
        EXPECT_THROW(sst_jfnk_propose_v22({1.,4.},{1.,0.},1,1e-12,1e-12,
            [](const auto& v){return v;},ctl),std::invalid_argument);
    });
    run_case("v22_fv_residual_is_exactly_the_production_reassembled_system", [] {
        using namespace cfdx::core;
        const auto mesh=one_cell_fv_cube();
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"flux","kg/s",1);
        flux.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),w(1,"w","1/s",1),
            strain(1,"strain","1/s",1),distance(1,"wall_distance","m",1);
        k(0)=0.4;w(0)=8.;strain(0)=2.;distance(0)=0.1;
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        ScalarBoundaryConditions kb,wb;
        kb["wall"]={ScalarBoundaryType::FIXED_VALUE,0.01,0.};
        wb["wall"]={ScalarBoundaryType::FIXED_VALUE,2.,0.};
        const auto full=sst_coupled_residual_v22(mesh,geom,flux,k,w,strain,
                                                 distance,c,kb,wb,nullptr);
        const auto sample=sst_reassembled_returned_state_residual(
            mesh,geom,flux,k,w,strain,distance,c,kb,wb,nullptr);
        EXPECT_TRUE(full.values.size()==2 && full.row_norms.size()==2);
        EXPECT_TRUE(std::isfinite(full.values[0]) && std::isfinite(full.values[1]));
        EXPECT_NEAR(std::abs(full.values[0]),sample.k_raw_inf,1e-11);
        EXPECT_NEAR(std::abs(full.values[1]),sample.omega_raw_inf,1e-11);
        EXPECT_NEAR(std::abs(full.values[0])/full.row_norms[0],sample.k_row_scaled,1e-11);
        EXPECT_NEAR(std::abs(full.values[1])/full.row_norms[1],sample.omega_row_scaled,1e-11);
        EXPECT_NEAR(std::abs(full.values[0])/(full.row_norms[0]*full.local_scales[0]),
                    sample.k_local_relative,1e-11);
        EXPECT_NEAR(std::abs(full.values[1])/(full.row_norms[1]*full.local_scales[1]),
                    sample.omega_local_relative,1e-11);
    });
    return run_all();
}
