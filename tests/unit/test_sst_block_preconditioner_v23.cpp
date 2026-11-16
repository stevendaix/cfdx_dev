#include "cfdx/physics/sst_solver.h"
#include "common/test_harness.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

static cfdx::core::Mesh one_cell_cube_v23() {
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double coords[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                               {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)
        mesh.points().set(i,coords[i][0],coords[i][1],coords[i][2]);
    const std::vector<std::vector<FaceIndex>> faces={{0,3,2,1},{4,5,6,7},
        {0,1,5,4},{3,7,6,2},{0,4,7,3},{1,2,6,5}};
    for(const auto& face:faces)mesh.faces().push_face(face);
    mesh.ownership().resize(faces.size());
    for(std::size_t f=0;f<faces.size();++f){
        mesh.ownership().set_owner(f,0);
        mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    Patch wall;wall.name="wall";wall.type=PatchType::WALL;
    wall.face_ids={0,1,2,3,4,5};mesh.boundary().add_patch(wall);
    return mesh;
}

int main() {
    run_case("v23_block_inverse_captures_cross_equation_terms", [] {
        const SSTBlockPreconditionerV23 p({{4.,1.,2.,3.},{1.,0.5,0.25,2.}});
        const auto x=p.apply({9.,2.,8.,4.}); // [k0,k1,w0,w1]
        EXPECT_TRUE(p.coupled_blocks()==2 && p.diagonal_fallbacks()==0);
        EXPECT_NEAR(4*x[0]+x[2],9.,1e-12);
        EXPECT_NEAR(2*x[0]+3*x[2],8.,1e-12);
        EXPECT_NEAR(x[1]+0.5*x[3],2.,1e-12);
        EXPECT_NEAR(0.25*x[1]+2*x[3],4.,1e-12);
    });
    run_case("v23_singular_block_falls_back_to_invertible_identity", [] {
        const SSTBlockPreconditionerV23 p({{1.,2.,2.,4.}});
        EXPECT_TRUE(p.diagonal_fallbacks()==1);
        const auto x=p.apply({3.,5.});
        EXPECT_NEAR(x[0],3.,0.);
        EXPECT_NEAR(x[1],5.,0.);
        EXPECT_THROW(p.apply({1.}),std::invalid_argument);
        EXPECT_THROW(SSTBlockPreconditionerV23({{0.,0.,0.,0.}}),std::invalid_argument);
    });
    run_case("v23_preconditioned_newton_checks_original_directional_residual", [] {
        // A stiff off-diagonal system with known positive physical root.
        // The block matrix is an approximate Jacobian in logarithmic variables
        // at the starting state, not a replacement for the residual F.
        const SSTPhysicalResidualV22 F=[](const auto& x) {
            return std::vector<double>{1000.*(x[0]-2.)+500.*(x[1]-3.),
                                       100.*(x[0]-2.)+1000.*(x[1]-3.)};
        };
        CoupledNewtonKrylovControls ctrl;
        ctrl.krylov_dimension=2;ctrl.linear_relative_tolerance=0.25;
        ctrl.finite_difference_step=1e-7;
        const SSTBlockPreconditionerV23 pre({{1000.,1000.,100.,2000.}});
        const auto step=sst_jfnk_propose_v22(
            {1.,2.},{1.,1.},1,1e-10,1e-10,F,ctrl,{},
            [&pre](const auto& x){return pre.apply(x);});
        EXPECT_TRUE(step.available);
        EXPECT_TRUE(step.linear_true_relative<=ctrl.linear_relative_tolerance);
        EXPECT_TRUE(step.frozen_merit_after<step.frozen_merit_before);
        EXPECT_TRUE(step.fields[0]>0. && step.fields[1]>0.);
        EXPECT_TRUE(step.function_evaluations>2);
    });
    run_case("v23_fv_source_blocks_use_existing_sst_assembly", [] {
        using namespace cfdx::core;
        const auto mesh=one_cell_cube_v23();
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"m","kg/s",1);
        flux.fill(0.);
        Field<double,Location::CELL> k(1,"k","",1),omega(1,"w","",1),
            strain(1,"s","",1),dist(1,"d","",1);
        k(0)=0.4;omega(0)=8.;strain(0)=2.;dist(0)=0.1;
        TurbulenceTransportControls cfg;cfg.model=TurbulenceModel::SST;
        ScalarBoundaryConditions kb,wb;
        kb["wall"]={ScalarBoundaryType::FIXED_VALUE,0.01,0.};
        wb["wall"]={ScalarBoundaryType::FIXED_VALUE,2.,0.};
        const auto full=sst_coupled_residual_v22(
            mesh,geom,flux,k,omega,strain,dist,cfg,kb,wb,nullptr);
        std::vector<double> scales(2),left(2);
        for(std::size_t i=0;i<2;++i){
            scales[i]=full.row_norms[i]*full.local_scales[i];
            left[i]=std::clamp(std::abs(
                full.diagonal[i]*(i==0?k(0):omega(0))/scales[i]),0.02,50.);
        }
        const auto pre=sst_local_source_blocks_v23(
            mesh,geom,k,omega,strain,dist,cfg,kb,wb,nullptr,nullptr,
            full,scales,left);
        EXPECT_TRUE(pre.cells()==1);
        EXPECT_TRUE(pre.coupled_blocks()==1);
        const auto result=pre.apply({1.,0.1});
        EXPECT_TRUE(result.size()==2 && std::isfinite(result[0]) &&
                    std::isfinite(result[1]));
    });
    run_case("v23_option_requires_v22_and_is_default_off", [] {
        TurbulenceTransportControls controls;
        EXPECT_TRUE(!controls.sst_block_preconditioner_v23_experimental);
        controls.sst_block_preconditioner_v23_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(controls),std::invalid_argument);
        controls.sst_jfnk_v22_experimental=true;
        validate_turbulence_controls(controls);
    });
    return run_all();
}
