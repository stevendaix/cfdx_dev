#include "cfdx/physics/sst_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

using namespace cfdx::physics;
using namespace cfdx::testing;
using CellField = cfdx::core::Field<double,cfdx::core::Location::CELL>;

int main() {
    run_case("v21_k_change_owner_omega_and_G_same_candidate", [] {
        CellField k(3,"k","m2/s2",1),omega(3,"omega","1/s",1);
        k(0)=4.;k(1)=1.;k(2)=9.;
        omega(0)=100.;omega(1)=15.;omega(2)=20.;
        ScalarBoundaryFaceValues faces;
        SSTWallProductionV17 source;
        std::size_t calls=0;
        const SSTWallRefreshV17 callback=[&](const auto& candidate_k,const auto& candidate_w,
                                              auto& output_faces, auto& output_source) {
            ++calls;
            output_source.g_by_cell={0.01*candidate_w(0),0.,0.01*candidate_w(2)};
            output_source.has_wall_cell={1,0,1};
            output_source.has_omega_owner_cell={1,0,1};
            output_source.omega_owner_by_cell={
                96.399 + 100.*candidate_k(0),0.,8. + 10.*candidate_k(2)};
            output_faces.values["wall"]={candidate_w(0),candidate_w(2)};
        };
        const auto result=sst_synchronize_wall_owner_v21(
            k,omega,callback,faces,source,1e-10);
        EXPECT_TRUE(calls==2);
        EXPECT_TRUE(result.owners==2);
        EXPECT_NEAR(result.maximum_before,396.399,1e-10);
        EXPECT_NEAR(result.maximum_after,0.,1e-12);
        EXPECT_NEAR(result.relative_after,0.,1e-12);
        EXPECT_NEAR(omega(0),496.399,1e-10);
        EXPECT_NEAR(omega(1),15.,1e-10);
        EXPECT_NEAR(omega(2),98.,1e-10);
        // G is recomputed from the accepted omega, not the stale pre-step one.
        EXPECT_NEAR(source.g_by_cell[0],4.96399,1e-10);
        EXPECT_NEAR(source.g_by_cell[2],0.98,1e-12);
        EXPECT_NEAR(faces.values.at("wall")[0],496.399,1e-10);
        EXPECT_TRUE(result.worst_cell==0);
    });
    run_case("v21_noncontractive_owner_target_rolls_back_every_output", [] {
        CellField k(1,"k","m2/s2",1),omega(1,"omega","1/s",1);
        k(0)=1.;omega(0)=3.;
        ScalarBoundaryFaceValues faces;
        faces.values["sentinel"]={7.};
        SSTWallProductionV17 wall;
        wall.g_by_cell={0.};wall.has_wall_cell={0};
        const SSTWallRefreshV17 callback=[](const auto&,const auto& w,auto&,auto& p) {
            p.g_by_cell={0.5};p.has_wall_cell={1};
            p.has_omega_owner_cell={1};p.omega_owner_by_cell={w(0)+1.};
        };
        EXPECT_THROW(sst_synchronize_wall_owner_v21(
            k,omega,callback,faces,wall,1e-8),std::domain_error);
        EXPECT_NEAR(omega(0),3.,1e-14);
        EXPECT_NEAR(wall.g_by_cell[0],0.,1e-14);
        EXPECT_NEAR(faces.values.at("sentinel")[0],7.,1e-14);
        EXPECT_TRUE(!wall.has_owner_omega_constraints());
    });
    run_case("v21_invalid_floor_or_target_rejected_transactionally", [] {
        CellField k(1,"k","m2/s2",1),omega(1,"omega","1/s",1);
        k(0)=1.;omega(0)=5.;
        ScalarBoundaryFaceValues faces;
        SSTWallProductionV17 wall;
        const SSTWallRefreshV17 callback=[](const auto&,const auto&,auto&,auto& p) {
            p.g_by_cell={0.};p.has_wall_cell={1};
            p.has_omega_owner_cell={1};p.omega_owner_by_cell={1e-30};
        };
        EXPECT_THROW(sst_synchronize_wall_owner_v21(
            k,omega,callback,faces,wall,1e-8),std::domain_error);
        EXPECT_THROW(sst_synchronize_wall_owner_v21(
            k,omega,callback,faces,wall,0.),std::invalid_argument);
        EXPECT_NEAR(omega(0),5.,1e-14);
    });
    run_case("v21_no_owner_rows_does_not_modify_omega", [] {
        CellField k(1,"k","m2/s2",1),omega(1,"omega","1/s",1);
        k(0)=1.;omega(0)=5.;
        ScalarBoundaryFaceValues faces;
        SSTWallProductionV17 wall;
        std::size_t calls=0;
        const SSTWallRefreshV17 callback=[&](const auto&,const auto&,auto&,auto& p) {
            ++calls;
            p.g_by_cell={2.};p.has_wall_cell={1};
        };
        const auto r=sst_synchronize_wall_owner_v21(
            k,omega,callback,faces,wall,1e-8);
        EXPECT_TRUE(calls==1);
        EXPECT_TRUE(r.owners==0);
        EXPECT_NEAR(omega(0),5.,1e-14);
        EXPECT_NEAR(wall.g_by_cell[0],2.,1e-14);
    });
    return run_all();
}
