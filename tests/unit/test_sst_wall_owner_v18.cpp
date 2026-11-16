#include "../audit/chain_fixture.h"
#include "cfdx/physics/sst_wall_operator_v17.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("v18_exact_owner_row_and_original_equation_preserved_elsewhere", [] {
        ScalarEquation equation;
        cfdx::core::SparseMatrix matrix(3,3);
        for (std::size_t i=0;i<3;++i) {
            matrix.push_back(i,i,4.);
            if(i>0) matrix.push_back(i,i-1,-1.);
            if(i+1<3) matrix.push_back(i,i+1,-1.);
        }
        matrix.finalize();
        equation.matrix=matrix;
        equation.rhs=cfdx::core::Vector(3,2.);
        equation.diagonal.assign(3,4.);
        const auto unchanged=equation;
        SSTWallProductionV17 wall;
        wall.g_by_cell.assign(3,0.25);
        wall.has_wall_cell.assign(3,1);
        wall.omega_owner_by_cell={7.,0.,9.};
        wall.has_omega_owner_cell={1,0,1};
        EXPECT_TRUE(sst_apply_omega_owner_constraints_v18(equation,&wall)==2);
        EXPECT_NEAR(equation.matrix(0,0),1.,1e-14);
        EXPECT_NEAR(equation.matrix(0,1),0.,1e-14);
        EXPECT_NEAR(equation.rhs(0),7.,1e-14);
        EXPECT_NEAR(equation.matrix(1,0),unchanged.matrix(1,0),1e-14);
        EXPECT_NEAR(equation.matrix(1,1),unchanged.matrix(1,1),1e-14);
        EXPECT_NEAR(equation.matrix(1,2),unchanged.matrix(1,2),1e-14);
        EXPECT_NEAR(equation.rhs(1),unchanged.rhs(1),1e-14);
        EXPECT_NEAR(equation.matrix(2,2),1.,1e-14);
        EXPECT_NEAR(equation.rhs(2),9.,1e-14);
        cfdx::core::Field<double,cfdx::core::Location::CELL> omega(3,"omega","1/s",1);
        omega(0)=7.;omega(1)=5.;omega(2)=9.;
        const auto audit=sst_wall_owner_audit_v18(omega,&wall);
        EXPECT_TRUE(audit.owners==2);
        EXPECT_NEAR(audit.max_absolute_mismatch,0.,1e-14);
        omega(2)=10.;
        const auto mismatch=sst_wall_owner_audit_v18(omega,&wall);
        EXPECT_NEAR(mismatch.max_absolute_mismatch,1.,1e-14);
        EXPECT_TRUE(mismatch.worst_cell==2);
        // No opt-in must leave all bytes of the equation unchanged.
        ScalarEquation legacy=unchanged;
        SSTWallProductionV17 base;
        base.g_by_cell.assign(3,0.);
        base.has_wall_cell.assign(3,0);
        EXPECT_TRUE(sst_apply_omega_owner_constraints_v18(legacy,&base)==0);
        EXPECT_NEAR(legacy.matrix(0,1),unchanged.matrix(0,1),1e-14);
    });
    run_case("v18_reject_nonfinite_malformed_or_missing_diagonal_atomically", [] {
        ScalarEquation equation;
        cfdx::core::SparseMatrix matrix(2,2);
        matrix.push_back(0,0,2.);matrix.push_back(0,1,-1.);
        matrix.push_back(1,0,-1.);matrix.finalize(); // row 1 has no diagonal
        equation.matrix=matrix;
        equation.rhs=cfdx::core::Vector(2,1.);
        equation.diagonal.assign(2,2.);
        SSTWallProductionV17 wall;
        wall.g_by_cell.assign(2,0.);
        wall.has_wall_cell.assign(2,1);
        wall.omega_owner_by_cell={3.,4.};
        wall.has_omega_owner_cell={1,1};
        EXPECT_THROW(sst_apply_omega_owner_constraints_v18(equation,&wall),std::invalid_argument);
        EXPECT_NEAR(equation.matrix(0,0),2.,1e-14);
        EXPECT_NEAR(equation.rhs(0),1.,1e-14);
        wall.has_omega_owner_cell={1};
        EXPECT_THROW(sst_apply_omega_owner_constraints_v18(equation,&wall),std::invalid_argument);
        wall.has_omega_owner_cell={1,0};
        wall.omega_owner_by_cell[0]=std::numeric_limits<double>::infinity();
        EXPECT_THROW(sst_apply_omega_owner_constraints_v18(equation,&wall),std::invalid_argument);
    });
    run_case("v18_shared_corner_is_weighted_once_and_staged_transactionally", [] {
        const auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        SSTProductionState state(3);
        state.k.fill(0.05);state.omega.fill(2.);
        state.nu_t.fill(0.005);state.strain_rate.fill(2.);
        state.wall_distance.fill(0.3);
        cfdx::core::Field<double,cfdx::core::Location::CELL> U(3,"U","m/s",3);
        U.fill(0.);
        for(std::size_t i=0;i<3;++i) U(i,0)=2.;
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;c.molecular_viscosity=1e-5;c.density=1.;
        const auto triple=prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"wall"},U,state,{{"wall",{0.,0.,0.}}},c);
        const auto disabled=make_sst_wall_operator_v17(mesh,geometry,{"wall"},triple,c);
        const auto enabled=make_sst_wall_operator_v17(mesh,geometry,{"wall"},triple,c,true);
        EXPECT_TRUE(!disabled.production.has_owner_omega_constraints());
        EXPECT_TRUE(enabled.production.has_owner_omega_constraints());
        EXPECT_TRUE(enabled.faces==triple.wall_faces);
        EXPECT_TRUE(enabled.owners==triple.wall_cells);
        EXPECT_TRUE(triple.corner_cells>=1);
        for (std::size_t i=0;i<3;++i) {
            EXPECT_NEAR(enabled.production.omega_owner_by_cell[i],triple.omega_by_cell[i],1e-14);
            EXPECT_TRUE(enabled.production.has_omega_owner_cell[i]==triple.has_wall_cell[i]);
        }
        cfdx::core::Field<double,cfdx::core::Location::CELL> omega(3,"omega","1/s",1);
        for (std::size_t i=0;i<3;++i) omega(i)=triple.omega_by_cell[i];
        EXPECT_NEAR(sst_wall_owner_audit_v18(omega,&enabled.production).max_absolute_mismatch,0.,1e-14);
    });
    run_case("v18_returned_certificate_uses_same_owner_rows_not_old_PDE", [] {
        const auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        cfdx::core::Field<double,cfdx::core::Location::FACE> flux(mesh.n_faces(),"flux","kg/s",1);
        flux.fill(0.);
        cfdx::core::Field<double,cfdx::core::Location::CELL>
            k(3,"k","m2/s2",1), omega(3,"omega","1/s",1),
            strain(3,"S","1/s",1),distance(3,"d","m",1);
        k.fill(0.1);omega.fill(2.);strain.fill(1.);distance.fill(0.3);
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        c.density=1.;c.molecular_viscosity=1e-5;
        ScalarBoundaryConditions kb{{"wall",{ScalarBoundaryType::ZERO_GRADIENT,0.,0.}}};
        ScalarBoundaryConditions wb{{"wall",{ScalarBoundaryType::FIXED_VALUE,5.,0.}}};
        ScalarBoundaryFaceValues values;
        values.values["wall"].assign(mesh.n_faces(),5.);
        SSTWallProductionV17 wall;
        wall.g_by_cell.assign(3,0.1);
        wall.has_wall_cell.assign(3,1);
        const auto old=sst_reassembled_returned_state_residual(
            mesh,geometry,flux,k,omega,strain,distance,c,kb,wb,&values,&wall);
        wall.omega_owner_by_cell.assign(3,2.);
        wall.has_omega_owner_cell.assign(3,1);
        const auto new_sample=sst_reassembled_returned_state_residual(
            mesh,geometry,flux,k,omega,strain,distance,c,kb,wb,&values,&wall);
        EXPECT_TRUE(old.omega_raw_inf>1e-5);
        EXPECT_NEAR(new_sample.omega_raw_inf,0.,1e-12);
        EXPECT_NEAR(new_sample.omega_row_scaled,0.,1e-12);
        const auto b=sst_reassembled_discrete_cell_budgets(
            mesh,geometry,flux,k,omega,strain,distance,c,kb,wb,&values,0,&wall);
        EXPECT_TRUE(b.second.owner_constraint_row);
        EXPECT_NEAR(b.second.owner_constraint_target,2.,1e-12);
        EXPECT_NEAR(b.second.full_residual,0.,1e-12);
        EXPECT_NEAR(b.second.closure_error,0.,1e-12);
    });
    return run_all();
}
