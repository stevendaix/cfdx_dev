#include "../audit/chain_fixture.h"
#include "cfdx/physics/sst_wall_operator_v17.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("v17_atomic_triplet_same_state_with_density_units_and_corner", [] {
        const auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        SSTProductionState state(3);
        state.k.fill(0.05); state.omega.fill(2.0);
        state.nu_t.fill(0.005); state.strain_rate.fill(2.0);
        state.wall_distance.fill(0.3);
        cfdx::core::Field<double,cfdx::core::Location::CELL> U(3,"U","m/s",3);
        U.fill(0.0);
        for(std::size_t i=0;i<3;++i) U(i,0)=2.;
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        c.molecular_viscosity=1e-5; c.density=2.;
        const auto triplet=prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"wall"},U,state,{{"wall",{0.,0.,0.}}},c);
        auto op=make_sst_wall_operator_v17(mesh,geometry,{"wall"},triplet,c);
        EXPECT_TRUE(op.faces==12 && op.owners==3);
        EXPECT_TRUE(!op.openfoam_equivalent);
        const auto& face_values=op.omega_faces.values.at("wall");
        for(const auto face:mesh.boundary().patch(mesh.boundary().find("wall")).face_ids){
            EXPECT_NEAR(face_values[face],triplet.omega_by_face[face],1e-12);
            EXPECT_NEAR(op.momentum_face_mu[face],c.density*(c.molecular_viscosity+triplet.nut_by_face[face]),1e-12);
        }
        for(std::size_t i=0;i<mesh.n_cells();++i) {
            EXPECT_TRUE(triplet.omega_by_cell[i]>0.);
            EXPECT_NEAR(op.production.g_by_cell[i],triplet.production_by_cell[i],1e-12);
        }
        c.density=1.;
        const auto ref_unit_density=prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"wall"},U,state,{{"wall",{0.,0.,0.}}},c);
        EXPECT_NEAR(ref_unit_density.production_by_cell[1],triplet.production_by_cell[1],1e-13);
        // The historical V8 wall *audit* must also use mass-specific G,
        // irrespective of rho; previously it accidentally fed rho*G.
        c.density=2.;
        const auto audit_rho2=audit_openfoam_wall_coupling_v8(
            mesh,geometry,{"wall"},U,state,{{"wall",{0.,0.,0.}}},c);
        c.density=1.;
        const auto audit_rho1=audit_openfoam_wall_coupling_v8(
            mesh,geometry,{"wall"},U,state,{{"wall",{0.,0.,0.}}},c);
        EXPECT_NEAR(audit_rho1.max_production_gap,audit_rho2.max_production_gap,1e-13);
        auto missing=triplet;
        missing.has_wall_face[mesh.boundary().patch(mesh.boundary().find("wall")).face_ids.front()]=0;
        EXPECT_THROW(make_sst_wall_operator_v17(mesh,geometry,{"wall"},missing,c),std::invalid_argument);
        EXPECT_THROW(make_sst_wall_operator_v17(mesh,geometry,{"wall","wall"},triplet,c),std::invalid_argument);
    });
    run_case("v17_wall_G_override_changes_both_sst_sources_conservatively", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST; c.density=2.;
        c.molecular_viscosity=1e-5;
        const auto baseline=sst_cell_sources(2.,5.,3.,1.,1.,0.,c);
        const auto wall=sst_cell_sources(2.,5.,3.,1.,1.,0.,c,0.2);
        EXPECT_NEAR(wall.sk,0.4,1e-14);
        EXPECT_TRUE(std::abs(baseline.sk-wall.sk)>0.01);
        EXPECT_TRUE(std::abs(baseline.sw-wall.sw)>0.01);
        EXPECT_THROW(sst_cell_sources(2.,5.,3.,1.,1.,0.,c,
            std::numeric_limits<double>::quiet_NaN()),std::invalid_argument);
    });
    run_case("v17_wall_nut_enters_production_momentum_FV_matrix", [] {
        const auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        cfdx::core::Field<double,cfdx::core::Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        phi.fill(0.0);
        cfdx::core::Field<double,cfdx::core::Location::CELL> grad(mesh.n_cells(),"gp","Pa/m",3);
        grad.fill(0.0);
        cfdx::core::Field<double,cfdx::core::Location::CELL> body(mesh.n_cells(),"b","N/m3",1);
        body.fill(0.0);
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        std::vector<double> cell_mu(mesh.n_cells(),0.002);
        const double mu=1e-5;
        const auto baseline=assemble_momentum_component(mesh,geometry,phi,grad,body,cell_mu,bc,
            0,true,ConvectionScheme::UPWIND,nullptr,nullptr,mu);
        std::vector<double> wall_mu(mesh.n_faces(),-1.);
        for(const auto face:mesh.boundary().patch(mesh.boundary().find("wall")).face_ids)
            wall_mu[face]=0.02;
        const auto coupled=assemble_momentum_component(mesh,geometry,phi,grad,body,cell_mu,bc,
            0,true,ConvectionScheme::UPWIND,nullptr,nullptr,mu,1.,nullptr,&wall_mu);
        double delta=0.;
        for(std::size_t i=0;i<mesh.n_cells();++i)
            for(auto j=baseline.matrix.row_offsets_data()[i];j<baseline.matrix.row_offsets_data()[i+1];++j)
                delta+=std::abs(baseline.matrix.values_data()[j]-coupled.matrix.values_data()[j]);
        EXPECT_TRUE(delta>0.01);
        wall_mu[mesh.boundary().patch(mesh.boundary().find("wall")).face_ids.front()]=-1.;
        EXPECT_THROW(assemble_momentum_component(mesh,geometry,phi,grad,body,cell_mu,bc,
            0,true,ConvectionScheme::UPWIND,nullptr,nullptr,mu,1.,nullptr,&wall_mu),std::invalid_argument);
    });
    run_case("v17_returned_SST_equations_and_cell_budget_use_dynamic_wall_values", [] {
        const auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        cfdx::core::Field<double,cfdx::core::Location::FACE> phi(
            mesh.n_faces(),"flux","kg/s",1);
        phi.fill(0.0);
        cfdx::core::Field<double,cfdx::core::Location::CELL>
            k(mesh.n_cells(),"k","m2/s2",1),
            omega(mesh.n_cells(),"omega","1/s",1),
            strain(mesh.n_cells(),"S","1/s",1),
            distance(mesh.n_cells(),"d","m",1);
        k.fill(0.05);omega.fill(2.);strain.fill(1.);distance.fill(0.3);
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST; c.density=2.;
        c.molecular_viscosity=1e-5;
        c.sst_molecular_wall_diffusion=true;
        ScalarBoundaryConditions kb{{"wall",{ScalarBoundaryType::ZERO_GRADIENT,0.,0.}}};
        ScalarBoundaryConditions wb{{"wall",{ScalarBoundaryType::FIXED_VALUE,1.,0.}}};
        const auto baseline=sst_reassembled_returned_state_residual(
            mesh,geometry,phi,k,omega,strain,distance,c,kb,wb,nullptr);
        ScalarBoundaryFaceValues face;
        face.values["wall"].assign(mesh.n_faces(),10.);
        SSTWallProductionV17 G;
        G.g_by_cell.assign(mesh.n_cells(),0.25);
        G.has_wall_cell.assign(mesh.n_cells(),1);
        const auto modified=sst_reassembled_returned_state_residual(
            mesh,geometry,phi,k,omega,strain,distance,c,kb,wb,&face,&G);
        EXPECT_TRUE(std::abs(modified.k_raw_inf-baseline.k_raw_inf)>1e-6);
        EXPECT_TRUE(std::abs(modified.omega_raw_inf-baseline.omega_raw_inf)>1e-6);
        const auto budgets=sst_reassembled_discrete_cell_budgets(
            mesh,geometry,phi,k,omega,strain,distance,c,kb,wb,&face,0,&G);
        EXPECT_NEAR(budgets.first.closure_error,0.0,1e-8);
        EXPECT_NEAR(budgets.second.closure_error,0.0,1e-8);
        EXPECT_THROW(sst_reassembled_discrete_cell_budgets(
            mesh,geometry,phi,k,omega,strain,distance,c,kb,wb,&face,
            mesh.n_cells(),&G),std::invalid_argument);
    });
    return run_all();
}
