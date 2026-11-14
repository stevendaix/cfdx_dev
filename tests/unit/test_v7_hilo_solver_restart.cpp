#include "../audit/chain_fixture.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <cmath>
#include <filesystem>
#include <iostream>

using namespace cfdx::physics;
using namespace cfdx::testing;
namespace {
struct Driver {
    Mesh mesh=make_1d_chain(3);
    Field<double,Location::CELL> U{3,"U","m/s",3}, p{3,"p","Pa",1};
    VelocityBoundaryConditions ubc;
    ScalarBoundaryConditions pbc;
    IncompressibleSolverControls controls;
    Driver() {
        mesh.boundary().patch(0).type=PatchType::INLET;
        mesh.boundary().patch(1).type=PatchType::OUTLET;
        U.fill(0);p.fill(0);
        controls.density=1.2;
        controls.kinematic_viscosity=0.01;
        controls.convergence.max_iterations=2;
        controls.convergence.relative_tolerance=1e-9;
        controls.linear_tolerance=1e-11;
        controls.convergence.continuity_tolerance=1e-11;
        controls.initialization.mode=InitializationMode::Uniform;
        controls.initialization.uniform_velocity={0.8,0,0};
        controls.initialization.uniform_pressure=0;
        controls.v18_sparse_flux_pressure_jacobian=true;
        controls.v19_full_pressure_flux_step=true;
        controls.pressure_high_low_production=true;
        ubc["inlet"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
        ubc["outlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
        ubc["wall"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        pbc["inlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        pbc["outlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
        pbc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    }
};
}
int main() {
    run_case("production_hilo_dat6_restart_matches_uninterrupted_two_steps", [] {
        Driver original;
        const auto path=std::filesystem::temp_directory_path()/"cfdx_v7_prod_hilo.dat";
        bool saved=false;
        original.controls.accepted_compensated_iteration_callback=[&](std::size_t iter,const auto& mesh,
                        const auto& U,const auto& p,const auto& phi,const auto& pressure_low) {
            EXPECT_TRUE(pressure_low.size()==mesh.n_cells());
            EXPECT_TRUE(phi.size()==mesh.n_faces());
            if(iter!=1) return;
            Field<double,Location::CELL> low(mesh.n_cells(),"p_low","Pa",1);
            for(std::size_t cell=0;cell<mesh.n_cells();++cell) low(cell)=pressure_low(cell);
            auto& phi_mutable=const_cast<Field<double,Location::FACE>&>(phi);
            cfdx::io::DatRestartFields fields;
            fields.pressure_low=&low;
            fields.face_mass_flux=&phi_mutable;
            cfdx::io::write_dat_restart_fields(path.string(),mesh,U,p,fields,iter,0.0);
            saved=true;
        };
        const auto continuous=solve_steady_incompressible(
            original.mesh,original.U,original.p,original.ubc,original.pbc,original.controls);
        EXPECT_TRUE(saved);
        std::cout << "continuous_iterations="<< continuous.iterations
                  << " history=" << continuous.history.size() << "\n";
        Driver restarted;
        restarted.controls.initialization.mode=InitializationMode::Restart;
        restarted.controls.convergence.max_iterations=1;
        const auto resumed=solve_steady_incompressible(
            restarted.mesh,restarted.U,restarted.p,restarted.ubc,restarted.pbc,
            restarted.controls,path.string());
        EXPECT_TRUE(resumed.history.size()==1);
        EXPECT_TRUE(resumed.authoritative_pressure_low.size()==3);
        double maxu=0,maxp=0,maxphi=0,maxlo=0;
        for(std::size_t i=0;i<3;++i) {
            for(std::size_t d=0;d<3;++d)
                maxu=std::max(maxu,std::abs(original.U(i,d)-restarted.U(i,d)));
            maxp=std::max(maxp,std::abs(original.p(i)-restarted.p(i)));
            maxlo=std::max(maxlo,std::abs(
                continuous.authoritative_pressure_low(i)-resumed.authoritative_pressure_low(i)));
        }
        for(std::size_t f=0;f<original.mesh.n_faces();++f)
            maxphi=std::max(maxphi,std::abs(continuous.authoritative_mass_flux(f)-
                                           resumed.authoritative_mass_flux(f)));
        std::cout << "replay_max_u=" << maxu << " p=" << maxp
                  << " phi=" << maxphi << " low=" << maxlo << "\n";
        EXPECT_TRUE(maxu < 1e-9 && maxp < 1e-9 && maxphi < 1e-9 && maxlo < 1e-9);
        std::filesystem::remove(path);
    });
    return run_all();
}
