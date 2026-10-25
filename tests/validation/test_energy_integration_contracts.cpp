// Independent energy oracles for implicit sources and variable conductivity.
// The baseline macro runs only contracts expressible by the previous API.
#include "audit/chain_fixture.h"
#include "cfdx/physics/production_energy_stage.h"
#include "common/test_harness.h"
#include <limits>
#ifndef CFDX_THERMAL_CONTRACTS_BASELINE
#include "audit/sst_channel_fixture.h"
#include "cfdx/physics/production_buoyancy_stage.h"
#endif

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static EnergySolverControls energy_controls()
{
    EnergySolverControls c;
    c.density=1.0; c.cp=1.0; c.conductivity=1.0;
    c.relaxation=1.0; c.tolerance=1e-12;
    return c;
}

static cfdx::io::CaseSetup thermal_case()
{
    cfdx::io::CaseSetup setup;
    setup.energy_model="conjugate";
    cfdx::io::MaterialSpec material;
    material.name="test_fluid"; material.density=2.0;
    material.dynamic_viscosity=0.001;
    material.thermal_conductivity=1.0; material.specific_heat=10.0;
    setup.materials.push_back(material);
    cfdx::io::BoundarySpec inlet;
    inlet.patch_name="inlet"; inlet.type=cfdx::io::BCType::WALL;
    inlet.value_type=cfdx::io::BCValueType::WALL_THERMAL;
    inlet.wall_temperature=300.0;
    setup.boundary_conditions.push_back(inlet);
    cfdx::io::BoundarySpec outlet;
    outlet.patch_name="outlet"; outlet.type=cfdx::io::BCType::WALL;
    setup.boundary_conditions.push_back(outlet);
    cfdx::io::BoundarySpec wall;
    wall.patch_name="wall"; wall.type=cfdx::io::BCType::WALL;
    wall.raw_params["heat_flux"]="10";
    setup.boundary_conditions.push_back(wall);
    return setup;
}

int main()
{
    run_case("implicit_source_equilibrium_uses_final_temperature", [] {
        const auto mesh=make_1d_chain(1);
        const auto g=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(1,"T","K",1),Su(1,"Su","W/m3",1),Sp(1,"Sp","W/m3/K",1);
        phi.fill(0.0); T(0)=200.0; Su(0)=600.0; Sp(0)=-2.0;
        auto c=energy_controls(); c.conductivity=0.0;
        const auto r=solve_energy(mesh,g,phi,T,Su,c,{},nullptr,&Sp);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(T(0),300.0,1e-12);
        EXPECT_NEAR(r.history.back().independent_conservation_residual,0.0,1e-12);
        EXPECT_NEAR(r.history.back().independent_conservation_normalized,0.0,1e-12);
    });

    run_case("implicit_source_underrelaxed_state_remains_out_of_balance", [] {
        const auto mesh=make_1d_chain(1);
        const auto g=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(1,"T","K",1),Su(1,"Su","W/m3",1),Sp(1,"Sp","W/m3/K",1);
        phi.fill(0.0); T(0)=200.0; Su(0)=600.0; Sp(0)=-2.0;
        auto c=energy_controls(); c.conductivity=0.0;
        c.relaxation=0.5; c.max_iterations=1;
        const auto r=solve_energy(mesh,g,phi,T,Su,c,{},nullptr,&Sp);
        EXPECT_FALSE(r.converged);
        EXPECT_NEAR(T(0),250.0,1e-12);
        // No boundary flux or accumulation: residual = -(600 - 2*250).
        EXPECT_NEAR(r.history.back().independent_conservation_residual,-100.0,1e-12);
        EXPECT_NEAR(r.history.back().independent_conservation_normalized,1.0,1e-12);
    });

#ifndef CFDX_THERMAL_CONTRACTS_BASELINE
    run_case("production_temperature_admissibility_is_separate_from_conservation", [] {
        const auto mesh=make_1d_chain(1);
        const auto g=build_fv_geometry(mesh);
        auto setup=thermal_case();
        for (auto& bc : setup.boundary_conditions) {
            bc.type=cfdx::io::BCType::WALL;
            bc.value_type=cfdx::io::BCValueType::WALL_THERMAL;
            bc.wall_temperature=300.0;
            bc.raw_params.clear();
        }
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> T(1,"T","K",1);
        // Six unit-area walls at distance .5: 12*T = 12*300 + Su*V.
        // Do not clip the mathematical solution: report it and fail the
        // physical stage if that absolute temperature is not positive.
        for (const int power : {-4800,-3600,-3588}) {
            setup.materials[0].extra_properties["volumetric_power"]=std::to_string(power);
            T(0)=300.0;
            const auto r=run_production_energy_stage(mesh,g,setup,phi,T);
            const double exact=300.0+power/12.0;
            std::cout << "    admissibility oracle: Su=" << power << " T=" << T(0)
                      << " linear=" << r.linear_solver_converged
                      << " global=" << r.global_energy_balance_passed
                      << " local=" << r.local_energy_balance_passed
                      << " residual=" << r.independent_conservation_residual
                      << " converged=" << r.converged << '\n';
            EXPECT_NEAR(T(0),exact,1e-12);
            EXPECT_NEAR(r.min_temperature,exact,1e-12);
            EXPECT_TRUE(r.linear_solver_converged);
            EXPECT_TRUE(r.global_energy_balance_passed && r.local_energy_balance_passed);
            EXPECT_NEAR(r.independent_conservation_residual,0.0,1e-12);
            EXPECT_TRUE(r.converged==(exact>0.0));
#ifndef CFDX_THERMAL_ADMISSIBILITY_BASELINE
            EXPECT_TRUE(r.temperature_admissibility_passed==(exact>0.0));
            EXPECT_TRUE(r.nonpositive_temperature==(exact>0.0 ? 0u : 1u));
#endif
        }
        // The lower-level mathematical solver still supports relative or
        // manufactured temperature fields, including non-positive values.
        auto controls=energy_controls();
        ScalarBoundaryConditions bcs;
        for (const auto* name : {"inlet","outlet","wall"})
            bcs[name]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        Field<double,Location::CELL> source(1,"Su","W/m3",1); source(0)=-4800.0;
        T(0)=300.0;
        const auto kernel=solve_energy(mesh,g,phi,T,source,controls,bcs);
        EXPECT_TRUE(kernel.converged);
        EXPECT_NEAR(T(0),-100.0,1e-12);
    });

    run_case("independent_source_audit_detects_perturbation_and_rejects_invalid_sp", [] {
        const auto mesh=make_1d_chain(1);
        const auto g=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(1,"T","K",1),Su(1,"Su","W/m3",1),Sp(1,"Sp","W/m3/K",1);
        phi.fill(0.0); T(0)=310.0; Su(0)=600.0; Sp(0)=-2.0;
        const auto c=energy_controls();
        const auto audit=reconstruct_energy_balance(mesh,g,phi,T,T,Su,c,{},nullptr,&Sp);
        EXPECT_NEAR(audit.source_integral,-20.0,1e-12);
        EXPECT_NEAR(audit.residual,20.0,1e-12);
        EXPECT_NEAR(audit.max_cell_residual,20.0,1e-12);
        EXPECT_NEAR(audit.normalized_residual,1.0,1e-12);
        EXPECT_NEAR(energy_balance_relative(mesh,g,phi,T,T,Su,c,{},nullptr,nullptr,&Sp),1.0,1e-12);
        for (const double invalid : {1.0, std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity()}) {
            Sp(0)=invalid;
            EXPECT_THROW_WITH(reconstruct_energy_balance(mesh,g,phi,T,T,Su,c,{},nullptr,&Sp),
                              std::invalid_argument,"implicit source");
        }
        Field<double,Location::CELL> short_sp;
        EXPECT_THROW_WITH(reconstruct_energy_balance(mesh,g,phi,T,T,Su,c,{},nullptr,&short_sp),
                          std::invalid_argument,"implicit source field");
    });

    run_case("two_material_conduction_matches_series_resistance_and_both_audits", [] {
        const auto mesh=make_1d_chain(2);
        const auto g=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(2,"T","K",1),Su(2,"Su","W/m3",1);
        phi.fill(0.0); T.fill(300.0); Su.fill(0.0);
        const auto old=T;
        const std::vector<double> k{1.0,4.0};
        auto c=energy_controls(); c.cp=17.0;
        ScalarBoundaryConditions bcs;
        bcs["inlet"]={ScalarBoundaryType::FIXED_VALUE,400.0,0.0};
        bcs["outlet"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        const auto r=solve_energy(mesh,g,phi,T,Su,c,bcs,nullptr,nullptr,nullptr,&k);
        EXPECT_TRUE(r.converged);
        // R=0.5/1+0.5/4, q=100/R=160; cell centres lie at x=.25,.75.
        EXPECT_NEAR(T(0),360.0,1e-10);
        EXPECT_NEAR(T(1),310.0,1e-10);
        const auto flux=reconstruct_energy_face_flux(mesh,g,phi,T,c,bcs,nullptr,&k);
        EXPECT_NEAR(flux(0),-160.0,1e-10);
        EXPECT_NEAR(flux(1),160.0,1e-10);
        const auto audit=reconstruct_energy_balance(mesh,g,phi,T,old,Su,c,bcs,nullptr,nullptr,&k);
        EXPECT_NEAR(audit.residual,0.0,1e-9);
        EXPECT_TRUE(audit.max_cell_residual<1e-9);
        EXPECT_TRUE(energy_balance_relative(mesh,g,phi,T,old,Su,c,bcs,nullptr,&k)<1e-9);
        T(1)+=1.0;
        const auto disturbed=reconstruct_energy_balance(mesh,g,phi,T,old,Su,c,bcs,nullptr,nullptr,&k);
        EXPECT_TRUE(disturbed.max_cell_residual>10.0);
    });

    run_case("turbulent_conductivity_is_laminar_limit_and_prandtl_scaled", [] {
        const auto ec=build_production_energy_case(thermal_case());
        Field<double,Location::CELL> nu(2,"nu_t","m2/s",1);
        nu(0)=0.0; nu(1)=0.15;
        const auto k=build_effective_thermal_conductivity(ec,nu,1.0);
        EXPECT_NEAR(k[0],1.0,0.0);
        EXPECT_NEAR(k[1],4.0,1e-14);
        const auto k_half_pr=build_effective_thermal_conductivity(ec,nu,0.5);
        EXPECT_NEAR(k_half_pr[1],7.0,1e-14);
        for (const double invalid : {0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity()}) {
            EXPECT_THROW_WITH(build_effective_thermal_conductivity(ec,nu,invalid),
                              std::invalid_argument,"turbulent Prandtl");
        }
        for (const double invalid : {-1.0,std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity()}) {
            nu(1)=invalid;
            EXPECT_THROW_WITH(build_effective_thermal_conductivity(ec,nu,1.0),
                              std::invalid_argument,"turbulent viscosity");
        }
        nu(1)=std::numeric_limits<double>::max();
        EXPECT_THROW_WITH(build_effective_thermal_conductivity(ec,nu,1.0),
                          std::invalid_argument,"effective conductivity");
    });

    run_case("production_flux_bc_uses_each_owner_conductivity", [] {
        const auto mesh=make_1d_chain(2);
        const auto g=build_fv_geometry(mesh);
        const auto setup=thermal_case();
        const auto ec=build_production_energy_case(setup);
        const std::vector<double> k{1.0,4.0};
        const auto bcs=build_production_temperature_bcs(mesh,setup,ec);
        const auto faces=build_production_temperature_face_conditions(mesh,bcs,ec,k);
        for (const auto f : mesh.boundary().patch(2).face_ids) {
            const auto owner=mesh.ownership().owner(f);
            EXPECT_TRUE(faces.has(f));
            EXPECT_NEAR(faces.conditions[f].value,10.0/k[owner],1e-14);
        }
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(2,"T","K",1);
        phi.fill(0.0); T.fill(300.0);
        const auto r=run_production_energy_stage(mesh,g,setup,phi,T,{},&k);
        EXPECT_TRUE(r.converged);
        // Each cell receives 20 W from four side faces. The left boundary
        // removes 40 W: T0=300+40/4, T1=T0+20/(1/(.25/1+.25/4)).
        EXPECT_NEAR(T(0),310.0,1e-10);
        EXPECT_NEAR(T(1),316.25,1e-10);
        EXPECT_NEAR(r.patch_flux_out.at("wall"),-40.0,1e-10);
        EXPECT_NEAR(r.patch_flux_out.at("inlet"),40.0,1e-9);
        EXPECT_TRUE(r.global_energy_balance_passed && r.local_energy_balance_passed);
    });

    run_case("variable_conductivity_rejected_before_temperature_mutation", [] {
        const auto mesh=make_1d_chain(2);
        const auto g=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> T(2,"T","K",1),Su(2,"Su","W/m3",1);
        phi.fill(0.0); T.fill(300.0); Su.fill(0.0);
        auto c=energy_controls(); c.dt=1.0;
        for (const std::vector<double>& invalid : {
                 std::vector<double>{1.0},std::vector<double>{1.0,-1.0},
                 std::vector<double>{1.0,std::numeric_limits<double>::quiet_NaN()},
                 std::vector<double>{1.0,std::numeric_limits<double>::infinity()}}) {
            EXPECT_THROW_WITH(solve_energy(mesh,g,phi,T,Su,c,{},nullptr,nullptr,nullptr,&invalid),
                              std::invalid_argument,"conductivity");
            EXPECT_NEAR(T(0),300.0,0.0); EXPECT_NEAR(T(1),300.0,0.0);
        }
        const std::vector<double> zero_boundary_k{1.0,0.0};
        EXPECT_THROW_WITH(run_production_energy_stage(mesh,g,thermal_case(),phi,T,{},&zero_boundary_k),
                          std::invalid_argument,"positive");
        EXPECT_NEAR(T(0),300.0,0.0); EXPECT_NEAR(T(1),300.0,0.0);
        // The flux audit also requires a scalar temperature. Reject an
        // incompatible input before a successful solve can change component 0.
        Field<double,Location::CELL> vector_T(2,"T","K",2); vector_T.fill(300.0);
        Su.fill(10.0);
        EXPECT_THROW_WITH(solve_energy(mesh,g,phi,vector_T,Su,c),
                          std::invalid_argument,"energy field size");
        for (std::size_t cell=0;cell<2;++cell)
            for (std::size_t d=0;d<2;++d) EXPECT_NEAR(vector_T(cell,d),300.0,0.0);
    });

    run_case("accepted_sst_state_drives_energy_with_authoritative_flux", [] {
        SSTChannelCase carrier;
        const auto flow=solve_steady_incompressible(
            carrier.mesh,carrier.U,carrier.p,carrier.ubc,carrier.pbc,carrier.controls);
        // One accepted coupled iteration proves this software seam. It is
        // deliberately not a claim of a fully converged turbulent channel.
        EXPECT_TRUE(flow.history.size()==1);
        EXPECT_TRUE(carrier.state.k(0)<0.01 && carrier.state.k(0)>carrier.state.k(2));
        EXPECT_NEAR(flow.authoritative_mass_flux(0),-1.2,1e-12);
        EXPECT_NEAR(flow.authoritative_mass_flux(1),1.2,1e-12);
        auto setup=thermal_case();
        setup.materials[0].density=1.2;
        setup.materials[0].dynamic_viscosity=0.012;
        setup.materials[0].extra_properties["volumetric_power"]="30";
        setup.boundary_conditions[0].type=cfdx::io::BCType::INLET;
        setup.boundary_conditions[0].value_type=cfdx::io::BCValueType::FIXED;
        setup.boundary_conditions[0].wall_temperature=0.0;
        setup.boundary_conditions[0].temperature=300.0;
        setup.boundary_conditions[1].type=cfdx::io::BCType::OUTLET;
        setup.boundary_conditions[2].raw_params.clear();
        const auto ec=build_production_energy_case(setup);
        const auto k=build_effective_thermal_conductivity(ec,carrier.state.nu_t,0.9);
        for (std::size_t cell=0;cell<3;++cell) {
            EXPECT_NEAR(k[cell],1.0+1.2*10.0*carrier.state.nu_t(cell)/0.9,1e-14);
            EXPECT_TRUE(k[cell]>1.0);
        }
        Field<double,Location::CELL> T(3,"T","K",1); T.fill(300.0);
        auto laminar=T;
        const auto g=build_fv_geometry(carrier.mesh);
        const auto r=run_production_energy_stage(carrier.mesh,g,setup,
            flow.authoritative_mass_flux,T,{},&k);
        const auto r_laminar=run_production_energy_stage(carrier.mesh,g,setup,
            flow.authoritative_mass_flux,laminar);
        EXPECT_TRUE(r.converged && r_laminar.converged);
        // Independent 3x3 upwind/series-resistance system: dx=1/3,
        // m_dot*cp=12 W/K and Su*V=10 W in every cell. No solver assembly
        // or energy face-reconstruction helper is used by this oracle.
        const double d_left=6.0*k[0];
        const double d01=1.0/(1.0/(6.0*k[0])+1.0/(6.0*k[1]));
        const double d12=1.0/(1.0/(6.0*k[1])+1.0/(6.0*k[2]));
        const double lower[3]={0.0,-(12.0+d01),-(12.0+d12)};
        double diagonal[3]={d_left+d01+12.0,d01+d12+12.0,d12+12.0};
        const double upper[3]={-d01,-d12,0.0};
        double exact[3]={10.0+(d_left+12.0)*300.0,10.0,10.0};
        for (std::size_t cell=1;cell<3;++cell) {
            const double factor=lower[cell]/diagonal[cell-1];
            diagonal[cell]-=factor*upper[cell-1];
            exact[cell]-=factor*exact[cell-1];
        }
        exact[2]/=diagonal[2];
        for (int cell=1;cell>=0;--cell)
            exact[cell]=(exact[cell]-upper[cell]*exact[cell+1])/diagonal[cell];
        double max_effect=0.0;
        for (std::size_t cell=0;cell<3;++cell) {
            EXPECT_NEAR(T(cell),exact[cell],1e-9);
            max_effect=std::max(max_effect,std::abs(T(cell)-laminar(cell)));
        }
        EXPECT_TRUE(max_effect>1e-4);
        EXPECT_NEAR(r.boundary_flux_out,30.0,1e-8);
        EXPECT_TRUE(r.local_energy_balance_passed && r.global_energy_balance_passed);
        std::cout << "    SST-energy seam: max thermal effect=" << max_effect
                  << " K; balance=" << r.independent_conservation_residual << " W\n";
    });

    run_case("restored_temperature_drives_the_first_buoyancy_source", [] {
        SSTChannelCase carrier;
        carrier.controls.sst_runtime=nullptr;
        auto setup=thermal_case();
        setup.energy_model="buoyant";
        setup.ref_temperature=300.0;
        setup.gravity_vector="0 -9.81 0";
        setup.materials[0].density=1.2;
        setup.materials[0].dynamic_viscosity=0.012;
        setup.materials[0].extra_properties["thermal_expansion_coefficient"]="0.001";
        Field<double,Location::CELL> loaded_T(3,"T","K",1);
        for (std::size_t cell=0;cell<3;++cell) loaded_T(cell)=310.0+10.0*cell;
        auto expected_U=carrier.U; auto expected_p=carrier.p;
        auto reference_controls=carrier.controls;
        // Build the oracle directly, independently of the thermal seeding
        // and of build_cell_buoyancy_acceleration.
        reference_controls.cell_body_force={
            {0.0,0.001*10.0*9.81,0.0},
            {0.0,0.001*20.0*9.81,0.0},
            {0.0,0.001*30.0*9.81,0.0}};
        const auto reference=solve_steady_incompressible(carrier.mesh,expected_U,expected_p,
            carrier.ubc,carrier.pbc,reference_controls);
        const auto result=run_production_buoyant_outer_loop(carrier.mesh,carrier.ubc,carrier.pbc,
            carrier.controls,setup,carrier.U,carrier.p,{},&loaded_T);
        EXPECT_TRUE(result.outer_iterations==1 && reference.history.size()==1);
        double vertical_motion=0.0;
        for (std::size_t cell=0;cell<3;++cell) {
            for (std::size_t d=0;d<3;++d)
                EXPECT_NEAR(carrier.U(cell,d),expected_U(cell,d),1e-13);
            EXPECT_NEAR(carrier.p(cell),expected_p(cell),1e-13);
            vertical_motion=std::max(vertical_motion,std::abs(carrier.U(cell,1)));
            EXPECT_NEAR(loaded_T(cell),310.0+10.0*cell,0.0);
        }
        EXPECT_TRUE(vertical_motion>1e-4);
        std::cout << "    restored buoyancy: first-iterate max Uy=" << vertical_motion << " m/s\n";
    });

    run_case("invalid_restored_temperature_is_rejected_before_carrier_mutation", [] {
        SSTChannelCase carrier;
        carrier.controls.sst_runtime=nullptr;
        auto setup=thermal_case();
        setup.energy_model="buoyant"; setup.ref_temperature=300.0;
        setup.gravity_vector="0 -9.81 0";
        setup.materials[0].extra_properties["thermal_expansion_coefficient"]="0.001";
        Field<double,Location::CELL> loaded_T(3,"T","K",1); loaded_T.fill(300.0);
        std::size_t predictor_calls=0;
        carrier.controls.momentum_predictor_callback=[&](std::size_t,const auto&,const auto&,const auto&) {
            ++predictor_calls;
        };
        for (const double invalid : {0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity()}) {
            loaded_T(2)=invalid;
            EXPECT_THROW_WITH(run_production_buoyant_outer_loop(carrier.mesh,carrier.ubc,carrier.pbc,
                carrier.controls,setup,carrier.U,carrier.p,{},&loaded_T),
                std::invalid_argument,"initial temperature field must be finite and positive");
        }
        Field<double,Location::CELL> short_T(2,"T","K",1); short_T.fill(300.0);
        Field<double,Location::CELL> vector_T(3,"T","K",2); vector_T.fill(300.0);
        for (const auto* invalid : {&short_T,&vector_T})
            EXPECT_THROW_WITH(run_production_buoyant_outer_loop(carrier.mesh,carrier.ubc,carrier.pbc,
                carrier.controls,setup,carrier.U,carrier.p,{},invalid),
                std::invalid_argument,"initial temperature field must be a scalar matching the mesh");
        EXPECT_TRUE(predictor_calls==0);
        for (std::size_t cell=0;cell<3;++cell) {
            for (std::size_t d=0;d<3;++d) EXPECT_NEAR(carrier.U(cell,d),-7.0,0.0);
            EXPECT_NEAR(carrier.p(cell),4.0,0.0);
        }
    });
#endif
    return run_all();
}
