#include "sst_channel_fixture.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_shear_chain(std::size_t n)
{
    auto mesh=make_1d_chain(n);
    const auto faces=mesh.boundary().patch(2).face_ids;
    mesh.boundary().remove_patch(2);
    const std::array<std::string,4> names{"y0","y1","z0","z1"};
    for (std::size_t direction=0; direction<4; ++direction) {
        Patch patch; patch.name=names[direction]; patch.type=PatchType::WALL;
        for (std::size_t cell=0; cell<n; ++cell)
            patch.face_ids.push_back(faces[direction*n+cell]);
        mesh.boundary().add_patch(patch);
    }
    return mesh;
}

static std::array<ScalarBoundaryConditions,3> shear_scalar_bcs(const Mesh& mesh)
{
    std::array<ScalarBoundaryConditions,3> bcs;
    for (std::size_t component=0; component<3; ++component)
        for (const auto& patch : mesh.boundary().patches()) {
            const double value=component==0 ? (patch.name=="y0" ? 0.0 :
                                               (patch.name=="y1" ? 1.0 : 0.5)) : 0.0;
            bcs[component][patch.name]={ScalarBoundaryType::FIXED_VALUE,value,0.0};
        }
    return bcs;
}

int main()
{
    run_case("sst_wall_omega_requires_every_mesh_face_before_carrier_mutation", [] {
        SSTChannelCase c;
        const auto geometry = build_fv_geometry(c.mesh);
        const SSTRuntimeBoundaryContract contract{
            c.runtime.k_boundary_conditions, c.runtime.omega_boundary_conditions, {"wall"}};
        require_supported_sst_runtime_boundaries(contract, c.mesh, geometry, &c.runtime.omega_face_values);
        const auto wall_faces = c.mesh.boundary().patch(2).face_ids;
        EXPECT_TRUE(wall_faces.size() == 12);
        const auto good = c.runtime.omega_face_values;
        for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity(), 0.0, -1.0}) {
            c.runtime.omega_face_values = good;
            c.runtime.omega_face_values.values.at("wall")[wall_faces.back()] = invalid;
            EXPECT_THROW_WITH(solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls),
                              std::invalid_argument, "missing or invalid value");
            for (std::size_t cell=0; cell<c.mesh.n_cells(); ++cell) {
                for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(c.U(cell,d),-7.0,0.0);
                EXPECT_NEAR(c.p(cell),4.0,0.0);
                EXPECT_NEAR(c.state.k(cell),0.01,0.0);
                EXPECT_NEAR(c.state.omega(cell),2.0,0.0);
                EXPECT_NEAR(c.state.nu_t(cell),9.0,0.0);
                EXPECT_NEAR(c.state.wall_distance(cell),17.0,0.0);
                EXPECT_NEAR(c.state.strain_rate(cell),11.0,0.0);
            }
        }
        auto short_values = good;
        short_values.values.at("wall").resize(1);
        EXPECT_THROW(require_supported_sst_runtime_boundaries(contract,c.mesh,geometry,&short_values),
                     std::invalid_argument);
        auto wrong_geometry = geometry;
        wrong_geometry.face_patch[wall_faces.front()] = 0;
        EXPECT_THROW(require_supported_sst_runtime_boundaries(contract,c.mesh,wrong_geometry,&good),
                     std::invalid_argument);
        auto flat_geometry = geometry;
        const auto face = wall_faces.front();
        flat_geometry.face_centres[face] = geometry.cell_centres[c.mesh.ownership().owner(face)];
        EXPECT_THROW_WITH(require_supported_sst_runtime_boundaries(contract,c.mesh,flat_geometry,&good),
                          std::invalid_argument,"positive geometric normal distance");
    });

    run_case("sst_cross_diffusion_uses_real_scalar_and_face_boundaries", [] {
        const auto mesh = make_1d_chain(3);
        const auto geometry = build_fv_geometry(mesh);
        Field<double,Location::CELL> k(3,"k","",1), omega(3,"omega","",1);
        for (std::size_t cell=0; cell<3; ++cell) {
            const auto x = geometry.cell_centres[cell];
            k(cell)=2+x.x;
            omega(cell)=3+2*x.x+3*x.y;
        }
        ScalarBoundaryConditions kbc{
            {"inlet",{ScalarBoundaryType::FIXED_VALUE,2,0}},
            {"outlet",{ScalarBoundaryType::FIXED_VALUE,3,0}}};
        ScalarBoundaryConditions wbc{
            {"inlet",{ScalarBoundaryType::FIXED_VALUE,4.5,0}},
            {"outlet",{ScalarBoundaryType::FIXED_VALUE,6.5,0}},
            {"wall",{ScalarBoundaryType::FIXED_VALUE,999,0}}};
        ScalarBoundaryFaceValues faces;
        faces.values["wall"].assign(mesh.n_faces(),std::numeric_limits<double>::quiet_NaN());
        for (const auto face : mesh.boundary().patch(2).face_ids) {
            const auto x=geometry.face_centres[face];
            faces.values["wall"][face]=3+2*x.x+3*x.y;
        }
        const auto contraction=sst_boundary_gradient_dot(k,omega,mesh,geometry,kbc,wbc,&faces);
        const auto gradw=scalar_gradient_with_transport_boundaries(omega,mesh,geometry,wbc,&faces);
        for (std::size_t cell=0; cell<3; ++cell) {
            EXPECT_NEAR(contraction[cell],2.0,2e-14);
            EXPECT_NEAR(gradw(cell,0),2.0,2e-14);
            EXPECT_NEAR(gradw(cell,1),3.0,2e-14);
            EXPECT_NEAR(gradw(cell,2),0.0,2e-14);
        }
        const auto old_contraction=cell_gradient_dot(k,omega,mesh);
        EXPECT_NEAR(old_contraction.front(),0.5,2e-14);
        EXPECT_NEAR(old_contraction.back(),0.5,2e-14);
        kbc["inlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,-1};
        kbc["outlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,1};
        const auto gradk=scalar_gradient_with_transport_boundaries(k,mesh,geometry,kbc);
        for (std::size_t cell=0; cell<3; ++cell) EXPECT_NEAR(gradk(cell,0),1.0,2e-14);
    });

    run_case("sst_residuals_and_convergence_describe_bounded_returned_fields", [] {
        const auto mesh=make_1d_chain(1);
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"phi","",1); flux.fill(0);
        Field<double,Location::CELL> k(1,"k","",1),w(1,"omega","",1),S(1,"S","",1),F(1,"F","",1);
        k.fill(1); w.fill(1); S.fill(0); F.fill(1);
        TurbulenceTransportControls controls;
        controls.model=TurbulenceModel::SST;
        controls.k_min=1; controls.omega_min=1;
        const auto result=solve_sst_transport(mesh,geometry,flux,k,w,S,F,F,controls,{}, {},3,1e-8);
        EXPECT_FALSE(result.converged);
        EXPECT_NEAR(k(0),1.0,0.0); EXPECT_NEAR(w(0),1.0,0.0);
        EXPECT_NEAR(result.k_residual,controls.density*controls.beta_star,1e-14);
        EXPECT_NEAR(result.second_residual,controls.density*controls.beta1,1e-14);
        EXPECT_TRUE(result.k_floor_events==3 && result.second_floor_events==3);
        std::cout << "    bounded-state residuals: k=" << result.k_residual
                  << " omega=" << result.second_residual << " converged=" << result.converged << '\n';
        ScalarEquation finite_equation{SparseMatrix(1,1),Vector(1,2.0),{1.0},0.0};
        finite_equation.matrix.push_back(0,0,1.0);
        finite_equation.matrix.finalize();
        EXPECT_NEAR(sst_linear_tolerance_for_solution_defect(finite_equation,1e-8),5e-9,1e-24);
        EXPECT_THROW_WITH(sst_linear_tolerance_for_solution_defect(
                              finite_equation,std::numeric_limits<double>::denorm_min()),
                          std::domain_error,"cannot be represented");
        finite_equation.rhs(0)=std::numeric_limits<double>::infinity();
        EXPECT_THROW(sst_linear_tolerance_for_solution_defect(finite_equation,1e-8),std::invalid_argument);
        ScalarEquation zero_row{SparseMatrix(1,1),Vector(1,1.0),{0.0},0.0};
        zero_row.matrix.push_back(0,0,0.0);
        zero_row.matrix.finalize();
        EXPECT_THROW_WITH(sst_linear_tolerance_for_solution_defect(zero_row,1e-8),
                          std::invalid_argument,"zero or non-finite matrix row norm");
    });

    run_case("sst_multicell_first_predictor_uses_refreshed_viscosity_and_nonzero_flux", [] {
        SSTChannelCase c;
        std::size_t calls=0;
        c.controls.momentum_predictor_callback = [&](std::size_t iteration,const auto& mu,const auto& U,const auto& flux) {
            ++calls;
            EXPECT_TRUE(iteration==1);
            for (std::size_t cell=0; cell<c.mesh.n_cells(); ++cell) {
                EXPECT_NEAR(U(cell,0),1.0,1e-14);
                EXPECT_NEAR(c.state.strain_rate(cell),0.0,1e-14);
                EXPECT_NEAR(c.state.wall_distance(cell),0.5,1e-14);
                EXPECT_NEAR(mu[cell],1.2*(0.01+0.01/2.0),1e-14);
            }
            EXPECT_NEAR(flux(0),-1.2,1e-14);
            EXPECT_NEAR(flux(1),1.2,1e-14);
            for (std::size_t face=0; face<c.mesh.n_faces(); ++face)
                if (c.mesh.ownership().neighbour(face)>=0) EXPECT_NEAR(flux(face),1.2,1e-14);
        };
        const auto result=solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(calls==1 && result.history.size()==1);
        EXPECT_TRUE(result.authoritative_mass_flux.size()==c.mesh.n_faces());
        EXPECT_NEAR(result.authoritative_mass_flux(0),-1.2,1e-12);
        EXPECT_NEAR(result.authoritative_mass_flux(1),1.2,1e-12);
        EXPECT_TRUE(c.state.k(0)>c.state.k(2));
        EXPECT_TRUE(c.state.k(0)<0.01);
        std::cout << "    nonzero-flux SST: k=[" << c.state.k(0) << ',' << c.state.k(1) << ',' << c.state.k(2)
                  << "] omega=[" << c.state.omega(0) << ',' << c.state.omega(1) << ',' << c.state.omega(2) << "]\n";
    });

    run_case("sst_multicell_retry_restores_every_field_and_authoritative_flux", [] {
        SSTChannelCase c;
        c.controls.nonlinear_retry.max_retries=1;
        std::size_t predictor_calls=0,output_calls=0;
        Field<double,Location::CELL> first_U;
        Field<double,Location::FACE> first_flux;
        SSTProductionState first_state(3);
        std::vector<double> first_mu;
        c.controls.momentum_predictor_callback = [&](std::size_t,const auto& mu,const auto& U,const auto& flux) {
            ++predictor_calls;
            if (predictor_calls==1) { first_U=U; first_flux=flux; first_state=c.state; first_mu=mu; }
            else {
                for (std::size_t cell=0; cell<3; ++cell) {
                    for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(U(cell,d),first_U(cell,d),0.0);
                    EXPECT_NEAR(c.state.k(cell),first_state.k(cell),0.0);
                    EXPECT_NEAR(c.state.omega(cell),first_state.omega(cell),0.0);
                    EXPECT_NEAR(c.state.nu_t(cell),first_state.nu_t(cell),0.0);
                    EXPECT_NEAR(c.state.wall_distance(cell),first_state.wall_distance(cell),0.0);
                    EXPECT_NEAR(c.state.strain_rate(cell),first_state.strain_rate(cell),0.0);
                    EXPECT_NEAR(mu[cell],first_mu[cell],0.0);
                }
                for (std::size_t face=0; face<flux.size(); ++face) EXPECT_NEAR(flux(face),first_flux(face),0.0);
            }
        };
        c.controls.iteration_output_callback = [&](std::size_t,double,const auto&,const auto&,const auto&) {
            if (++output_calls==1) {
                EXPECT_TRUE(c.state.k(0)!=first_state.k(0));
                throw NonlinearRetryableFailure("test rejection after SST commit with nonzero flux");
            }
            return false;
        };
        const auto result=solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(predictor_calls==2 && output_calls==2 && result.history.size()==1);
        SSTChannelCase reference;
        reference.controls.coupling.alpha_u=result.history[0].effective_alpha_u;
        reference.controls.coupling.alpha_p=result.history[0].effective_alpha_p;
        const auto clean=solve_steady_incompressible(reference.mesh,reference.U,reference.p,
                                                   reference.ubc,reference.pbc,reference.controls);
        for (std::size_t cell=0; cell<3; ++cell) {
            for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(c.U(cell,d),reference.U(cell,d),1e-12);
            EXPECT_NEAR(c.p(cell),reference.p(cell),1e-12);
            EXPECT_NEAR(c.state.k(cell),reference.state.k(cell),1e-12);
            EXPECT_NEAR(c.state.omega(cell),reference.state.omega(cell),1e-12);
            EXPECT_NEAR(c.state.nu_t(cell),reference.state.nu_t(cell),1e-12);
        }
        for (std::size_t face=0; face<c.mesh.n_faces(); ++face)
            EXPECT_NEAR(result.authoritative_mass_flux(face),clean.authoritative_mass_flux(face),1e-12);
    });

    run_case("sst_multicell_dat_rehydrates_primary_fields_before_first_predictor", [] {
        SSTChannelCase saved;
        const auto accepted=solve_steady_incompressible(
            saved.mesh,saved.U,saved.p,saved.ubc,saved.pbc,saved.controls);
        EXPECT_NEAR(accepted.authoritative_mass_flux(0),-1.2,1e-12);
        const auto path=unique_temp_path("_sst_multicell_restore.dat");
        cfdx::io::DatRestartFields fields;
        fields.k=&saved.state.k;
        fields.omega=&saved.state.omega;
        cfdx::io::write_dat_restart_fields(path.string(),saved.mesh,saved.U,saved.p,fields,1,0.0);
        SSTChannelCase restored;
        const std::map<std::string,std::array<double,3>> fixed_velocity{
            {"inlet",{1,0,0}},{"wall",{1,0,0}}};
        const auto metadata=restore_sst_production_state_from_dat(
            path.string(),restored.mesh,restored.U,restored.p,restored.state,
            {0.5,0.5,0.5},restored.runtime.transport,fixed_velocity);
        EXPECT_TRUE(metadata.iteration==1);
        for (std::size_t cell=0; cell<3; ++cell) {
            for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(restored.U(cell,d),saved.U(cell,d),0.0);
            EXPECT_NEAR(restored.p(cell),saved.p(cell),0.0);
            EXPECT_NEAR(restored.state.k(cell),saved.state.k(cell),0.0);
            EXPECT_NEAR(restored.state.omega(cell),saved.state.omega(cell),0.0);
            EXPECT_NEAR(restored.state.nu_t(cell),saved.state.nu_t(cell),1e-15);
        }
        restored.controls.initialization.mode=InitializationMode::Provided;
        std::size_t predictor_calls=0;
        restored.controls.momentum_predictor_callback=[&](std::size_t,const auto& mu,const auto& U,const auto&) {
            ++predictor_calls;
            for (std::size_t cell=0; cell<3; ++cell) {
                EXPECT_NEAR(mu[cell],1.2*(0.01+saved.state.nu_t(cell)),1e-15);
                for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(U(cell,d),saved.U(cell,d),0.0);
            }
        };
        solve_steady_incompressible(restored.mesh,restored.U,restored.p,
                                   restored.ubc,restored.pbc,restored.controls);
        EXPECT_TRUE(predictor_calls==1);
        std::filesystem::remove(path);
        // DAT v4 does not store authoritative face flux. This test certifies
        // primary-field hydration and derived viscosity only, not equivalence
        // of a general interrupted coupled trajectory.
    });

    run_case("variable_viscosity_transpose_stress_mms_converges_and_is_conservative", [] {
        // Divergence-free U=(y,0,0), mu=1+x^2(1-x)^2. The scalar diffusion
        // alone gives zero, whereas div(mu*grad(U)^T)=(0,mu'(x),0).
        // mu'(0)=mu'(1)=0 makes the declared owner boundary-mu convention
        // consistent to second order on this orthogonal family.
        double previous_error=0.0;
        for (const std::size_t n : {16u,32u,64u}) {
            const auto mesh=make_shear_chain(n);
            const auto geometry=build_fv_geometry(mesh);
            Field<double,Location::CELL> U(n,"U","m/s",3); U.fill(0.0);
            std::vector<double> mu(n);
            for (std::size_t cell=0; cell<n; ++cell) {
                U(cell,0)=geometry.cell_centres[cell].y;
                const double x=geometry.cell_centres[cell].x;
                mu[cell]=1+x*x*(1-x)*(1-x);
            }
            const auto flux=transpose_viscous_stress_flux(mesh,geometry,U,mu,shear_scalar_bcs(mesh));
            const auto force=viscous_face_flux_divergence(mesh,geometry,flux);
            double error_sq=0.0,omitted_sq=0.0,volume_sum=0.0,boundary_sum=0.0;
            for (std::size_t cell=0; cell<n; ++cell) {
                const double x=geometry.cell_centres[cell].x;
                const double exact=2*x*(1-x)*(1-2*x);
                error_sq+=std::pow(force(cell,1)-exact,2)*geometry.cell_volumes[cell];
                omitted_sq+=exact*exact*geometry.cell_volumes[cell];
                volume_sum+=force(cell,1)*geometry.cell_volumes[cell];
                EXPECT_NEAR(force(cell,0),0.0,2e-12);
                EXPECT_NEAR(force(cell,2),0.0,2e-12);
            }
            for (std::size_t face=0; face<mesh.n_faces(); ++face)
                if (mesh.ownership().neighbour(face)<0) boundary_sum+=flux(face,1);
            EXPECT_NEAR(volume_sum,boundary_sum,2e-13);
            const double error=std::sqrt(error_sq);
            EXPECT_TRUE(error<std::sqrt(omitted_sq)/20.0);
            double order=0.0;
            if (previous_error>0.0) {
                order=std::log(previous_error/error)/std::log(2.0);
                EXPECT_TRUE(order>1.8);
            }
            std::cout << "    transpose-stress MMS N=" << n << " L2=" << error
                      << " omitted=" << std::sqrt(omitted_sq) << " order=" << order << '\n';
            previous_error=error;
            std::fill(mu.begin(),mu.end(),1.0);
            const auto constant_flux=transpose_viscous_stress_flux(mesh,geometry,U,mu,shear_scalar_bcs(mesh));
            const auto constant_force=viscous_face_flux_divergence(mesh,geometry,constant_flux);
            for (std::size_t cell=0; cell<n; ++cell)
                for (std::size_t component=0; component<3; ++component)
                    EXPECT_NEAR(constant_force(cell,component),0.0,2e-12);
        }
    });

    run_case("variable_viscosity_production_predictor_consumes_transpose_stress", [] {
        const auto mesh=make_shear_chain(8);
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> U(8,"U","m/s",3),p(8,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        const auto scalar_bcs=shear_scalar_bcs(mesh);
        for (const auto& patch : mesh.boundary().patches()) {
            ubc[patch.name]={VelocityBoundaryCondition::Type::FIXED_VALUE,
                            {scalar_bcs[0].at(patch.name).value,0,0}};
            pbc[patch.name]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        }
        IncompressibleSolverControls controls;
        controls.kinematic_viscosity=1.0;
        controls.convergence.max_iterations=1;
        controls.linear_tolerance=1e-12;
        controls.cell_turbulent_viscosity.resize(8);
        for (std::size_t cell=0; cell<8; ++cell) {
            const double x=geometry.cell_centres[cell].x;
            U(cell,0)=0.5;
            controls.cell_turbulent_viscosity[cell]=x*x*(1-x)*(1-x);
        }
        const auto result=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
        EXPECT_TRUE(result.history.size()==1);
        // With the old scalar-only operator Uy stays exactly zero. The new
        // transverse motion must have the sign of mu'(x), not its opposite.
        EXPECT_TRUE(U(0,1)>1e-5);
        EXPECT_TRUE(U(7,1)<-1e-5);
        std::cout << "    transverse predictor Uy=" << U(0,1) << ',' << U(7,1) << '\n';
    });

    run_case("momentum_balance_uses_physical_source_sign_and_density_in_laminar_limit", [] {
        const auto mesh=make_1d_chain(1);
        Field<double,Location::CELL> U(1,"U","m/s",3),p(1,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ScalarBoundaryConditions pbc;
        for (const auto& patch : mesh.boundary().patches()) {
            ubc[patch.name]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
            pbc[patch.name]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        }
        IncompressibleSolverControls controls;
        controls.density=2.5;
        controls.kinematic_viscosity=0.4;
        controls.cell_turbulent_viscosity={0.0};
        controls.body_force={3,0,0};
        controls.convergence.max_iterations=1;
        controls.coupling.alpha_u=0.7;
        const auto result=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
        // Six unit faces at distance 1/2 give conductance 12*mu=12.
        // Outward diffusive flux is 12*Ux, source integral is -rho*f=-7.5.
        EXPECT_NEAR(U(0,0),0.7*7.5/12.0,1e-12);
        EXPECT_NEAR(result.history.at(0).momentum_conservation_residual[0],-2.25,1e-12);
        const auto cellwise_U=U;
        U.fill(0.0); p.fill(0.0);
        controls.cell_turbulent_viscosity.clear();
        const auto uniform=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
        EXPECT_NEAR(U(0,0),cellwise_U(0,0),1e-12);
        EXPECT_NEAR(uniform.history.at(0).momentum_conservation_residual[0],-2.25,1e-12);
    });
    return run_all();
}
