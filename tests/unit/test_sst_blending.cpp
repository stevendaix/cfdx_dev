#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/sst_strain_rate.h"
#include "cfdx/physics/sst_case_contract.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include <array>
#include "common/test_harness.h"
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace cfdx::testing;
using namespace cfdx::physics;

static cfdx::core::Mesh make_sst_wall_operator_test_mesh()
{
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{4,0,0},{4,2,0},{0,2,0},
        {0,0,2},{4,0,2},{4,2,2},{0,2,2}
    };
    for (std::size_t i = 0; i < 8; ++i)
        mesh.points().set(i, p[i][0], p[i][1], p[i][2]);
    const std::vector<std::vector<FaceIndex>> faces = {
        {0,4,7,3}, {1,2,6,5}, {0,1,5,4},
        {3,7,6,2}, {0,3,2,1}, {4,5,6,7}
    };
    for (const auto& face : faces)
        mesh.faces().push_face(face);
    mesh.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        mesh.ownership().set_owner(f, 0);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    const std::array<std::pair<const char*, PatchType>, 6> patches{{
        {"wall", PatchType::WALL}, {"inlet", PatchType::INLET},
        {"side_low", PatchType::SYMMETRY}, {"side_high", PatchType::SYMMETRY},
        {"bottom", PatchType::SYMMETRY}, {"top", PatchType::SYMMETRY}
    }};
    for (std::size_t f = 0; f < patches.size(); ++f) {
        Patch patch;
        patch.name = patches[f].first;
        patch.type = patches[f].second;
        patch.face_ids = {static_cast<FaceIndex>(f)};
        mesh.boundary().add_patch(patch);
    }
    return mesh;
}



static cfdx::core::Mesh make_sst_transaction_unit_cube()
{
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double points[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        mesh.points().set(i, points[i][0], points[i][1], points[i][2]);

    const std::vector<std::vector<FaceIndex>> faces = {
        {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
        {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
    };
    for (const auto& face : faces)
        mesh.faces().push_face(face);

    mesh.ownership().resize(faces.size());
    for (std::size_t face = 0; face < faces.size(); ++face) {
        mesh.ownership().set_owner(face, 0);
        mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});

    Patch x0;
    x0.name = "x0";
    x0.type = PatchType::WALL;
    x0.face_ids = {4};
    mesh.boundary().add_patch(x0);
    Patch x1;
    x1.name = "x1";
    x1.type = PatchType::WALL;
    x1.face_ids = {5};
    mesh.boundary().add_patch(x1);
    Patch walls;
    walls.name = "walls";
    walls.type = PatchType::WALL;
    walls.face_ids = {0,1,2,3};
    mesh.boundary().add_patch(walls);
    return mesh;
}

int main()
{
    run_case("v11_sst_post_damping_residual_is_reassembled_on_returned_fields", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_transaction_unit_cube();
        const auto geometry = build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"mdot","kg/s",1);
        flux.fill(0.0);
        Field<double,Location::CELL> k(mesh.n_cells(),"k","m2/s2",1),
            omega(mesh.n_cells(),"w","1/s",1),
            strain(mesh.n_cells(),"strain","1/s",1),
            wall(mesh.n_cells(),"y","m",1);
        k.fill(0.4); omega.fill(8.0); strain.fill(2.0); wall.fill(0.1);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.sst_picard_damping = 0.2;
        ScalarBoundaryConditions k_bcs, w_bcs;
        for (const char* patch : {"x0", "x1", "walls"}) {
            k_bcs[patch] = {ScalarBoundaryType::FIXED_VALUE, 0.01, 0.0};
            w_bcs[patch] = {ScalarBoundaryType::FIXED_VALUE, 2.0, 0.0};
        }
        const auto r = solve_sst_transport_dynamic_blending(
            mesh,geometry,flux,k,omega,strain,wall,controls,
            k_bcs,w_bcs,1,1e-12);
        const auto expected = sst_reassembled_returned_state_residual(
            mesh,geometry,flux,k,omega,strain,wall,controls,k_bcs,w_bcs,nullptr);
        EXPECT_TRUE(!r.converged);
        EXPECT_NEAR(r.k_residual,expected.k_raw_inf,1e-12);
        EXPECT_NEAR(r.second_residual,expected.omega_raw_inf,1e-12);
        EXPECT_TRUE(std::isfinite(r.k_residual) && std::isfinite(r.second_residual));
    });
    run_case("v15_extended_residual_detects_fp64_cancellation_on_same_matrix", [] {
        using namespace cfdx::core;
        ScalarEquation equation;
        equation.matrix=SparseMatrix(3,3);
        equation.matrix.push_back(0,0,1e16);
        equation.matrix.push_back(0,1,-1e16);
        equation.matrix.push_back(0,2,1.0);
        equation.matrix.push_back(1,1,1.0);
        equation.matrix.push_back(2,2,1.0);
        equation.matrix.finalize();
        equation.rhs=Vector(3,0.0);
        equation.rhs(0)=1.0;
        equation.rhs(1)=1.0;
        equation.rhs(2)=1.0;
        Vector x(3,1.0);
        const auto r=sst_local_residual_precision_audit(equation,x,1e-20);
        EXPECT_TRUE(r.native.maximum_relative>0.0);
        EXPECT_NEAR(r.compensated.maximum_relative,0.0,1e-25);
        EXPECT_TRUE(r.maximum_row_disagreement>0.0);
        EXPECT_TRUE(r.disagreement_cell==0);
        equation.rhs(0)=0.0;
        const auto actual=sst_local_residual_precision_audit(equation,x,1e-20);
        EXPECT_TRUE(actual.compensated.maximum_relative>0.0);
        EXPECT_NEAR(actual.compensated.maximum_relative,
                    1.0/(2e16+1.0),1e-28);
    });
    run_case("v15_wall_audit_distinguishes_face_dirichlet_and_first_cell_omega", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_wall_operator_test_mesh();
        const auto geom=build_fv_geometry(mesh);
        TurbulenceTransportControls controls;
        controls.model=TurbulenceModel::SST;
        controls.molecular_viscosity=1e-5;
        controls.sst_wall_omega_factor=6.0;
        const auto values=make_sst_wall_omega_face_values(
            mesh,geom,std::vector<std::string>{"wall"},controls);
        Field<double,Location::CELL> U(1,"velocity","m/s",3);
        U(0,0)=0.0; U(0,1)=2.0; U(0,2)=0.0;
        Field<double,Location::CELL> w(1,"omega","1/s",1);
        const double y=2.0;
        const double nasa=60.0*controls.molecular_viscosity/(controls.beta1*y*y);
        w(0)=nasa;
        const auto audit=audit_sst_wall_omega_discretization(
            mesh,geom,{"wall"},values,U,w,
            {{"wall",{0.0,0.0,0.0}}},controls);
        EXPECT_TRUE(audit.faces==1 && audit.velocity_proxy_faces==1);
        EXPECT_NEAR(audit.min_distance,y,1e-12);
        EXPECT_NEAR(audit.min_nasa_first_cell_value,nasa,1e-12);
        EXPECT_NEAR(audit.min_face_value,nasa/10.0,1e-12);
        EXPECT_NEAR(audit.min_owner_to_first_cell_ratio,1.0,1e-12);
        EXPECT_NEAR(audit.max_viscous_yplus_proxy,std::sqrt(y*2.0/1e-5),1e-9);
        EXPECT_NEAR(audit.max_tangential_geometry_ratio,0.0,1e-12);
        EXPECT_THROW(audit_sst_wall_omega_discretization(
            mesh,geom,{"missing"},values,U,w,{},controls),std::invalid_argument);
    });
    run_case("v15_scalar_relaxation_validates_domain_and_keeps_default", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        EXPECT_NEAR(c.sst_linear_relaxation,0.7,0.0);
        c.sst_linear_relaxation=1.0;
        validate_turbulence_controls(c);
        c.sst_linear_relaxation=0.0;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_linear_relaxation=1.1;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
    });
    run_case("v16_sst_source_jacobian_and_cross_diffusion_identity", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        c.density=1.2;
        c.sst_omega_newton_sink=true;
        const double k=0.31,omega=25.0,S=4.0,F1=0.35,F2=0.8,cross=-1.2;
        const auto q=sst_cell_sources(k,omega,S,F1,F2,cross,c);
        const double beta=F1*c.beta1+(1.0-F1)*c.beta2;
        const double gamma=F1*c.gamma1+(1.0-F1)*c.gamma2;
        const double limitedPk=std::min(c.density*q.nut*S*S,
            c.sst_production_limiter*c.beta_star*c.density*k*omega);
        const double cross_source=2.0*(1.0-F1)*c.density*c.sst_sigma_w2/
                                  omega*cross;
        const double omega_physical=gamma*limitedPk/q.nut -
            c.density*beta*omega*omega+cross_source;
        EXPECT_NEAR(q.sk+q.spk*k,limitedPk-c.density*c.beta_star*k*omega,1e-12);
        EXPECT_NEAR(q.sw+q.spw*omega,omega_physical,1e-11);
        // The Newton option linearizes only the quadratic destruction.
        const double h=1e-4;
        const auto sink=[&](double w) {return -c.density*beta*w*w;};
        const double fd=(sink(omega+h)-sink(omega-h))/(2.0*h);
        EXPECT_NEAR(fd,-2.0*c.density*beta*omega,1e-8);
        // With negative cross diffusion, Spw also contains a stability
        // split. Its slope is not claimed to be the full nonlinear Jacobian.
        EXPECT_NEAR(q.spw-std::min(cross_source,0.0)/omega,fd,1e-8);
        c.sst_omega_newton_sink=false;
        const auto legacy=sst_cell_sources(k,omega,S,F1,F2,cross,c);
        EXPECT_NEAR(legacy.sw+legacy.spw*omega,omega_physical,1e-11);
    });
    run_case("v16_wall_diffusion_only_molecular_at_fixed_value_wall", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geom=build_fv_geometry(mesh);
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        c.density=1.2;
        c.molecular_viscosity=1e-5;
        const double molecular=c.density*c.molecular_viscosity;
        std::vector<double> gamma(mesh.n_cells(),20.0*molecular);
        ScalarBoundaryConditions bc;
        for (const char* patch: {"x0","x1","walls"})
            bc[patch]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        auto face_gamma=sst_face_diffusion_coefficients(mesh,geom,gamma,bc,nullptr,c);
        for (const auto value : face_gamma)
            EXPECT_NEAR(value,molecular,1e-16);
        Field<double,Location::FACE> zero_flux(mesh.n_faces(),"flux","kg/s",1);
        zero_flux.fill(0.0);
        Field<double,Location::CELL> zero(mesh.n_cells(),"z","",1);
        zero.fill(0.0);
        Field<double,Location::CELL> k(mesh.n_cells(),"k","m2/s2",1);
        k.fill(0.4);
        const auto legacy=assemble_scalar_equation(
            mesh,geom,zero_flux,0.0,zero,zero,bc,true,
            nullptr,nullptr,nullptr,&gamma,ConvectionScheme::UPWIND,&k);
        const auto physical=assemble_scalar_equation(
            mesh,geom,zero_flux,0.0,zero,zero,bc,true,
            nullptr,nullptr,nullptr,&gamma,ConvectionScheme::UPWIND,&k,
            nullptr,true,&face_gamma);
        EXPECT_NEAR(legacy.diagonal[0]/physical.diagonal[0],20.0,1e-12);
        const auto budget=sst_cell_discrete_budget(
            mesh,geom,zero_flux,k,zero,zero,gamma,bc,nullptr,c,0,c.k_min);
        EXPECT_NEAR(budget.full_residual,budget.recomposed,1e-13);
        EXPECT_NEAR(budget.wall_owner_diffusion_flux,budget.wall_orthogonal_flux,1e-13);
        EXPECT_TRUE(budget.dirichlet_wall_faces==mesh.n_faces());
        c.sst_molecular_wall_diffusion=true;
        const auto budget_modified=sst_cell_discrete_budget(
            mesh,geom,zero_flux,k,zero,zero,gamma,bc,nullptr,c,0,c.k_min);
        EXPECT_NEAR(budget_modified.full_residual,budget_modified.recomposed,1e-13);
        EXPECT_NEAR(budget_modified.wall_molecular_diffusion_flux,
                    budget_modified.wall_orthogonal_flux,1e-13);
        EXPECT_NEAR(budget.full_residual/budget_modified.full_residual,20.0,1e-12);
    });
    run_case("v16_wall_molecular_diffusion_does_not_change_nonwall_dirichlet", [] {
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto geometry = build_fv_geometry(mesh);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.density = 1.2;
        controls.molecular_viscosity = 1e-5;
        const double wall_gamma = controls.density * controls.molecular_viscosity;
        std::vector<double> gamma(mesh.n_cells(), 15.0 * wall_gamma);
        ScalarBoundaryConditions bcs;
        bcs["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        bcs["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.2, 0.0};
        for (const char* name: {"side_low", "side_high", "bottom", "top"})
            bcs[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        const auto result = sst_face_diffusion_coefficients(
            mesh, geometry, gamma, bcs, nullptr, controls);
        EXPECT_TRUE(result.size() == mesh.n_faces());
        EXPECT_NEAR(result[0], wall_gamma, 1e-16);
        EXPECT_NEAR(result[1], gamma[0], 1e-16);
        for (std::size_t face=2; face<mesh.n_faces(); ++face)
            EXPECT_NEAR(result[face], gamma[0], 1e-16);
    });
    run_case("v16_sst_nonlinear_row_budget_has_no_hidden_flux_terms", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"flux","kg/s",1);
        flux.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),omega(1,"omega","1/s",1),
            strain(1,"strain","1/s",1),wall(1,"wall","m",1);
        k(0)=0.4;omega(0)=20.0;strain(0)=4.0;wall(0)=0.2;
        ScalarBoundaryConditions kb,wb;
        for (const char* patch: {"x0","x1","walls"}) {
            kb[patch]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
            wb[patch]={ScalarBoundaryType::FIXED_VALUE,2.0,0.0};
        }
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        for (const bool molecular : {false,true}) {
            c.sst_molecular_wall_diffusion=molecular;
            const auto [bk,bw]=sst_reassembled_discrete_cell_budgets(
                mesh,geom,flux,k,omega,strain,wall,c,kb,wb,nullptr,0);
            const double k_scale=std::max({1.0,std::abs(bk.full_residual),std::abs(bk.recomposed)});
            const double w_scale=std::max({1.0,std::abs(bw.full_residual),std::abs(bw.recomposed)});
            EXPECT_TRUE(std::abs(bk.closure_error) <= 1e-12*k_scale);
            EXPECT_TRUE(std::abs(bw.closure_error) <= 1e-12*w_scale);
            EXPECT_TRUE(bk.diffusion>0.0 && bw.diffusion>0.0);
            EXPECT_TRUE(bk.implicit_source>0.0 && bw.implicit_source>0.0);
            EXPECT_TRUE(bk.explicit_source<0.0 && bw.explicit_source<0.0);
            EXPECT_TRUE(std::abs(bk.nonlinear_source_identity_error)<=1e-12*k_scale);
            EXPECT_TRUE(std::abs(bw.nonlinear_source_identity_error)<=1e-12*w_scale);
        }
    });
    run_case("v15_nonlinear_line_search_merit_is_strict_and_fail_closed", [] {
        const double old=0.2;
        EXPECT_TRUE(sst_nonlinear_step_reduces_merit(old,0.19,0.5));
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(old,old,1.0));
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(old,0.3,0.1));
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(old,-1.0,1.0));
        EXPECT_TRUE(!sst_nonlinear_step_reduces_merit(old,0.01,0.0));
        TurbulenceTransportControls controls;
        controls.model=TurbulenceModel::SST;
        controls.sst_nonlinear_line_search=true;
        controls.sst_picard_damping=0.5;
        EXPECT_THROW(validate_turbulence_controls(controls),std::invalid_argument);
        controls.sst_picard_damping=1.0;
        validate_turbulence_controls(controls);
    });
    run_case("v15_line_search_one_cell_cannot_increase_reassembled_merit", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"flux","kg/s",1);
        flux.fill(0.0);
        Field<double,Location::CELL> k(1,"k","m2/s2",1),w(1,"w","1/s",1),
            strain(1,"strain","1/s",1),wall(1,"y","m",1);
        k(0)=0.4;w(0)=8.0;strain(0)=2.0;wall(0)=0.1;
        ScalarBoundaryConditions kb,ob;
        for(const char* patch:{"x0","x1","walls"}) {
            kb[patch]={ScalarBoundaryType::FIXED_VALUE,0.01,0.0};
            ob[patch]={ScalarBoundaryType::FIXED_VALUE,2.0,0.0};
        }
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        c.sst_nonlinear_line_search=true;
        const auto before=sst_reassembled_returned_state_residual(
            mesh,geom,flux,k,w,strain,wall,c,kb,ob,nullptr);
        const auto result=solve_sst_transport_dynamic_blending(
            mesh,geom,flux,k,w,strain,wall,c,kb,ob,1,1e-12);
        const auto after=sst_reassembled_returned_state_residual(
            mesh,geom,flux,k,w,strain,wall,c,kb,ob,nullptr);
        EXPECT_TRUE(std::isfinite(result.k_residual) &&
                    std::isfinite(result.second_residual));
        EXPECT_TRUE(sst_nonlinear_merit(after) <= sst_nonlinear_merit(before));
        EXPECT_TRUE(!result.converged);
    });
    run_case("v14_sst_records_linear_solution_and_relaxed_state_separately", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"mass_flux","kg/s",1);
        flux.fill(0.0);
        Field<double,Location::CELL> k(mesh.n_cells(),"k","m2/s2",1),
            omega(mesh.n_cells(),"omega","1/s",1),
            strain(mesh.n_cells(),"strain","1/s",1),
            F1(mesh.n_cells(),"F1","1",1),F2(mesh.n_cells(),"F2","1",1);
        k.fill(0.4); omega.fill(20.0); strain.fill(10.0);
        F1.fill(0.8); F2.fill(0.7);
        ScalarBoundaryConditions kb, ob;
        for (const char* patch : {"x0", "x1", "walls"}) {
            kb[patch]={ScalarBoundaryType::FIXED_VALUE,0.03,0.0};
            ob[patch]={ScalarBoundaryType::FIXED_VALUE,2.0,0.0};
        }
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::SST;
        const auto r=solve_sst_transport(mesh,geometry,flux,k,omega,strain,
                                          F1,F2,c,kb,ob,1,1.e-10);
        EXPECT_TRUE(r.last_k_linear_status==SolverStatus::CONVERGED);
        EXPECT_TRUE(r.last_second_linear_status==SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(r.last_k_linear_before));
        EXPECT_TRUE(std::isfinite(r.last_second_linear_before));
        EXPECT_TRUE(std::isfinite(r.last_k_linear_after));
        EXPECT_TRUE(std::isfinite(r.last_second_linear_after));
        EXPECT_TRUE(r.last_k_linear_before<=1.e-10);
        EXPECT_TRUE(r.last_second_linear_before<=1.e-10);
        // A certified algebraic solve can still leave a nonzero residual
        // after relaxation; neither status certifies the SST nonlinear state.
        EXPECT_TRUE(r.last_k_linear_after>=r.last_k_linear_before);
        EXPECT_TRUE(r.last_second_linear_after>=r.last_second_linear_before);
        EXPECT_TRUE(r.k_residual>0 || r.second_residual>0);
    });
    run_case("v13_sst_newton_sink_preserves_the_exact_omega_source", [] {
        TurbulenceTransportControls old;
        old.model = TurbulenceModel::SST;
        old.density = 1.2;
        old.sst_omega_newton_sink = false;
        auto improved=old;
        improved.sst_omega_newton_sink = true;
        const double k=0.35, w=15.0, strain=4.0, f1=0.7, f2=0.4;
        const double cross=-1.3;
        const auto old_s=sst_cell_sources(k,w,strain,f1,f2,cross,old);
        const auto new_s=sst_cell_sources(k,w,strain,f1,f2,cross,improved);
        EXPECT_NEAR(old_s.sk, new_s.sk, 1e-13);
        EXPECT_NEAR(old_s.spk, new_s.spk, 1e-13);
        EXPECT_NEAR(old_s.gamma_k, new_s.gamma_k, 1e-13);
        EXPECT_NEAR(old_s.sw + old_s.spw*w,
                    new_s.sw + new_s.spw*w, 1e-12);
        const double beta= f1*old.beta1+(1.0-f1)*old.beta2;
        EXPECT_NEAR(new_s.sw-old_s.sw,old.density*beta*w*w,1e-11);
        EXPECT_NEAR(new_s.spw-old_s.spw,-old.density*beta*w,1e-12);
        EXPECT_TRUE(new_s.spw < old_s.spw);
    });
    run_case("v13_wall_omega_requires_explicit_treatment_choice", [] {
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto geometry = build_fv_geometry(mesh);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity=1e-5;
        controls.sst_wall_omega_factor=6.0;
        const auto legacy=make_sst_wall_omega_face_values(
            mesh,geometry,std::vector<std::string>{"wall"},controls);
        controls.sst_wall_omega_factor=60.0;
        const auto menter=make_sst_wall_omega_face_values(
            mesh,geometry,std::vector<std::string>{"wall"},controls);
        EXPECT_NEAR(menter.values.at("wall")[0],
                    10.0*legacy.values.at("wall")[0],1e-12);
        controls.sst_wall_omega_factor=9.0;
        EXPECT_THROW(validate_turbulence_controls(controls),std::invalid_argument);
    });
    run_case("v13_omega_first_coupling_uses_updated_omega_in_k_equation", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_transaction_unit_cube();
        const auto geometry = build_fv_geometry(mesh);
        Field<double,Location::FACE> flux(mesh.n_faces(),"phi","kg/s",1);
        flux.fill(0.0);
        Field<double,Location::CELL> k(mesh.n_cells(),"k","m2/s2",1),
            omega(mesh.n_cells(),"omega","1/s",1),
            strain(mesh.n_cells(),"s","1/s",1),
            F1(mesh.n_cells(),"f1","1",1), F2(mesh.n_cells(),"f2","1",1);
        k.fill(0.4); omega.fill(20.0); strain.fill(10.0);
        F1.fill(0.8); F2.fill(0.7);
        ScalarBoundaryConditions kb,ob;
        for (const char* patch : {"x0", "x1", "walls"}) {
            kb[patch]={ScalarBoundaryType::FIXED_VALUE,0.03,0.0};
            ob[patch]={ScalarBoundaryType::FIXED_VALUE,2.0,0.0};
        }
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.sst_omega_newton_sink = true;
        auto parallel_k=k,parallel_omega=omega;
        auto sequential_k=k,sequential_omega=omega;
        const auto parallel_result=solve_sst_transport(
            mesh,geometry,flux,parallel_k,parallel_omega,strain,F1,F2,
            controls,kb,ob,1,1e-10);
        controls.sst_sequential_omega_first=true;
        const auto sequential_result=solve_sst_transport(
            mesh,geometry,flux,sequential_k,sequential_omega,strain,F1,F2,
            controls,kb,ob,1,1e-10);
        EXPECT_TRUE(std::isfinite(parallel_result.k_residual));
        EXPECT_TRUE(std::isfinite(sequential_result.k_residual));
        EXPECT_TRUE(sequential_k(0)>0.0 && sequential_omega(0)>0.0);
        EXPECT_NEAR(parallel_omega(0),sequential_omega(0),1e-11);
        EXPECT_TRUE(std::abs(parallel_k(0)-sequential_k(0))>1e-7);
    });
    run_case("sst_dynamic_blending_diagnostics_accumulate_all_bound_events", [] {
        TurbulenceTransportResult aggregate;
        aggregate.k_floor_events = 2;
        aggregate.second_floor_events = 3;
        aggregate.nonfinite_bound_events = 1;
        aggregate.k_total_floor_correction = 0.25;
        aggregate.second_total_floor_correction = 0.5;
        aggregate.k_max_floor_correction = 0.2;
        aggregate.second_max_floor_correction = 0.4;

        TurbulenceTransportResult step;
        step.k_floor_events = 4;
        step.second_floor_events = 5;
        step.nonfinite_bound_events = 2;
        step.k_total_floor_correction = 0.75;
        step.second_total_floor_correction = 1.5;
        step.k_max_floor_correction = 0.6;
        step.second_max_floor_correction = 1.2;

        accumulate_turbulence_diagnostics(aggregate, step);
        EXPECT_TRUE(aggregate.k_floor_events == 6);
        EXPECT_TRUE(aggregate.second_floor_events == 8);
        EXPECT_TRUE(aggregate.nonfinite_bound_events == 3);
        EXPECT_NEAR(aggregate.k_total_floor_correction, 1.0, 0.0);
        EXPECT_NEAR(aggregate.second_total_floor_correction, 2.0, 0.0);
        EXPECT_NEAR(aggregate.k_max_floor_correction, 0.6, 0.0);
        EXPECT_NEAR(aggregate.second_max_floor_correction, 1.2, 0.0);
    });

    run_case("sst_blending_is_bounded", [] {
        const auto [f1, f2] = compute_sst_blending(0.5, 10.0, 0.01, 1.5e-5, 0.09);
        EXPECT_TRUE(f1 >= 0.0 && f1 <= 1.0);
        EXPECT_TRUE(f2 >= 0.0 && f2 <= 1.0);
    });

    run_case("sst_blending_responds_to_wall_distance", [] {
        const auto near_wall = compute_sst_blending(0.5, 10.0, 1.0e-3, 1.5e-5, 0.09);
        const auto far_wall = compute_sst_blending(0.5, 10.0, 10.0, 1.5e-5, 0.09);
        EXPECT_TRUE(near_wall.first > far_wall.first);
        EXPECT_TRUE(near_wall.second > far_wall.second);
    });

    run_case("sst_production_viscosity_uses_cellwise_sst_closure", [] {
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(2, "k", "m2/s2", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> omega(2, "omega", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> strain(2, "strain", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> f2(2, "F2", "1", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> nu_t(2, "nu_t", "m2/s", 1);
        k(0) = 0.5; k(1) = 1.0;
        omega(0) = 10.0; omega(1) = 20.0;
        strain(0) = 2.0; strain(1) = 8.0;
        f2(0) = 0.5; f2(1) = 1.0;
        nu_t.fill(7.0);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.a1 = 0.31;

        refresh_sst_production_eddy_viscosity(k, omega, strain, f2, controls, nu_t);
        EXPECT_NEAR(nu_t(0), 0.31 * 0.5 / std::max(0.31 * 10.0, 2.0 * 0.5), 1e-14);
        EXPECT_NEAR(nu_t(1), 0.31 * 1.0 / std::max(0.31 * 20.0, 8.0), 1e-14);
        EXPECT_TRUE(nu_t(0) != nu_t(1));
    });

    run_case("sst_production_viscosity_rejects_invalid_state_atomically", [] {
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(2), omega(2), strain(2), f2(2), nu_t(2);
        k(0) = 0.5; k(1) = 1.0;
        omega(0) = 10.0; omega(1) = 20.0;
        strain(0) = 2.0; strain(1) = 4.0;
        f2(0) = 0.5; f2(1) = std::numeric_limits<double>::quiet_NaN();
        nu_t(0) = 3.0; nu_t(1) = 4.0;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(refresh_sst_production_eddy_viscosity(k, omega, strain, f2, controls, nu_t),
                     std::invalid_argument);
        EXPECT_NEAR(nu_t(0), 3.0, 0.0);
        EXPECT_NEAR(nu_t(1), 4.0, 0.0);
    });

    run_case("sst_transaction_preflight_rejects_field_mismatch_without_mutation", [] {
        cfdx::core::Mesh mesh;
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(1), omega(1), strain(1), wall(1), nu_t(1);
        cfdx::core::Field<double, cfdx::core::Location::FACE> flux(0);
        k(0) = 0.4; omega(0) = 8.0; strain(0) = 1.0; wall(0) = 0.01; nu_t(0) = 0.02;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, FvGeometry{}, flux, k, omega, strain, wall, controls, nu_t),
            std::invalid_argument);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });


    run_case("sst_transaction_rolls_back_fields_after_nonconverged_transport", [] {
        using cfdx::core::FaceOwnership;
        using cfdx::core::Mesh;
        using cfdx::core::Patch;
        using cfdx::physics::build_fv_geometry;

        Mesh mesh;
        mesh.points().resize(8);
        const double points[8][3] = {
            {0,0,0},{1,0,0},{1,1,0},{0,1,0},
            {0,0,1},{1,0,1},{1,1,1},{0,1,1}
        };
        for (std::size_t i = 0; i < 8; ++i)
            mesh.points().set(i, points[i][0], points[i][1], points[i][2]);
        mesh.faces().push_face({0,3,2,1});
        mesh.faces().push_face({4,5,6,7});
        mesh.faces().push_face({0,1,5,4});
        mesh.faces().push_face({3,7,6,2});
        mesh.faces().push_face({0,4,7,3});
        mesh.faces().push_face({1,2,6,5});
        mesh.ownership().resize(6);
        for (std::size_t f = 0; f < 6; ++f) {
            mesh.ownership().set_owner(f, 0);
            mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
        mesh.cells().push_cell({0,1,2,3,4,5});
        Patch wall;
        wall.name = "wall";
        wall.type = cfdx::core::PatchType::WALL;
        wall.face_ids = {0,1,2,3,4,5};
        mesh.boundary().add_patch(wall);
        const auto geometry = build_fv_geometry(mesh);

        cfdx::core::Field<double, cfdx::core::Location::FACE> flux(
            mesh.n_faces(), "phi", "kg/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(
            mesh.n_cells(), "k", "m2/s2", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> omega(
            mesh.n_cells(), "omega", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> strain(
            mesh.n_cells(), "strain", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> wall_distance(
            mesh.n_cells(), "wall_distance", "m", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> nu_t(
            mesh.n_cells(), "nu_t", "m2/s", 1);
        flux.fill(0.0);
        k.fill(0.4);
        omega.fill(8.0);
        strain.fill(1.0);
        wall_distance.fill(0.1);
        nu_t.fill(0.02);

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        ScalarBoundaryConditions k_bcs;
        k_bcs["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.01, 0.0};
        ScalarBoundaryConditions omega_bcs;
        omega_bcs["wall"] = {ScalarBoundaryType::FIXED_VALUE, 2.0, 0.0};

        auto direct_k = k;
        auto direct_omega = omega;
        const auto direct = solve_sst_transport_dynamic_blending(
            mesh, geometry, flux, direct_k, direct_omega, strain, wall_distance,
            controls, k_bcs, omega_bcs, 1, 1e-12);
        EXPECT_TRUE(!direct.converged);
        EXPECT_TRUE(direct_k(0) != k(0) || direct_omega(0) != omega(0));

        const auto step = advance_sst_transport_transactionally(
            mesh, geometry, flux, k, omega, strain, wall_distance, controls, nu_t,
            k_bcs, omega_bcs, 1, 1e-12);
        EXPECT_TRUE(!step.committed);
        EXPECT_TRUE(!step.transport.converged);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);

        // Malformed dimensions must now be rejected before direct transport
        // mutates k/omega, as well as by the transactional wrapper. Exception
        // rollback after a real inner commit is exercised separately below.
        cfdx::core::Field<double, cfdx::core::Location::CELL> malformed_strain(
            mesh.n_cells(), "strain_bad_dimension", "1/s", 2);
        malformed_strain.fill(1.0);
        auto exception_direct_k = k;
        auto exception_direct_omega = omega;
        EXPECT_THROW(solve_sst_transport_dynamic_blending(
            mesh, geometry, flux, exception_direct_k, exception_direct_omega,
            malformed_strain, wall_distance, controls, k_bcs, omega_bcs, 100, 1e-4),
            std::invalid_argument);
        EXPECT_NEAR(exception_direct_k(0), k(0), 0.0);
        EXPECT_NEAR(exception_direct_omega(0), omega(0), 0.0);

        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, geometry, flux, k, omega, malformed_strain, wall_distance,
            controls, nu_t, k_bcs, omega_bcs, 100, 1e-4), std::invalid_argument);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });

    run_case("sst_strain_rate_magnitude_matches_simple_shear", [] {
        // Ux = gamma*y, Uy = Uz = 0 gives |S| = |gamma|.
        const double gamma = 3.25;
        const std::array<double, 9> grad{
            0.0, gamma, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0};
        EXPECT_NEAR(sst_strain_rate_magnitude(grad), gamma, 1e-14);
    });

    run_case("sst_strain_rate_magnitude_matches_extension", [] {
        // Ux = a*x and Uy = -a*y gives |S| = 2*|a|.
        const double a = 1.75;
        const std::array<double, 9> grad{
            a, 0.0, 0.0,
            0.0, -a, 0.0,
            0.0, 0.0, 0.0};
        EXPECT_NEAR(sst_strain_rate_magnitude(grad), 2.0 * a, 1e-14);
    });

    run_case("sst_strain_rate_magnitude_rejects_nonfinite_gradient", [] {
        const std::array<double, 9> grad{
            0.0, 0.0, 0.0,
            0.0, std::numeric_limits<double>::infinity(), 0.0,
            0.0, 0.0, 0.0};
        EXPECT_THROW(sst_strain_rate_magnitude(grad), std::invalid_argument);
    });

    run_case("sst_strain_rate_uses_explicit_fixed_velocity_boundaries", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);

        const auto zero_gradient = compute_sst_strain_rate(velocity, mesh);
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"inlet", {1.0, 0.0, 0.0}}
        };
        const auto boundary_aware =
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values);

        EXPECT_NEAR(zero_gradient(0), 0.0, 1e-14);
        EXPECT_TRUE(std::isfinite(boundary_aware(0)));
        EXPECT_TRUE(boundary_aware(0) > 0.0);
    });

    run_case("sst_strain_rate_reconstructs_affine_field_with_fixed_boundaries", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);
        velocity.component_data(0)[0] = 2.0;

        // Ux=x is exactly representable on this orthogonal box. Prescribing
        // the affine field at every boundary face should recover dUx/dx=1.
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"wall", {0.0, 0.0, 0.0}},
            {"inlet", {4.0, 0.0, 0.0}},
            {"side_low", {2.0, 0.0, 0.0}},
            {"side_high", {2.0, 0.0, 0.0}},
            {"bottom", {2.0, 0.0, 0.0}},
            {"top", {2.0, 0.0, 0.0}}
        };

        const auto strain =
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values);
        EXPECT_NEAR(strain(0), std::sqrt(2.0), 1e-13);
    });

    run_case("sst_strain_rate_boundary_values_reject_nonfinite_input", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"inlet", {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}}
        };
        EXPECT_THROW(
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values),
            std::invalid_argument);
    });

    run_case("sst_strain_rate_rejects_unknown_fixed_value_patch", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"missing_patch", {1.0, 0.0, 0.0}}
        };
        EXPECT_THROW(
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values),
            std::invalid_argument);
    });

    run_case("sst_blending_rejects_invalid_distance", [] {
        EXPECT_THROW(compute_sst_blending(0.5, 10.0, 0.0, 1.5e-5, 0.09), std::invalid_argument);
    });


    run_case("sst_inlet_turbulence_state_uses_velocity_intensity_and_length_scale", [] {
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.velocity_vector = {3.0, 4.0, 0.0};
        inlet.turbulence_intensity = 0.05;
        inlet.turbulence_length_scale = 0.2;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        const auto state = derive_sst_inlet_turbulence_state(inlet, controls);
        const double expected_k = 1.5 * std::pow(5.0 * 0.05, 2);
        const double expected_omega = std::sqrt(expected_k) /
            (std::pow(controls.beta_star, 0.25) * 0.2);
        EXPECT_NEAR(state.k, expected_k, 1e-14);
        EXPECT_NEAR(state.omega, expected_omega, 1e-14);
    });

    run_case("sst_inlet_turbulence_state_rejects_percentage_intensity", [] {
        cfdx::io::BoundarySpec inlet;
        inlet.type = cfdx::io::BCType::INLET;
        inlet.velocity_magnitude = 10.0;
        inlet.turbulence_intensity = 5.0;
        inlet.turbulence_length_scale = 0.1;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(derive_sst_inlet_turbulence_state(inlet, controls),
                     std::invalid_argument);
    });

    run_case("sst_case_prerequisites_require_explicit_initial_k_omega", [] {
        cfdx::io::CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.01;
        setup.initial_condition.scalar_fields["omega"] = 2.0;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.velocity_magnitude = 10.0;
        inlet.turbulence_intensity = 0.05;
        inlet.turbulence_length_scale = 0.1;
        setup.boundary_conditions.push_back(inlet);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        validate_sst_case_prerequisites(setup, controls);

        setup.initial_condition.scalar_fields.erase("omega");
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
    });

    run_case("sst_case_prerequisites_reject_nonpositive_initial_fields_and_duplicate_patches", [] {
        cfdx::io::CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.01;
        setup.initial_condition.scalar_fields["omega"] = 2.0;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        setup.boundary_conditions.push_back(wall);
        setup.boundary_conditions.push_back(wall);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
        setup.boundary_conditions.clear();
        setup.initial_condition.scalar_fields["k"] = 0.0;
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
        setup.initial_condition.scalar_fields["k"] = 0.01;
        cfdx::io::BoundarySpec unknown;
        unknown.patch_name = "mystery";
        unknown.type = cfdx::io::BCType::UNKNOWN;
        setup.boundary_conditions.push_back(unknown);
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
    });

    run_case("sst_scalar_boundary_contract_accepts_explicit_wall_inlet_outlet", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.01},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 2.0}
        };
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        cfdx::io::BoundarySpec outlet;
        outlet.patch_name = "outlet";
        outlet.type = cfdx::io::BCType::PRESSURE_OUTLET;
        outlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0}
        };
        setup.boundary_conditions = {inlet, wall, outlet};
        validate_sst_scalar_boundary_contract(setup);

        setup.boundary_conditions.push_back(outlet);
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
        setup.boundary_conditions.clear();
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
    });

    run_case("sst_scalar_boundary_contract_rejects_zero_fixed_k", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 2.0}
        };
        setup.boundary_conditions = {inlet};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
    });

    run_case("sst_scalar_boundary_contract_rejects_missing_duplicate_and_incompatible_conditions", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0}
        };
        setup.boundary_conditions = {wall};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);

        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions = {wall};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);

        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.01},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 2.0}
        };
        setup.boundary_conditions = {wall};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
    });


    run_case("sst_production_state_validates_owned_fields", [] {
        SSTProductionState state(2);
        state.k(0) = 0.1; state.k(1) = 0.2;
        state.omega(0) = 10.0; state.omega(1) = 20.0;
        state.nu_t(0) = 0.001; state.nu_t(1) = 0.002;
        state.wall_distance(0) = 0.01; state.wall_distance(1) = 0.02;
        state.strain_rate(0) = 1.0; state.strain_rate(1) = 2.0;

        validate_sst_production_state(state, 2);
        EXPECT_THROW(validate_sst_production_state(state, 3), std::invalid_argument);

        state.wall_distance(1) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(validate_sst_production_state(state, 2), std::invalid_argument);
    });

    run_case("sst_production_state_snapshot_restores_all_owned_fields", [] {
        SSTProductionState state(2);
        state.k(0) = 0.1; state.k(1) = 0.2;
        state.omega(0) = 10.0; state.omega(1) = 20.0;
        state.nu_t(0) = 0.001; state.nu_t(1) = 0.002;
        state.wall_distance(0) = 0.01; state.wall_distance(1) = 0.02;
        state.strain_rate(0) = 1.0; state.strain_rate(1) = 2.0;

        const auto snapshot = snapshot_sst_production_state(state);
        state.k(0) = 9.0;
        state.omega(1) = 99.0;
        state.nu_t(0) = 0.9;
        state.wall_distance(1) = 9.0;
        state.strain_rate(0) = 99.0;

        restore_sst_production_state(state, snapshot);
        EXPECT_NEAR(state.k(0), 0.1, 0.0);
        EXPECT_NEAR(state.k(1), 0.2, 0.0);
        EXPECT_NEAR(state.omega(0), 10.0, 0.0);
        EXPECT_NEAR(state.omega(1), 20.0, 0.0);
        EXPECT_NEAR(state.nu_t(0), 0.001, 0.0);
        EXPECT_NEAR(state.nu_t(1), 0.002, 0.0);
        EXPECT_NEAR(state.wall_distance(0), 0.01, 0.0);
        EXPECT_NEAR(state.wall_distance(1), 0.02, 0.0);
        EXPECT_NEAR(state.strain_rate(0), 1.0, 0.0);
        EXPECT_NEAR(state.strain_rate(1), 2.0, 0.0);
    });

    run_case("sst_coupled_transaction_rolls_back_carrier_flux_and_all_turbulence_fields", [] {
        using cfdx::core::Field;
        Field<double, cfdx::core::Location::CELL> U(2, "U", "m/s", 3);
        Field<double, cfdx::core::Location::CELL> p(2, "p", "Pa", 1);
        Field<double, cfdx::core::Location::FACE> phi(3, "phi", "kg/s", 1);
        SSTProductionState state(2);
        for (std::size_t c = 0; c < 2; ++c) {
            U.component_data(0)[c] = 1.0 + static_cast<double>(c);
            U.component_data(1)[c] = 2.0 + static_cast<double>(c);
            U.component_data(2)[c] = 3.0 + static_cast<double>(c);
            p(c) = 10.0 + static_cast<double>(c);
            state.k(c) = 0.1 + 0.1 * static_cast<double>(c);
            state.omega(c) = 10.0 + 10.0 * static_cast<double>(c);
            state.nu_t(c) = 0.001 + 0.001 * static_cast<double>(c);
            state.wall_distance(c) = 0.01 + 0.01 * static_cast<double>(c);
            state.strain_rate(c) = 1.0 + static_cast<double>(c);
        }
        phi(0) = 0.25; phi(1) = -0.5; phi(2) = 0.75;

        SSTCoupledStateTransaction transaction(U, p, phi, state);
        transaction.begin();
        U.component_data(0)[0] = -100.0;
        p(1) = -200.0;
        phi(2) = 300.0;
        state.k(0) = 400.0;
        state.omega(1) = 500.0;
        state.nu_t(0) = 600.0;
        state.wall_distance(1) = 700.0;
        state.strain_rate(0) = 800.0;
        transaction.rollback();

        EXPECT_TRUE(!transaction.active());
        EXPECT_NEAR(U.component_data(0)[0], 1.0, 0.0);
        EXPECT_NEAR(p(1), 11.0, 0.0);
        EXPECT_NEAR(phi(2), 0.75, 0.0);
        EXPECT_NEAR(state.k(0), 0.1, 0.0);
        EXPECT_NEAR(state.omega(1), 20.0, 0.0);
        EXPECT_NEAR(state.nu_t(0), 0.001, 0.0);
        EXPECT_NEAR(state.wall_distance(1), 0.02, 0.0);
        EXPECT_NEAR(state.strain_rate(0), 1.0, 0.0);
    });

    run_case("sst_coupled_transaction_commit_keeps_all_candidate_fields", [] {
        using cfdx::core::Field;
        Field<double, cfdx::core::Location::CELL> U(1, "U", "m/s", 3);
        Field<double, cfdx::core::Location::CELL> p(1, "p", "Pa", 1);
        Field<double, cfdx::core::Location::FACE> phi(1, "phi", "kg/s", 1);
        SSTProductionState state(1);
        U.fill(1.0); p.fill(2.0); phi.fill(3.0);
        state.k(0) = 0.1; state.omega(0) = 10.0;
        state.nu_t(0) = 0.01; state.wall_distance(0) = 0.1;
        state.strain_rate(0) = 1.0;

        SSTCoupledStateTransaction transaction(U, p, phi, state);
        transaction.begin();
        U.component_data(0)[0] = 4.0;
        p(0) = 5.0;
        phi(0) = 6.0;
        state.k(0) = 0.2;
        transaction.commit();

        EXPECT_TRUE(!transaction.active());
        EXPECT_NEAR(U.component_data(0)[0], 4.0, 0.0);
        EXPECT_NEAR(p(0), 5.0, 0.0);
        EXPECT_NEAR(phi(0), 6.0, 0.0);
        EXPECT_NEAR(state.k(0), 0.2, 0.0);
        EXPECT_THROW(transaction.rollback(), std::logic_error);
    });


    run_case("sst_runtime_boundary_adapter_preserves_wall_omega_obligation", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions.push_back(wall);

        const auto runtime = make_sst_runtime_boundary_contract(setup);
        EXPECT_TRUE(runtime.k.count("wall") == 1);
        EXPECT_TRUE(runtime.k.at("wall").type == ScalarBoundaryType::FIXED_VALUE);
        EXPECT_NEAR(runtime.k.at("wall").value, 0.0, 0.0);
        EXPECT_TRUE(runtime.omega.count("wall") == 0);
        EXPECT_TRUE(runtime.wall_omega_sst_patches.size() == 1);
        EXPECT_TRUE(runtime.wall_omega_sst_patches.front() == "wall");
        EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime),
                     std::invalid_argument);
    });

    run_case("sst_runtime_boundary_adapter_maps_fixed_and_zero_gradient", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.4},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 8.0}
        };
        cfdx::io::BoundarySpec outlet;
        outlet.patch_name = "outlet";
        outlet.type = cfdx::io::BCType::OUTLET;
        outlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0}
        };
        setup.boundary_conditions = {inlet, outlet};

        const auto runtime = make_sst_runtime_boundary_contract(setup);
        EXPECT_NEAR(runtime.k.at("inlet").value, 0.4, 0.0);
        EXPECT_NEAR(runtime.omega.at("inlet").value, 8.0, 0.0);
        EXPECT_TRUE(runtime.k.at("outlet").type == ScalarBoundaryType::ZERO_GRADIENT);
        EXPECT_TRUE(runtime.omega.at("outlet").type == ScalarBoundaryType::ZERO_GRADIENT);
        EXPECT_TRUE(runtime.wall_omega_sst_patches.empty());
        bool supported = true;
        try {
            require_supported_sst_runtime_boundaries(runtime);
        } catch (...) {
            supported = false;
        }
        EXPECT_TRUE(supported);
    });

    run_case("sst_wall_omega_face_operator_only_populates_requested_wall_faces", [] {
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto geometry = build_fv_geometry(mesh);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        controls.beta1 = 0.075;
        const auto face_values = make_sst_wall_omega_face_values(
            mesh, geometry, std::vector<std::string>{"wall"}, controls);
        EXPECT_TRUE(face_values.values.count("wall") == 1);
        EXPECT_TRUE(face_values.values.count("inlet") == 0);
        EXPECT_NEAR(face_values.values.at("wall")[0], 2.0e-4, 1.0e-15);
        for (std::size_t face = 1; face < mesh.n_faces(); ++face)
            EXPECT_TRUE(!std::isfinite(face_values.values.at("wall")[face]));
        const auto runtime = SSTRuntimeBoundaryContract{{}, {}, {"wall"}};
        EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime),
                     std::invalid_argument);
        EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime, &face_values),
                     std::invalid_argument);
        require_supported_sst_runtime_boundaries(runtime, mesh, geometry, &face_values);
        // Inverting the owner-outward area vector makes the geometry
        // inconsistent. An absolute normal projection previously hid it.
        auto inward_geometry=geometry;
        inward_geometry.face_area_vectors[0] = inward_geometry.face_area_vectors[0] * -1.0;
        EXPECT_THROW(make_sst_wall_omega_face_values(
            mesh,inward_geometry,std::vector<std::string>{"wall"},controls),
            std::invalid_argument);
        EXPECT_THROW(require_supported_sst_runtime_boundaries(
            runtime,mesh,inward_geometry,&face_values),
            std::invalid_argument);
    });

    run_case("sst_wall_omega_value_uses_resolved_wall_formula", [] {
        const double value = sst_wall_omega_value(1.0e-5, 0.075, 2.0e-3);
        EXPECT_NEAR(value, 200.0, 1.0e-10);
        EXPECT_THROW(sst_wall_omega_value(1.0e-5, 0.075, 0.0),
                     std::invalid_argument);
        EXPECT_THROW(sst_wall_omega_value(0.0, 0.075, 2.0e-3),
                     std::invalid_argument);
    });

    run_case("sst_production_state_initializes_from_explicit_case_and_mesh_fields", [] {
        auto mesh = make_sst_wall_operator_test_mesh();
        cfdx::io::CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.12;
        setup.initial_condition.scalar_fields["omega"] = 18.0;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions = {wall};

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        const auto state = initialize_sst_production_state(
            mesh, setup, controls, std::vector<double>{0.01},
            std::vector<double>{2.0});
        validate_sst_production_state(state, mesh.n_cells());
        EXPECT_NEAR(state.k(0), 0.12, 0.0);
        EXPECT_NEAR(state.omega(0), 18.0, 0.0);
        EXPECT_NEAR(state.wall_distance(0), 0.01, 0.0);
        EXPECT_NEAR(state.strain_rate(0), 2.0, 0.0);
        EXPECT_TRUE(std::isfinite(state.nu_t(0)) && state.nu_t(0) >= 0.0);

        EXPECT_THROW(initialize_sst_production_state(
            mesh, setup, controls, std::vector<double>{0.0},
            std::vector<double>{2.0}), std::invalid_argument);
        EXPECT_THROW(initialize_sst_production_state(
            mesh, setup, controls, std::vector<double>{0.01},
            std::vector<double>{std::numeric_limits<double>::quiet_NaN()}),
            std::invalid_argument);
    });

    run_case("sst_production_state_builds_wall_distance_and_strain_from_mesh", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);

        cfdx::io::CaseSetup setup;
        setup.physics_model = "incompressible";
        setup.turbulence_model = "SST";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.12;
        setup.initial_condition.scalar_fields["omega"] = 18.0;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions = {wall};

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"inlet", {1.0, 0.0, 0.0}}
        };

        auto state = initialize_sst_production_state_from_mesh(
            mesh, setup, controls, velocity, fixed_values);
        validate_sst_production_state(state, mesh.n_cells());
        EXPECT_NEAR(state.k(0), 0.12, 0.0);
        EXPECT_NEAR(state.omega(0), 18.0, 0.0);
        EXPECT_TRUE(std::isfinite(state.wall_distance(0)) &&
                    state.wall_distance(0) > 0.0);
        EXPECT_TRUE(std::isfinite(state.strain_rate(0)) &&
                    state.strain_rate(0) > 0.0);
        EXPECT_TRUE(std::isfinite(state.nu_t(0)) && state.nu_t(0) >= 0.0);

        const double initial_strain = state.strain_rate(0);
        const std::map<std::string, std::array<double, 3>> updated_values{
            {"inlet", {2.0, 0.0, 0.0}}
        };
        refresh_sst_production_derived_state_from_mesh(
            state, mesh, controls, velocity, updated_values);
        EXPECT_TRUE(state.strain_rate(0) > initial_strain);
        validate_sst_production_state(state, mesh.n_cells());

        const double valid_strain = state.strain_rate(0);
        const std::map<std::string, std::array<double, 3>> invalid_values{
            {"missing_patch", {1.0, 0.0, 0.0}}
        };
        EXPECT_THROW(refresh_sst_production_derived_state_from_mesh(
            state, mesh, controls, velocity, invalid_values), std::invalid_argument);
        EXPECT_NEAR(state.strain_rate(0), valid_strain, 0.0);

        Field<double, Location::CELL> invalid_velocity(mesh.n_cells(), "bad", "m/s", 1);
        EXPECT_THROW(initialize_sst_production_state_from_mesh(
            mesh, setup, controls, invalid_velocity, fixed_values),
            std::invalid_argument);
    });

    run_case("sst_dat_restart_restores_named_scalars_and_rebuilds_derived_state", [] {
        using cfdx::core::Field;
        using cfdx::core::Location;
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto path = std::filesystem::temp_directory_path() /
            "cfdx_sst_restart_hydration.dat";
        std::filesystem::remove(path);
        std::filesystem::remove(path.string() + ".tmp");

        Field<double, Location::CELL> saved_u(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> saved_p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::CELL> saved_k(mesh.n_cells(), "k", "m2/s2", 1);
        Field<double, Location::CELL> saved_omega(mesh.n_cells(), "omega", "1/s", 1);
        saved_u.component_data(0)[0] = 0.0;
        saved_u.component_data(1)[0] = 0.0;
        saved_u.component_data(2)[0] = 0.0;
        saved_p(0) = 14.0;
        saved_k(0) = 0.42;
        saved_omega(0) = 125.0;
        cfdx::io::DatRestartFields write_fields{};
        write_fields.k = &saved_k;
        write_fields.omega = &saved_omega;
        cfdx::io::write_dat_restart_fields(
            path.string(), mesh, saved_u, saved_p, write_fields, 17, 2.5);

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.component_data(0)[0] = 3.0;
        U.component_data(1)[0] = -2.0;
        U.component_data(2)[0] = 1.0;
        p(0) = -9.0;
        SSTProductionState state(mesh.n_cells());
        state.k(0) = 9.0;
        state.omega(0) = 9.0;
        state.nu_t(0) = 9.0;
        state.wall_distance(0) = 0.5;

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.5e-5;
        const std::map<std::string, std::array<double, 3>> fixed_velocity_boundary_values{
            {"inlet", {4.0, 0.0, 0.0}}
        };
        const auto metadata = restore_sst_production_state_from_dat(
            path.string(), mesh, U, p, state, std::vector<double>{0.01}, controls,
            fixed_velocity_boundary_values);

        EXPECT_TRUE(metadata.iteration == 17);
        EXPECT_NEAR(metadata.time, 2.5, 0.0);
        EXPECT_NEAR(U.component_data(0)[0], 0.0, 0.0);
        EXPECT_NEAR(p(0), 14.0, 0.0);
        EXPECT_NEAR(state.k(0), 0.42, 0.0);
        EXPECT_NEAR(state.omega(0), 125.0, 0.0);
        EXPECT_NEAR(state.wall_distance(0), 0.01, 0.0);
        EXPECT_TRUE(std::isfinite(state.strain_rate(0)) && state.strain_rate(0) > 0.0);
        EXPECT_TRUE(std::isfinite(state.nu_t(0)) && state.nu_t(0) >= 0.0);
        validate_sst_production_state(state, mesh.n_cells());
        std::filesystem::remove(path);
    });

    run_case("sst_dat_restart_missing_omega_leaves_carrier_and_state_unchanged", [] {
        using cfdx::core::Field;
        using cfdx::core::Location;
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto path = std::filesystem::temp_directory_path() /
            "cfdx_sst_restart_missing_omega.dat";
        std::filesystem::remove(path);
        std::filesystem::remove(path.string() + ".tmp");

        Field<double, Location::CELL> saved_u(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> saved_p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::CELL> saved_k(mesh.n_cells(), "k", "m2/s2", 1);
        saved_u.component_data(0)[0] = 1.0;
        saved_u.component_data(1)[0] = 2.0;
        saved_u.component_data(2)[0] = 3.0;
        saved_p(0) = 4.0;
        saved_k(0) = 0.2;
        cfdx::io::DatRestartFields write_fields{};
        write_fields.k = &saved_k;
        cfdx::io::write_dat_restart_fields(
            path.string(), mesh, saved_u, saved_p, write_fields, 4, 1.0);

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.component_data(0)[0] = 8.0;
        U.component_data(1)[0] = 7.0;
        U.component_data(2)[0] = 6.0;
        p(0) = 5.0;
        SSTProductionState state(mesh.n_cells());
        state.k(0) = 0.3;
        state.omega(0) = 6.0;
        state.nu_t(0) = 0.04;
        state.wall_distance(0) = 0.01;
        state.strain_rate(0) = 2.0;
        const auto old_state = state;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;

        EXPECT_THROW(restore_sst_production_state_from_dat(
            path.string(), mesh, U, p, state, std::vector<double>{0.01}, controls),
            std::runtime_error);
        EXPECT_NEAR(U.component_data(0)[0], 8.0, 0.0);
        EXPECT_NEAR(U.component_data(1)[0], 7.0, 0.0);
        EXPECT_NEAR(U.component_data(2)[0], 6.0, 0.0);
        EXPECT_NEAR(p(0), 5.0, 0.0);
        EXPECT_NEAR(state.k(0), old_state.k(0), 0.0);
        EXPECT_NEAR(state.omega(0), old_state.omega(0), 0.0);
        EXPECT_NEAR(state.nu_t(0), old_state.nu_t(0), 0.0);
        EXPECT_NEAR(state.wall_distance(0), old_state.wall_distance(0), 0.0);
        EXPECT_NEAR(state.strain_rate(0), old_state.strain_rate(0), 0.0);
        std::filesystem::remove(path);
    });

    run_case("sst_transaction_rejects_non_scalar_transport_fields", [] {
        using namespace cfdx::core;
        const Mesh mesh = make_sst_transaction_unit_cube();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
        Field<double, Location::CELL> k(mesh.n_cells(), "k", "m2/s2", 1);
        Field<double, Location::CELL> omega(mesh.n_cells(), "omega", "1/s", 1);
        Field<double, Location::CELL> strain(mesh.n_cells(), "strain", "1/s", 1);
        Field<double, Location::CELL> wall(mesh.n_cells(), "wall_distance", "m", 1);
        Field<double, Location::CELL> nu_t(mesh.n_cells(), "nu_t", "m2/s", 1);
        Field<double, Location::CELL> malformed_k(mesh.n_cells(), "k_bad", "m2/s2", 2);
        flux.fill(0.0);
        k(0) = 0.4;
        omega(0) = 8.0;
        strain(0) = 1.0;
        wall(0) = 0.01;
        nu_t(0) = 0.02;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;

        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, geometry, flux, malformed_k, omega, strain, wall,
            controls, nu_t, {}, {}, 1, 1e-8), std::invalid_argument);

        Field<double, Location::FACE> malformed_flux(
            mesh.n_faces(), "phi_bad", "kg/s", 2);
        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, geometry, malformed_flux, k, omega, strain, wall,
            controls, nu_t, {}, {}, 1, 1e-8), std::invalid_argument);
    });

    run_case("sst_transaction_rolls_back_k_omega_after_nonconverged_transport", [] {
        using namespace cfdx::core;
        const Mesh mesh = make_sst_transaction_unit_cube();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
        Field<double, Location::CELL> k(mesh.n_cells(), "k", "m2/s2", 1);
        Field<double, Location::CELL> omega(mesh.n_cells(), "omega", "1/s", 1);
        Field<double, Location::CELL> strain(mesh.n_cells(), "strain", "1/s", 1);
        Field<double, Location::CELL> wall(mesh.n_cells(), "wall_distance", "m", 1);
        Field<double, Location::CELL> nu_t(mesh.n_cells(), "nu_t", "m2/s", 1);
        flux.fill(0.0);
        k(0) = 0.4;
        omega(0) = 8.0;
        strain(0) = 1.0;
        wall(0) = 0.01;
        nu_t(0) = 0.02;

        ScalarBoundaryConditions k_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 0.9, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 0.2, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };
        ScalarBoundaryConditions omega_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 20.0, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 4.0, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;

        const auto step = advance_sst_transport_transactionally(
            mesh, geometry, flux, k, omega, strain, wall, controls, nu_t,
            k_bcs, omega_bcs, 1, 1e-30);

        EXPECT_TRUE(!step.committed);
        EXPECT_TRUE(!step.transport.converged);
        EXPECT_TRUE(step.transport.iterations == 1);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });


    run_case("sst_production_iterate_transaction_rolls_back_derived_and_carrier_state", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_transaction_unit_cube();
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::FACE> phi(mesh.n_faces(), "phi", "kg/s", 1);
        U.fill(0.0);
        U.component_data(0)[0] = 0.5;
        p.fill(2.0);
        phi.fill(0.0);

        SSTProductionState state(mesh.n_cells());
        state.k(0) = 0.4;
        state.omega(0) = 8.0;
        state.nu_t(0) = 0.02;
        state.wall_distance(0) = 0.25;
        state.strain_rate(0) = 3.0;
        const auto original_state = snapshot_sst_production_state(state);
        const auto original_U = U;
        const auto original_p = p;
        const auto original_phi = phi;

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"x0", {1.0, 0.0, 0.0}}
        };
        const ScalarBoundaryConditions k_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 0.9, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 0.2, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };
        const ScalarBoundaryConditions omega_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 20.0, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 4.0, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };

        const auto result = advance_sst_production_iterate_transactionally(
            mesh, U, p, phi, state, controls, fixed_values, {}, k_bcs,
            omega_bcs, 1, 1e-30);

        EXPECT_TRUE(!result.committed);
        EXPECT_TRUE(!result.transport.converged);
        EXPECT_TRUE(result.reason.find("fully rolled back") != std::string::npos);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            EXPECT_NEAR(U.component_data(0)[c], original_U.component_data(0)[c], 0.0);
            EXPECT_NEAR(U.component_data(1)[c], original_U.component_data(1)[c], 0.0);
            EXPECT_NEAR(U.component_data(2)[c], original_U.component_data(2)[c], 0.0);
            EXPECT_NEAR(p(c), original_p(c), 0.0);
            EXPECT_NEAR(state.k(c), original_state.k(c), 0.0);
            EXPECT_NEAR(state.omega(c), original_state.omega(c), 0.0);
            EXPECT_NEAR(state.nu_t(c), original_state.nu_t(c), 0.0);
            EXPECT_NEAR(state.wall_distance(c), original_state.wall_distance(c), 0.0);
            EXPECT_NEAR(state.strain_rate(c), original_state.strain_rate(c), 0.0);
        }
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi(f), original_phi(f), 0.0);
    });


    run_case("sst_relative_eddy_viscosity_change_is_dimensionless", [] {
        SSTProductionState previous(2);
        SSTProductionState current(2);
        for (std::size_t cell = 0; cell < 2; ++cell) {
            previous.k(cell) = current.k(cell) = 0.1;
            previous.omega(cell) = current.omega(cell) = 1.0;
            previous.wall_distance(cell) = current.wall_distance(cell) = 0.5;
            previous.strain_rate(cell) = current.strain_rate(cell) = 0.0;
        }
        previous.nu_t(0) = 1e-5;
        current.nu_t(0) = 2e-5;
        previous.nu_t(1) = current.nu_t(1) = 1e-3;

        EXPECT_NEAR(
            sst_relative_eddy_viscosity_change(previous, current, 1e-5),
            0.5, 1e-14);
        EXPECT_NEAR(
            sst_relative_eddy_viscosity_change(current, current, 1e-5),
            0.0, 0.0);
        EXPECT_THROW(
            sst_relative_eddy_viscosity_change(previous, current, 0.0),
            std::invalid_argument);
    });

    run_case("sst_production_runtime_requires_explicit_valid_state_and_excludes_prescribed_nu_t", [] {
        IncompressibleSolverControls controls;
        SSTProductionRuntime runtime;
        controls.sst_runtime = &runtime;

        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);

        SSTProductionState state(1);
        state.k(0) = 0.1;
        state.omega(0) = 1.0;
        state.nu_t(0) = 0.01;
        state.wall_distance(0) = 0.5;
        state.strain_rate(0) = 0.2;
        runtime.state = &state;
        runtime.transport.model = TurbulenceModel::SST;
        runtime.transport.density = controls.density;
        runtime.transport.molecular_viscosity = controls.kinematic_viscosity;
        validate_incompressible_controls(controls, 1);

        runtime.transport.density *= 2.0;
        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);
        runtime.transport.density = controls.density;

        controls.turbulent_viscosity = 1.0e-5;
        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);
    });

    run_case("sst_production_runtime_requires_explicit_valid_state_and_excludes_prescribed_nu_t", [] {
        IncompressibleSolverControls controls;
        SSTProductionRuntime runtime;
        controls.sst_runtime = &runtime;

        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);

        SSTProductionState state(1);
        state.k(0) = 0.1;
        state.omega(0) = 1.0;
        state.nu_t(0) = 0.01;
        state.wall_distance(0) = 0.5;
        state.strain_rate(0) = 0.2;
        runtime.state = &state;
        runtime.transport.model = TurbulenceModel::SST;
        runtime.transport.density = controls.density;
        runtime.transport.molecular_viscosity = controls.kinematic_viscosity;
        validate_incompressible_controls(controls, 1);

        runtime.transport.density *= 2.0;
        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);
        runtime.transport.density = controls.density;

        controls.turbulent_viscosity = 1.0e-5;
        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);
    });

    run_case("sst_outer_iteration_transaction_rolls_back_after_inner_sst_commit", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_transaction_unit_cube();
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::FACE> phi(mesh.n_faces(), "phi", "kg/s", 1);
        U.fill(0.0);
        U.component_data(0)[0] = 0.25;
        p.fill(3.0);
        phi.fill(0.0);

        SSTProductionState state(mesh.n_cells());
        state.k.fill(0.4);
        state.omega.fill(8.0);
        state.nu_t.fill(0.02);
        state.wall_distance.fill(0.1);
        state.strain_rate.fill(1.0);

        const auto original_U = U;
        const auto original_p = p;
        const auto original_phi = phi;
        const auto original_sst = snapshot_sst_production_state(state);

        // The outer transaction remains active across a successfully committed
        // inner SST transaction. A later carrier rejection must still restore
        // the entire coupled iterate to the pre-outer-iteration snapshot.
        SSTCoupledStateTransaction outer(U, p, phi, state);
        outer.begin();
        {
            SSTCoupledStateTransaction inner(U, p, phi, state);
            inner.begin();
            U.component_data(0)[0] = 1.25;
            p(0) = 9.0;
            phi(0) = 2.0;
            state.k(0) = 0.8;
            state.omega(0) = 12.0;
            state.nu_t(0) = 0.05;
            state.wall_distance(0) = 0.2;
            state.strain_rate(0) = 4.0;
            inner.commit();
        }

        EXPECT_NEAR(U.component_data(0)[0], 1.25, 0.0);
        EXPECT_NEAR(state.k(0), 0.8, 0.0);
        outer.rollback();

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(U.component_data(d)[c], original_U.component_data(d)[c], 0.0);
            EXPECT_NEAR(p(c), original_p(c), 0.0);
            EXPECT_NEAR(state.k(c), original_sst.k(c), 0.0);
            EXPECT_NEAR(state.omega(c), original_sst.omega(c), 0.0);
            EXPECT_NEAR(state.nu_t(c), original_sst.nu_t(c), 0.0);
            EXPECT_NEAR(state.wall_distance(c), original_sst.wall_distance(c), 0.0);
            EXPECT_NEAR(state.strain_rate(c), original_sst.strain_rate(c), 0.0);
        }
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi(f), original_phi(f), 0.0);
    });

    run_case("v17_pseudo_momentum_preserves_old_state_residual_and_steady_fixed_point", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        phi.fill(0.0);
        const auto masses=make_v17_pseudo_momentum_masses(
            mesh,geometry,phi,std::vector<double>{1.0},2.0,0.5);
        EXPECT_TRUE(masses.size()==1 && masses[0]>0.0);
        Field<double,Location::CELL> old(mesh.n_cells(),"old","m/s",3);
        old.fill(0.0); old.component_data(0)[0]=3.0;
        ScalarEquation eq;
        eq.matrix=SparseMatrix(1,1);
        eq.matrix.push_back(0,0,5.0); eq.matrix.finalize();
        eq.rhs=Vector(1,4.0);
        eq.diagonal={5.0};
        add_v17_pseudo_momentum_term(eq,old,masses,0);
        const double m=masses[0];
        EXPECT_NEAR(eq.diagonal[0],5.0+m,1e-12);
        EXPECT_NEAR(eq.rhs(0),4.0+3.0*m,1e-12);
        EXPECT_NEAR((5.0+m)*3.0-eq.rhs(0),5.0*3.0-4.0,1e-12);
        EXPECT_NEAR((5.0+m)*4.0-eq.rhs(0),5.0*4.0-4.0+m,1e-12);
        EXPECT_THROW(add_v17_pseudo_momentum_term(eq,old,{},0),std::invalid_argument);
    });
    run_case("v17_adaptive_CFL_falls_on_bad_continuity_and_grows_on_joint_progress", [] {
        V17PseudoMomentumControls cfg;
        EXPECT_NEAR(next_v17_pseudo_cfl(0.3,0.2,0.1,0.1,0.02,cfg),0.375,1e-14);
        EXPECT_NEAR(next_v17_pseudo_cfl(0.3,0.2,0.1,0.1,0.13,cfg),0.15,1e-14);
        EXPECT_NEAR(next_v17_pseudo_cfl(0.03,0.2,0.3,0.1,0.2,cfg),0.03,1e-14);
        cfg.minimum_cfl=2.0;
        EXPECT_THROW(validate_v17_pseudo_momentum_controls(cfg),std::invalid_argument);
    });
    run_case("v17_pseudo_is_rejected_for_physical_time_and_monolithic_coupling", [] {
        IncompressibleSolverControls cfg;
        cfg.pseudo_momentum.enabled=true;
        cfg.algorithm=PressureVelocityAlgorithm::COUPLED;
        EXPECT_THROW(validate_incompressible_controls(cfg,1),std::invalid_argument);
        cfg.algorithm=PressureVelocityAlgorithm::SIMPLE;
        cfg.transient.enabled=true;
        EXPECT_THROW(validate_incompressible_controls(cfg,1),std::invalid_argument);
    });
    run_case("v17_authoritative_face_flux_patch_and_cell_closure", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        phi.fill(0.0); phi(4)=-3.0; phi(5)=+3.25;
        const auto balance=audit_v17_authoritative_flux(mesh,geometry,phi);
        EXPECT_NEAR(balance.boundary_signed,0.25,1e-14);
        EXPECT_NEAR(balance.cell_signed,0.25,1e-14);
        EXPECT_NEAR(balance.cell_linf,0.25,1e-14);
        EXPECT_NEAR(balance.topology_error,0.0,1e-14);
        EXPECT_TRUE(balance.patches.size()==3);
        EXPECT_NEAR(balance.patches.at("x0").signed_flux,-3.0,1e-14);
        EXPECT_NEAR(balance.patches.at("x1").signed_flux,+3.25,1e-14);
        phi(0)=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(audit_v17_authoritative_flux(mesh,geometry,phi),std::invalid_argument);
    });

    run_case("v17_projected_Rhie_Chow_HbyA_equals_physical_U_contract_on_stretched_face", [] {
        using namespace cfdx::core;
        Mesh mesh;
        mesh.points().resize(12);
        const double coords[12][3] = {
            {0,0,0},{2,0,0},{3,0,0},{0,1,0},{2,1,0},{3,1,0},
            {0,0,1},{2,0,1},{3,0,1},{0,1,1},{2,1,1},{3,1,1}};
        for(std::size_t i=0;i<12;++i)
            mesh.points().set(i,coords[i][0],coords[i][1],coords[i][2]);
        for(const auto& f : std::vector<std::vector<FaceIndex>>{
            {0,6,9,3},{0,1,7,6},{3,9,10,4},{0,3,4,1},{6,7,10,9},
            {7,10,4,1},{2,5,11,8},{1,2,8,7},{4,10,11,5},
            {1,4,5,2},{7,8,11,10}})
            mesh.faces().push_face(f);
        mesh.ownership().resize(11);
        for(std::size_t f=0;f<11;++f){
            mesh.ownership().set_owner(f,f<=5?0:1);
            mesh.ownership().set_neighbour(f,f==5?1:FaceOwnership::BOUNDARY);
        }
        mesh.cells().push_cell({0,1,2,3,4,5});
        mesh.cells().push_cell({5,6,7,8,9,10});
        Patch a,b;
        a.name="left"; a.type=PatchType::WALL; a.face_ids={0,1,2,3,4};
        b.name="right"; b.type=PatchType::WALL; b.face_ids={6,7,8,9,10};
        mesh.boundary().add_patch(a);
        mesh.boundary().add_patch(b);
        const auto geom=build_fv_geometry(mesh);
        const double wn=pressure_face_neighbour_weight(geom,5,0,1,true);
        EXPECT_TRUE(wn>0.55 && wn<0.75);
        Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
        p(0)=1.0; p(1)=5.0;
        const auto grad=gauss_gradient_with_boundary(p,mesh,geom,{},false,true);
        Field<double,Location::CELL> physical(mesh.n_cells(),"U","m/s",3);
        physical.fill(0.0);
        std::array<std::vector<double>,3> rAU{{{0.2,0.6},{0.2,0.6},{0.2,0.6}}};
        Field<double,Location::CELL> hbya=physical;
        for(std::size_t c=0;c<2;++c)
            for(std::size_t d=0;d<3;++d)
                hbya.component_data(d)[c] =
                    geom.cell_volumes[c]*rAU[d][c]*grad.component_data(d)[c];
        const auto phi_h=make_rhie_chow_mass_flux(mesh,geom,hbya,p,rAU,1.0,
            VelocityBoundaryConditions{},ScalarBoundaryConditions{},
            RhieChowVelocityContract::HBYA_PREDICTOR,false,true);
        const auto phi_u=make_rhie_chow_mass_flux(mesh,geom,physical,p,rAU,1.0,
            VelocityBoundaryConditions{},ScalarBoundaryConditions{},
            RhieChowVelocityContract::PHYSICAL_VELOCITY,false,true);
        EXPECT_NEAR(phi_h(5),phi_u(5),1e-11);
    });

    run_case("v17_pseudo_spectral_excludes_empty_extrusion_and_symmetry", [] {
        using namespace cfdx::core;
        auto mesh=make_sst_transaction_unit_cube();
        const auto original_geom=build_fv_geometry(mesh);
        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
        phi.fill(0.0);
        const auto baseline=make_v17_pseudo_momentum_masses(
            mesh,original_geom,phi,{1.0},1.0,0.3);
        // The last four boundary faces are a 2D extrusion/symmetry;
        // their small normal distance must not create artificial diffusion.
        mesh.boundary().patch(2).type=PatchType::EMPTY;
        const auto geometry=build_fv_geometry(mesh);
        const auto without_empty=make_v17_pseudo_momentum_masses(
            mesh,geometry,phi,{1.0},1.0,0.3);
        EXPECT_TRUE(without_empty[0] < baseline[0]);
        EXPECT_TRUE(without_empty[0] > 0.0);
        // Applying molecular wall viscosity further reduces artificial
        // stiffness without changing the underlying steady equation.
        const auto molecular=make_v17_pseudo_momentum_masses(
            mesh,geometry,phi,{1.0},1.0,0.3,0.01);
        EXPECT_TRUE(molecular[0] < without_empty[0]);
    });

    run_case("v18_sparse_pressure_jacobian_matches_full_rhie_chow_oracle", [] {
        using namespace cfdx::core;
        Mesh mesh;
        mesh.points().resize(12);
        const double coord[12][3]={{0,0,0},{2,0,0},{3,0,0},
            {0,1,0},{2,1,0},{3,1,0},
            {0,0,1},{2,0,1},{3,0,1},
            {0,1,1},{2,1,1},{3,1,1}};
        for(std::size_t i=0;i<12;++i) {
            const double shift=(i==1 || i==4 || i==7 || i==10) ? 0.15 :
                ((i==2 || i==5 || i==8 || i==11) ? 0.55 : 0.0);
            mesh.points().set(i,coord[i][0],coord[i][1]+shift,coord[i][2]);
        }
        for(const auto& faces:std::vector<std::vector<FaceIndex>>{
            {0,6,9,3},{0,1,7,6},{3,9,10,4},{0,3,4,1},{6,7,10,9},
            {7,10,4,1},{2,5,11,8},{1,2,8,7},{4,10,11,5},
            {1,4,5,2},{7,8,11,10}})
            mesh.faces().push_face(faces);
        mesh.ownership().resize(11);
        for(std::size_t f=0;f<11;++f) {
            mesh.ownership().set_owner(f,f<=5?0:1);
            mesh.ownership().set_neighbour(f,f==5?1:FaceOwnership::BOUNDARY);
        }
        mesh.cells().push_cell({0,1,2,3,4,5});
        mesh.cells().push_cell({5,6,7,8,9,10});
        Patch left,right;
        left.name="left";left.type=PatchType::WALL;left.face_ids={0,1,2,3,4};
        right.name="right";right.type=PatchType::WALL;right.face_ids={6,7,8,9,10};
        mesh.boundary().add_patch(left);
        mesh.boundary().add_patch(right);
        const auto geom=build_fv_geometry(mesh);
        ScalarBoundaryConditions pbcs;
        pbcs["left"]={ScalarBoundaryType::FIXED_VALUE,3.0,0.0};
        pbcs["right"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        const std::array<std::vector<double>,3> rAU{{{0.2,0.3},{0.4,0.5},{0.6,0.8}}};
        Field<double,Location::CELL> p(2,"p","Pa",1);
        Field<double,Location::CELL> perturb(2,"perturb","Pa",1);
        Field<double,Location::CELL> H(2,"H","m/s",3);
        H.fill(0.0);
        p(0)=0.7;p(1)=1.9;
        perturb(0)=0.29;perturb(1)=-0.61;
        for(const bool projected: {false,true}) {
            const auto sparse=make_full_rhie_chow_pressure_rows_sparse(
                mesh,geom,rAU,1.2,pbcs,projected);
            const auto reference=make_full_rhie_chow_pressure_rows_reference(
                mesh,geom,rAU,1.2,{},pbcs,projected);
            for(std::size_t i=0;i<2;++i) {
                for(std::size_t j=0;j<2;++j) {
                    const auto value=[](const auto& row,std::size_t index){
                        const auto it=row.find(index);
                        return it==row.end()?0.0:it->second;
                    };
                    EXPECT_NEAR(value(sparse[i],j),value(reference[i],j),1e-12);
                }
            }
            const auto f0=make_rhie_chow_mass_flux(mesh,geom,H,p,rAU,1.2,
                VelocityBoundaryConditions{},pbcs,RhieChowVelocityContract::HBYA_PREDICTOR,
                false,projected);
            for(std::size_t c=0;c<2;++c) p(c)+=perturb(c);
            const auto f1=make_rhie_chow_mass_flux(mesh,geom,H,p,rAU,1.2,
                VelocityBoundaryConditions{},pbcs,RhieChowVelocityContract::HBYA_PREDICTOR,
                false,projected);
            for(std::size_t c=0;c<2;++c) p(c)-=perturb(c);
            std::array<double,2> flux_diff{0.0,0.0};
            for(std::size_t f=0;f<mesh.n_faces();++f) {
                const double d=f1(f)-f0(f);
                flux_diff[mesh.ownership().owner(f)]+=d;
                const auto nr=mesh.ownership().neighbour(f);
                if(nr>=0) flux_diff[static_cast<std::size_t>(nr)]-=d;
            }
            for(std::size_t i=0;i<2;++i) {
                double expected=0.0;
                for(const auto& [j,a]:sparse[i]) expected+=a*perturb(j);
                EXPECT_NEAR(flux_diff[i],expected,1e-11);
            }
        }
    });

    run_case("v19_pressure_flux_step_is_not_underrelaxed_cell_pressure", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> p(1,"p","Pa",1), pflux(1,"pflux","Pa",1);
        p(0)=2.0;
        Vector dp(1,0.0);
        dp(0)=5.0;
        apply_v19_pressure_flux_policy(p,pflux,dp,0.1,true,geometry,
            PressureGaugePolicy::NONE,0,0.0,true);
        EXPECT_NEAR(p(0),2.5,1e-14);
        EXPECT_NEAR(pflux(0),7.0,1e-14);
        // Repeating the legacy policy keeps the pressure and flux states equal.
        p(0)=2.0;
        apply_v19_pressure_flux_policy(p,pflux,dp,0.1,false,geometry,
            PressureGaugePolicy::NONE,0,0.0,true);
        EXPECT_NEAR(p(0),2.5,1e-14);
        EXPECT_NEAR(pflux(0),2.5,1e-14);
        p(0)=2.0;
        dp(0)=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(apply_v19_pressure_flux_policy(p,pflux,dp,0.1,true,
            geometry,PressureGaugePolicy::NONE,0,0.0,true),std::invalid_argument);
        EXPECT_NEAR(p(0),2.0,1e-14);
        // A failed update must not partially overwrite either state.
        p(0)=1.e308;
        pflux(0)=3.0;
        dp(0)=1.e308;
        EXPECT_THROW(apply_v19_pressure_flux_policy(p,pflux,dp,0.5,true,
            geometry,PressureGaugePolicy::NONE,0,0.0,true),std::runtime_error);
        EXPECT_NEAR(p(0),1.e308,1.e292);
        EXPECT_NEAR(pflux(0),3.0,1e-14);
    });
    run_case("v19_full_flux_step_uses_pressure_BC_and_not_relaxed_flux", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> p(1,"p","Pa",1),flux_p(1,"pf","Pa",1);
        Field<double,Location::CELL> H(1,"H","m/s",3);
        H.fill(0.0); p(0)=0.0;
        Vector dp(1,0.0); dp(0)=0.25;
        apply_v19_pressure_flux_policy(p,flux_p,dp,0.1,true,geometry,
            PressureGaugePolicy::NONE,0,0.0,true);
        const std::array<std::vector<double>,3> rAU{{{0.2},{0.2},{0.2}}};
        ScalarBoundaryConditions pbcs;
        pbcs["x1"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        const auto full=make_rhie_chow_mass_flux(mesh,geometry,H,flux_p,rAU,1.2,
            VelocityBoundaryConditions{},pbcs,RhieChowVelocityContract::HBYA_PREDICTOR,
            false,false);
        const auto relaxed=make_rhie_chow_mass_flux(mesh,geometry,H,p,rAU,1.2,
            VelocityBoundaryConditions{},pbcs,RhieChowVelocityContract::HBYA_PREDICTOR,
            false,false);
        double largest=0.0;
        for (std::size_t f=0;f<mesh.n_faces();++f) {
            EXPECT_NEAR(full(f),10.0*relaxed(f),1e-12);
            largest=std::max(largest,std::abs(full(f)));
        }
        EXPECT_TRUE(largest>1e-6);
    });
    run_case("v19_pressure_flux_policy_refuses_incompatible_solver_modes", [] {
        IncompressibleSolverControls controls;
        EXPECT_TRUE(!controls.v19_full_pressure_flux_step);
        controls.v19_full_pressure_flux_step=true;
        EXPECT_THROW(validate_incompressible_controls(controls,1),std::invalid_argument);
        controls.v18_sparse_flux_pressure_jacobian=true;
        validate_incompressible_controls(controls,1);
        controls.pressure_continuity_line_search.enabled=true;
        EXPECT_THROW(validate_incompressible_controls(controls,1),std::invalid_argument);
        controls.pressure_continuity_line_search.enabled=false;
        controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        EXPECT_THROW(validate_incompressible_controls(controls,1),std::invalid_argument);
    });

    run_case("v20_flux_decomposition_separates_split_from_rc", [] {
        using namespace cfdx::core;
        const auto mesh=make_sst_transaction_unit_cube();
        Field<double,Location::FACE> full(mesh.n_faces(),"full","kg/s",1);
        Field<double,Location::FACE> relaxed(mesh.n_faces(),"relaxed","kg/s",1);
        Field<double,Location::FACE> from_u(mesh.n_faces(),"U","kg/s",1);
        full.fill(0.0); relaxed.fill(0.0); from_u.fill(0.0);
        full(0)=3.0; relaxed(0)=2.0; from_u(0)=0.5;
        full(1)=4.0; relaxed(1)=5.0; from_u(1)=2.0;
        const auto audit=audit_v20_flux_velocity_decomposition(
            mesh,full,relaxed,from_u);
        EXPECT_NEAR(audit.total_linf,2.5,1e-14);
        EXPECT_NEAR(audit.pressure_split_linf,1.0,1e-14);
        EXPECT_NEAR(audit.rhie_chow_linf,3.0,1e-14);
        EXPECT_NEAR(audit.total_l1,4.5,1e-14);
        EXPECT_NEAR(audit.pressure_split_l1,2.0,1e-14);
        EXPECT_NEAR(audit.rhie_chow_l1,4.5,1e-14);
        EXPECT_NEAR(audit.triangle_error_linf,0.0,1e-14);
        EXPECT_TRUE(audit.total_worst_boundary);
        // With no split, the mismatch is solely the physical collocated
        // face interpolation difference, not a missing transport flux.
        relaxed=full;
        const auto nosplit=audit_v20_flux_velocity_decomposition(
            mesh,full,relaxed,from_u);
        EXPECT_NEAR(nosplit.pressure_split_linf,0.0,1e-14);
        EXPECT_NEAR(nosplit.rhie_chow_linf,nosplit.total_linf,1e-14);
        from_u(0)=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(audit_v20_flux_velocity_decomposition(
            mesh,full,relaxed,from_u),std::domain_error);
    });

    run_case("v21_pressure_adaptation_uses_coupling_not_roundoff_continuity", [] {
        // V20 full-flux continuity is ~1e-12 even when momentum is large.
        EXPECT_NEAR(v21_pressure_relaxation_merit(true, 1e-12, .14, .02), .14, 1e-15);
        EXPECT_NEAR(v21_pressure_relaxation_merit(true, 1e-16, .03, .08), .08, 1e-15);
        EXPECT_NEAR(v21_pressure_relaxation_merit(false, .11, .03, .08), .11, 1e-15);
        EXPECT_NEAR(v21_pressure_relaxation_merit(false, .11,
                    std::numeric_limits<double>::quiet_NaN(), .08), .11, 1e-15);
        EXPECT_THROW(v21_pressure_relaxation_merit(true, 0.0, -1.0, 0.0),
                     std::invalid_argument);
        EXPECT_THROW(v21_pressure_relaxation_merit(true,
            std::numeric_limits<double>::quiet_NaN(), 0.1, 0.0), std::invalid_argument);
    });

    return run_all();
}
