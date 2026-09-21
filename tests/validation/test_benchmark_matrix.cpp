#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/numerics/divergence.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/numerics/laplacian.h"
#include "cfdx/core/numerics/temporal.h"
#include "cfdx/physics/equation_of_state.h"
#include "cfdx/physics/incompressible.h"
#include "cfdx/physics/m1_m4_models.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/turbulence.h"
#include "cfdx/transport/viscosity/viscosity.h"
#include "cfdx/transport/conductivity/conductivity.h"
#include "common/test_harness.h"

#include <cmath>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;
using namespace cfdx::physics::m1m4;
namespace transport = cfdx::transport;

namespace {
Mesh unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL; wall.face_ids={0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}
SparseMatrix spd_matrix()
{
    SparseMatrix A(3,3);
    A.push_back(0,0,4); A.push_back(0,1,1);
    A.push_back(1,0,1); A.push_back(1,1,3); A.push_back(1,2,1);
    A.push_back(2,1,1); A.push_back(2,2,2); A.finalize(); return A;
}
SparseMatrix nonsym_matrix()
{
    SparseMatrix A(3,3);
    A.push_back(0,0,4); A.push_back(0,1,1);
    A.push_back(1,0,2); A.push_back(1,1,3); A.push_back(1,2,1);
    A.push_back(2,1,1); A.push_back(2,2,2); A.finalize(); return A;
}
}

int main()
{
    // M1 — canonical analytical and numerical verification.
    run_case("M1_Couette_linear_profile", [] {
        const double H=2.0, U=3.0;
        EXPECT_NEAR(U*0.0/H,0.0,1e-14);
        EXPECT_NEAR(U*1.0/H,1.5,1e-14);
        EXPECT_NEAR(U*2.0/H,3.0,1e-14);
    });
    run_case("M1_Poiseuille_profile_and_flow_rate", [] {
        const double H=2.0, dpdx=-4.0, mu=2.0;
        EXPECT_NEAR(plane_poiseuille_velocity(1.0,H,dpdx,mu),1.0,1e-14);
        EXPECT_NEAR(plane_poiseuille_flow_rate_per_width(H,dpdx,mu),8.0/3.0,1e-14);
    });
    run_case("M1_Taylor_Green_2D_decay", [] {
        const double nu=0.01, t=0.7, k=1.0;
        const double decay=std::exp(-2.0*nu*k*k*t);
        EXPECT_NEAR(decay,0.9860975442628619,1e-14);
        EXPECT_NEAR(decay*decay,0.9723883668012469,1e-14);
    });
    run_case("M1_continuity_closed_cube_constant_velocity", [] {
        auto m=unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3); U(0,0)=2;U(0,1)=-1;U(0,2)=0.5;
        auto r=compute_continuity_residual(U,m); EXPECT_NEAR(r(0),0.0,1e-14);
    });

    // M2 — turbulence closures and canonical reference correlations.
    run_case("M2_k_epsilon_eddy_viscosity", [] {
        EXPECT_NEAR(turbulent_kinematic_viscosity_kepsilon(1.0,2.0),0.045,1e-14);
    });
    run_case("M2_SST_eddy_viscosity_limits", [] {
        EXPECT_NEAR(turbulent_kinematic_viscosity_komega_sst(0.5,2.0,0.0),0.25,1e-14);
        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(0.5,2.0,1.0)>0.0);
    });
    run_case("M2_Smagorinsky_scaling", [] {
        const double n1=smagorinsky_eddy_viscosity(0.1,4.0);
        const double n2=smagorinsky_eddy_viscosity(0.2,4.0);
        EXPECT_NEAR(n2/n1,4.0,1e-14);
    });
    run_case("M2_DES_length_scale_limiting", [] {
        EXPECT_NEAR(des_length_scale(1.0,0.1),0.1,1e-14);
        EXPECT_NEAR(des_length_scale(1.0,2.0),0.65,1e-14);
    });
    run_case("M2_wall_function_viscous_and_log_limits", [] {
        EXPECT_NEAR(wall_function_u_plus(5.0),5.0,1e-14);
        EXPECT_NEAR(wall_function_u_plus(100.0),std::log(9.793*100.0)/0.41,1e-14);
    });
    run_case("M2_flat_plate_laminar_Cf_correlation", [] {
        const double Re=1e6; EXPECT_NEAR(0.664/std::sqrt(Re),6.64e-4,1e-12);
    });
    run_case("M2_flat_plate_turbulent_Cf_correlation", [] {
        const double Re=1e6; EXPECT_NEAR(0.0592/std::pow(Re,0.2),0.003735,1e-6);
    });
    // M3 — heat transfer / CHT.
    run_case("M3_thermal_diffusivity", [] {
        EXPECT_NEAR(thermal_diffusivity(2.0,4.0,1000.0),5e-4,1e-14);
    });
    run_case("M3_two_layer_conduction_resistance", [] {
        const double q=(400.0-300.0)/(0.1/2.0+0.2/4.0);
        EXPECT_NEAR(q,1000.0,1e-12);
    });
    run_case("M3_interface_conductance_and_flux", [] {
        const double G=interface_conductance(2,4,0.1,0.2,1.0);
        EXPECT_NEAR(G,10.0,1e-14);
        EXPECT_NEAR(interface_heat_flux(G,400,300),1000.0,1e-12);
    });
    run_case("M3_fully_developed_channel_temperature_rise", [] {
        EXPECT_NEAR(fully_developed_temperature(2.0,300,1000,1,1000,2,1,4),304.0,1e-12);
    });
    run_case("M3_transient_energy_exponential_oracle", [] {
        const double alpha=2e-4, t=10.0, L=1.0;
        const double Fo=alpha*t/(L*L); EXPECT_NEAR(Fo,0.002,1e-14);
        EXPECT_NEAR(std::exp(-Fo),std::exp(-0.002),1e-14);
    });
    run_case("M3_radiative_conduction_balance", [] {
        const double qrad=blackbody(1000)-blackbody(500);
        EXPECT_TRUE(qrad>0.0);
        EXPECT_NEAR(qrad,-(blackbody(500)-blackbody(1000)),1e-12*qrad);
    });
    run_case("M3_source_linearization_exact_at_reference", [] {
        double Su=0,Sp=0; const double S=100-2*350;
        const double rebuilt=energy_source_linearization(S,-2,350,Su,Sp);
        EXPECT_NEAR(rebuilt,S,1e-14); EXPECT_NEAR(Su+Sp*350,S,1e-14);
    });

    // M4 — radiation canonical limits and enclosure identities.
    run_case("M4_Stefan_Boltzmann", [] {
        EXPECT_NEAR(blackbody(1000),56703.74419,1e-5);
    });
    run_case("M4_gray_surface_emission", [] {
        EXPECT_NEAR(gray_emission(0.8,1000),0.8*blackbody(1000),1e-10);
    });
    run_case("M4_black_two_surface_exchange", [] {
        EXPECT_NEAR(two_surface_exchange(1,1,1000,500,1),blackbody(1000)-blackbody(500),1e-8);
    });
    run_case("M4_gray_two_surface_exchange_symmetry", [] {
        const double q12=two_surface_exchange(0.7,0.4,1000,600,1);
        const double q21=two_surface_exchange(0.7,0.4,600,1000,1);
        EXPECT_NEAR(q12,-q21,1e-10*std::abs(q12));
    });
    run_case("M4_view_factor_closure", [] {
        validate_view_factors({0,1,1,0},2);
        validate_reciprocity({0,1,1,0},{1,1},2);
    });
    run_case("M4_P1_equilibrium_source", [] {
        EXPECT_NEAR(p1_source(2.0,p1_blackbody_intensity(900),900),0.0,1e-10*p1_blackbody_intensity(900));
    });
    run_case("M4_DOM_quadrature_closure", [] {
        std::vector<Direction> d={{1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},
                                  {0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},
                                  {0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};
        validate_dom(d);
    });
    run_case("M4_optically_thin_zero_absorption", [] {
        EXPECT_NEAR(p1_source(0,0,1000),0.0,1e-14);
    });

    // Numerics / infrastructure — manufactured and solver-level checks.
    run_case("N_gradient_constant_field", [] {
        auto m=unit_cube(); Field<double,Location::CELL> f(1,"f","1",1); f(0)=7;
        auto g=compute_gradient_gauss(f,m);
        EXPECT_NEAR(g(0,0),0,1e-14); EXPECT_NEAR(g(0,1),0,1e-14); EXPECT_NEAR(g(0,2),0,1e-14);
    });
    run_case("N_laplacian_constant_field", [] {
        auto m=unit_cube(); Field<double,Location::CELL> f(1,"f","1",1); f(0)=7;
        auto l=compute_laplacian(f,m); EXPECT_NEAR(l(0),0,1e-14);
    });
    run_case("N_CG_manufactured_linear_system", [] {
        auto A=spd_matrix(); Vector b(3,0),x(3,0); b(0)=6;b(1)=10;b(2)=8;
        auto r=solve_cg(A,b,x,100,1e-12); EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0),1.0,1e-10); EXPECT_NEAR(x(1),2.0,1e-10); EXPECT_NEAR(x(2),3.0,1e-10);
    });
    run_case("N_BiCGStab_manufactured_nonsymmetric_system", [] {
        auto A=nonsym_matrix(); Vector b(3,0),x(3,0); b(0)=6;b(1)=11;b(2)=8;
        auto r=solve_bicgstab(A,b,x,200,1e-12); EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0),1.0,1e-10); EXPECT_NEAR(x(1),2.0,1e-10); EXPECT_NEAR(x(2),3.0,1e-10);
    });
    run_case("N_GMRES_linear_operator_interface", [] {
        auto A=nonsym_matrix();
        LinearOperator op;
        op.size=3;
        op.apply=[&A](const Vector& in, Vector& out) {
            const auto y=A.matvec(in);
            for(std::size_t i=0;i<y.size();++i) out(i)=y[i];
        };
        Vector b(3,0),x(3,0); b(0)=6;b(1)=11;b(2)=8;
        auto r=solve_gmres(op,b,x,3,200,1e-12);
        EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0),1.0,1e-10);
        EXPECT_NEAR(x(1),2.0,1e-10);
        EXPECT_NEAR(x(2),3.0,1e-10);
    });
    run_case("N_GMRES_manufactured_nonsymmetric_system", [] {
        auto A=nonsym_matrix(); Vector b(3,0),x(3,0); b(0)=6;b(1)=11;b(2)=8;
        auto r=solve_gmres(A,b,x,3,200,1e-12); EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0),1.0,1e-10); EXPECT_NEAR(x(1),2.0,1e-10); EXPECT_NEAR(x(2),3.0,1e-10);
    });
    run_case("N_explicit_Euler_temporal_oracle", [] {
        Field<double,Location::CELL> phi(1,"phi","1",1); phi(0)=1;
        const double lambda=2,dt=1e-3;
        const auto rhs=[&](const auto& in,auto& out){out(0)=-lambda*in(0);};
        auto next=advance_time(phi,dt,rhs,TimeScheme::EULER_EXPLICIT);
        EXPECT_NEAR(next(0),1-lambda*dt,1e-14);
    });
    run_case("N_Crank_Nicolson_temporal_oracle", [] {
        Field<double,Location::CELL> phi(1,"phi","1",1); phi(0)=1;
        const double lambda=2,dt=1e-2;
        const auto rhs=[&](const auto& in,auto& out){out(0)=-lambda*in(0);};
        auto next=advance_time(phi,dt,rhs,TimeScheme::CRANK_NICOLSON);
        const double exact=(1-lambda*dt/2)/(1+lambda*dt/2);
        EXPECT_NEAR(next(0),exact,1e-10);
    });
    run_case("N_EOS_ideal_gas_round_trip", [] {
        IdealGasEOS eos; const double p=101325,T=300; const double rho=eos.density(p,T);
        EXPECT_NEAR(eos.pressure_from_density_temp(rho,T),p,1e-10);
        EXPECT_NEAR(eos.speed_of_sound(p,T),
                    std::sqrt(1.4*(8.314462618/0.02896546)*300),1e-10);
    });
    run_case("N_Sutherland_transport_reference", [] {
        const double mu=transport::viscosity::sutherland(273.15);
        EXPECT_NEAR(mu,1.716e-5,1e-12);
        EXPECT_NEAR(transport::conductivity::prandtl(1e-3,1000,1),1.0,1e-14);
    });
    run_case("N_momentum_zero_state_fixed_point", [] {
        auto m=unit_cube(); Field<double,Location::CELL> U(1,"U","m/s",3),p(1,"p","Pa",1);
        U.fill(0); p.fill(0);
        auto terms=evaluate_momentum_terms(U,p,m,1e-3);
        EXPECT_NEAR(terms.convection(0,0),0,1e-14);
        EXPECT_NEAR(terms.diffusion(0,0),0,1e-14);
        EXPECT_NEAR(terms.pressure_gradient(0,0),0,1e-14);
    });
    run_case("N_turbulence_production_zero_strain", [] {
        Vec3 s{0,0,0}; EXPECT_NEAR(turbulence_production(s,0.5),0,1e-14);
    });

    return run_all();
}
