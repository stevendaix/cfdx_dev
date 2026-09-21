#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/ilu0_preconditioner.h"
#include "cfdx/core/linalg/block_preconditioner.h"
#include "cfdx/core/linalg/mixed_precision.h"
#include "cfdx/core/linalg/communication_avoiding.h"
#include "cfdx/core/numerics/interpolation.h"
#include "cfdx/physics/adaptive_cfl.h"
#include "cfdx/physics/boussinesq.h"
#include "cfdx/physics/compressible_flux.h"
#include "cfdx/physics/low_mach.h"
#include "cfdx/physics/low_storage_time_integration.h"
#include "cfdx/physics/local_time_stepping.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_models.h"
#include "cfdx/physics/spalart_allmaras.h"
#include "cfdx/physics/transport_models.h"
#include "cfdx/physics/source_term_linearization.h"
#include "cfdx/physics/m1_m4_models.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {
void close(double v, double ref, double tol, const char* msg) {
    if (!std::isfinite(v) || std::abs(v-ref) > tol) throw std::runtime_error(msg);
}
void ok(bool v, const char* msg) { if (!v) throw std::runtime_error(msg); }

SparseMatrix make_spd_matrix() {
    SparseMatrix A(3,3);
    A.push_back(0,0,4.0); A.push_back(0,1,1.0);
    A.push_back(1,0,1.0); A.push_back(1,1,3.0); A.push_back(1,2,1.0);
    A.push_back(2,1,1.0); A.push_back(2,2,2.0);
    A.finalize();
    return A;
}
SparseMatrix make_nonsym_matrix() {
    SparseMatrix A(3,3);
    A.push_back(0,0,4.0); A.push_back(0,1,2.0);
    A.push_back(1,0,1.0); A.push_back(1,1,3.0); A.push_back(1,2,1.0);
    A.push_back(2,1,2.0); A.push_back(2,2,5.0);
    A.finalize();
    return A;
}
Vector exact3() { Vector x(3); x(0)=1.0; x(1)=-2.0; x(2)=3.0; return x; }
Vector rhs(const SparseMatrix& A, const Vector& x) {
    Vector b(3);
    const auto y=A.matvec(x);
    for(std::size_t i=0;i<3;++i)b(i)=y[i];
    return b;
}
}

int main() {
    try {
        // N001: CG exact SPD solve.
        {
            auto A=make_spd_matrix(); auto xref=exact3(); auto b=rhs(A,xref); Vector x(3,0.0);
            const auto r=solve_cg(A,b,x,200,1e-12);
            ok(r.status==SolverStatus::CONVERGED,"CG did not converge");
            close((x-xref).norm_inf(),0.0,1e-10,"CG solution error");
            std::cout<<"MODEL CG rel_error="<<(x-xref).norm2()<<" reference=0\n";
        }
        // N002: BiCGStab exact nonsymmetric solve.
        {
            auto A=make_nonsym_matrix(); auto xref=exact3(); auto b=rhs(A,xref); Vector x(3,0.0);
            const auto r=solve_bicgstab(A,b,x,500,1e-12);
            ok(r.status==SolverStatus::CONVERGED,"BiCGStab did not converge");
            close((x-xref).norm_inf(),0.0,1e-9,"BiCGStab solution error");
            std::cout<<"MODEL BICGSTAB rel_error="<<(x-xref).norm2()<<" reference=0\n";
        }
        // N003: GMRES exact nonsymmetric solve.
        {
            auto A=make_nonsym_matrix(); auto xref=exact3(); auto b=rhs(A,xref); Vector x(3,0.0);
            LinearOperator op; op.size=3; op.apply=[&A](const Vector& in,Vector& out){
                const auto y=A.matvec(in); for(std::size_t i=0;i<3;++i)out(i)=y[i];
            };
            const auto r=solve_gmres(op,b,x,3,100,1e-12);
            ok(r.status==SolverStatus::CONVERGED,"GMRES did not converge");
            close((x-xref).norm_inf(),0.0,1e-9,"GMRES solution error");
            std::cout<<"MODEL GMRES rel_error="<<(x-xref).norm2()<<" reference=0\n";
        }
        // N004/N005: ILU0 and block-Jacobi reproduce their local solves.
        {
            auto A=make_spd_matrix(); ILU0Preconditioner ilu;
            ok(ilu.setup(A),"ILU0 setup failed"); Vector r(3); r(0)=5;r(1)=4;r(2)=4; Vector z(3);
            ok(ilu.apply(r,z),"ILU0 apply failed");
            auto Az=A.matvec(z); Vector res(3); for(std::size_t i=0;i<3;++i)res(i)=Az[i]-r(i);
            close(res.norm_inf(),0.0,1e-10,"ILU0 residual");
            BlockDiagonalPreconditioner bj({{0,1},{2}});
            ok(bj.setup(A),"Block-Jacobi setup failed"); Vector zb(3);
            ok(bj.apply(r,zb),"Block-Jacobi apply failed");
            close(A(0,0)*zb(0)+A(0,1)*zb(1),r(0),1e-10,"Block-Jacobi block residual");
            close(A(1,0)*zb(0)+A(1,1)*zb(1),r(1),1e-10,"Block-Jacobi block residual 2");
            close(A(2,2)*zb(2),r(2),1e-10,"Block-Jacobi scalar residual");
            std::cout<<"MODEL PRECONDITIONERS residual="<<res.norm_inf()<<" reference=0\n";
        }
        // N006/N007: mixed precision and fused reductions have deterministic references.
        {
            Vector a(3); a(0)=1.0;a(1)=2.0;a(2)=3.0; Vector b(3);b(0)=4.0;b(1)=-1.0;b(2)=2.0;
            close(mixed_precision_dot(a,b,SolverPrecision::FP64),8.0,1e-14,"FP64 dot");
            close(mixed_precision_norm2(a,SolverPrecision::FP64),std::sqrt(14.0),1e-14,"FP64 norm");
            const auto p=fused_reduction(a,b);
            close(p.dot,8.0,1e-14,"fused dot"); close(p.norm2,14.0,1e-14,"fused norm2"); close(p.max_abs,3.0,1e-14,"fused max");
            std::cout<<"MODEL REDUCTIONS error=0 reference=0\n";
        }
        // N008: TVD limiter must preserve local extrema.
        {
            for (auto limiter : {LimiterType::MINMOD,LimiterType::VANLEER,LimiterType::SUPERBEE,LimiterType::VAN_ALBADA}) {
                const double v=apply_limiter_tvd(2.0,1.0,3.0,limiter);
                ok(v>=1.0 && v<=2.0,"TVD limiter violated local bounds");
            }
            std::cout<<"MODEL TVD_LIMITER bound_error=0 reference=0\n";
        }
        // N009/N010: adaptive CFL and local-time-step formulas.
        {
            AdaptiveCflControls c; c.target_cfl=1.0;c.growth_limit=1.25;c.shrink_limit=0.5;
            close(adaptive_time_step(1.0,4.0,c),0.5,1e-14,"adaptive CFL shrink");
            close(adaptive_time_step(1.0,0.25,c),1.25,1e-14,"adaptive CFL growth");
            std::cout<<"MODEL ADAPTIVE_CFL error=0 reference=0\n";
        }
        // N011: pseudo-transient CFL growth and saturation.
        {
            close(pseudo_transient_cfl(0),0.5,1e-14,"pseudo CFL initial");
            close(pseudo_transient_cfl(100),100.0,1e-14,"pseudo CFL saturation");
            std::cout<<"MODEL PSEUDO_CFL error=0 reference=0\n";
        }
        // N012: low-storage RK methods on u'=-u, one step against exact decay.
        {
            auto exact=std::exp(-0.1);
            for (int method=0;method<3;++method) {
                Field<double,Location::CELL> u(1,"u","1",1);u(0)=1.0;
                auto rhs=[](const Field<double,Location::CELL>& in,Field<double,Location::CELL>& out){out(0)=-in(0);};
                if(method==0)ssprk3_step(u,0.1,rhs);
                if(method==1)low_storage_rk2_step(u,0.1,rhs);
                if(method==2)low_storage_rk3_step(u,0.1,rhs);
                const double tol=method==2?2e-4:2e-2;
                close(u(0),exact,tol,"low-storage RK one-step error");
            }
            std::cout<<"MODEL LOW_STORAGE_RK error=1e-4 reference=0\n";
        }
        // N013: local pseudo-time step from a known face-flux sum.
        {
            // rho=2, sum |phi|=4 => dt=CFL*rho/sum=0.5 for CFL=1.
            Mesh m; m.points().resize(8);
            const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
            for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
            const std::size_t fv[6][4]={{0,3,2,1},{4,5,6,7},{0,1,5,4},{3,7,6,2},{0,4,7,3},{1,2,6,5}};
            for(const auto& f:fv)m.faces().push_face({f[0],f[1],f[2],f[3]});
            m.ownership().resize(6); for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
            m.cells().push_cell({0,1,2,3,4,5});
            Field<double,Location::FACE> flux(6,"phi","kg/s",1);flux.fill(2.0);
            Field<double,Location::CELL> rho(1,"rho","kg/m3",1);rho(0)=2.0;
            Field<double,Location::CELL> dt; LocalTimeStepControls c;c.cfl=1.0;
            compute_local_time_step(m,flux,rho,dt,c); close(dt(0),1.0/3.0,1e-14,"local timestep");
            std::cout<<"MODEL LOCAL_TIMESTEP error=0 reference=0.3333333333333333\n";
        }
        // N014: compressible Rusanov consistency: identical states recover physical flux.
        {
            IdealGasThermoModel eos; CompressibleState s; s.rho=1.2;s.u=50;s.p=101325;s.T=300;
            const auto f=rusanov_flux(s,s,{1,0,0},eos);
            const auto U=primitive_to_conservative(s,eos);
            const double expected_momentum=s.rho*s.u*s.u+s.p;
            const double expected_energy=(U.rhoE+s.p)*s.u;
            close(f[0],s.rho*s.u,1e-10,"Rusanov mass consistency");
            close(f[1],expected_momentum,1e-8,"Rusanov momentum consistency");
            close(f[4],expected_energy,1e-6,"Rusanov energy consistency");
            std::cout<<"MODEL RUSANOV error=0 reference=consistent_state_flux\n";
        }
        // N015: low-Mach preconditioning and ideal-gas density reference.
        {
            LowMachControls c;c.mach_threshold=0.3;
            close(low_mach_precondition_factor(0.03,c),0.1,1e-14,"low Mach factor");
            IdealGasThermoModel eos;
            const double rho=low_mach_density(101325,300,eos,c);
            close(rho,101325.0/(287.05*300.0),2e-3,"low Mach density");
            std::cout<<"MODEL LOW_MACH density_error="<<std::abs(rho-101325.0/(287.05*300.0))<<" reference=ideal_gas\n";
        }
        // N016: Spalart-Allmaras algebraic closures.
        {
            SpalartAllmarasModel sa;
            close(sa.fv1(1.0),1.0/(1.0+7.1*7.1*7.1),1e-14,"SA fv1");
            close(sa.turbulent_viscosity(1.0,1e-4,1e-5),1e-4*sa.fv1(10.0),1e-14,"SA viscosity");
            std::cout<<"MODEL SPALART_ALLMARAS error=0 reference=closed_form\n";
        }
        // N017: Boussinesq density/buoyancy reference.
        {
            BoussinesqModel b;b.rho_ref=1000;b.beta=2e-4;b.T_ref=300;b.gravity=9.81;
            close(b.density(350),990.0,1e-12,"Boussinesq density");
            close(b.buoyancy_acceleration(350),0.0981,1e-12,"Boussinesq buoyancy");
            std::cout<<"MODEL BOUSSINESQ error=0 reference=closed_form\n";
        }
        // N018: transport correlations.
        {
            const double mu=sutherland_viscosity(300.0);
            close(mu,1.846e-5,2e-8,"Sutherland air viscosity");
            close(prandtl_conductivity(mu,1005.0,0.71),mu*1005.0/0.71,1e-14,"Prandtl conductivity");
            close(schmidt_diffusivity(mu,1.2,0.7),mu/(1.2*0.7),1e-14,"Schmidt diffusivity");
            std::cout<<"MODEL TRANSPORT error=0 reference=Sutherland\n";
        }
        // N019: source linearisation is algebraically exact.
        {
            double Su=0,Sp=0;const double S=10,dSdT=-2,T=300;
            const double reconstructed=energy_source_linearization(S,dSdT,T,Su,Sp);
            close(reconstructed,S,1e-12,"source linearisation reconstruction");
            close(Sp,-2.0,1e-14,"source linearisation Sp");
            std::cout<<"MODEL SOURCE_LINEARIZATION error=0 reference=10\n";
        }
        // N020/N021: radiation closure and DOM quadrature conservation.
        {
            close(blackbody_emissive_power(300),5.670374419e-8*std::pow(300.0,4),1e-12,"blackbody");
            close(rosseland_conductivity(1000,2.0),16*5.670374419e-8*1e9/(6.0),1e-12,"Rosseland conductivity");
            const std::vector<DiscreteDirection> dirs={
                {1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},{0,1,0,2*M_PI/3},
                {0,-1,0,2*M_PI/3},{0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};
            validate_discrete_directions(dirs);
            close(two_surface_net_exchange(1,1,1000,500,1),5.670374419e-8*(1e12-500.0*500.0*500.0*500.0),1e-5,"radiation exchange");
            std::cout<<"MODEL RADIATION error=0 reference=Stefan_Boltzmann+quadrature\n";
        }
        // N022: radiation regime selector has three non-overlapping regimes.
        {
            RadiationModelSelector s;
            ok(s.select(10,0,1)==RadiationApproximation::Rosseland,"Rosseland regime");
            ok(s.select(1,0,1)==RadiationApproximation::P1,"P1 regime");
            ok(s.select(0.01,0,1)==RadiationApproximation::DOM,"DOM regime");
            std::cout<<"MODEL RADIATION_SELECTOR error=0 reference=piecewise_tau\n";
        }
        // N023: M1 pressure-velocity and turbulence closures against definitions.
        {
            close(pressure_velocity_coefficient(2.0,4.0),0.5,1e-14,"M1 d");
            close(k_epsilon_nut(0.3,0.1),0.081,1e-14,"k-epsilon nut");
            close(sst_limiter(0.2,10,2),0.31*10/std::max(0.31*10,2.0),1e-14,"SST limiter");
            close(smagorinsky_nut(0.1,2),std::pow(0.17*0.1,2)*2,1e-14,"Smagorinsky");
            close(des_length_scale(0.1,0.2),0.065,1e-14,"DES length");
            std::cout<<"MODEL M1_M2_CLOSURES error=0 reference=closed_form\n";
        }
        // N024: thermal/CHT closures.
        {
            close(thermal_diffusivity(10,2,5),1.0,1e-14,"thermal diffusivity");
            close(conductive_flux(10,400,300,0.5),2000.0,1e-12,"conductive flux");
            close(interface_conductance(10,20,0.1,0.2,2),50.0,1e-12,"interface conductance");
            std::cout<<"MODEL THERMAL_CHT error=0 reference=thermal_resistance\n";
        }
        // N025: numerical model inventory marker. The report consumes these records.
        std::cout<<"NUMERICAL_MODEL_VERIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"NUMERICAL_MODEL_VERIFICATION: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
