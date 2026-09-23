#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/turbulence_transport.h"
#include "cfdx/physics/spalart_allmaras.h"
#include "cfdx/core/numerics/gradient.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct TurbulenceTransportResult {
    bool converged=false;
    std::size_t iterations=0;
    double k_residual=0.0;
    double second_residual=0.0;
};

inline void compute_kepsilon_production_field(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& epsilon,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& production,
    const TurbulenceTransportControls& c)
{
    for(std::size_t i=0;i<mesh.n_cells();++i) {
        const double nut=c.C_mu*std::max(k(i),c.k_min)*std::max(k(i),c.k_min)/
                         std::max(epsilon(i),c.epsilon_min);
        production(i)=c.density*2.0*nut*strain_rate(i)*strain_rate(i);
    }
}

inline TurbulenceTransportResult solve_kepsilon_transport(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& epsilon,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& epsilon_bcs = {},
    std::size_t max_iterations = 100,
    double tolerance = 1e-8)
{
    validate_turbulence_controls(controls);
    if(strain_rate.size()!=mesh.n_cells() || k.size()!=mesh.n_cells() ||
       epsilon.size()!=mesh.n_cells() || mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("k-epsilon field size mismatch");
    if(controls.model!=TurbulenceModel::KEPSILON)
        throw std::invalid_argument("solve_kepsilon_transport requires KEPSILON model");

    const std::size_t n=mesh.n_cells();
    cfdx::core::Field<double,cfdx::core::Location::CELL> P(
        n,"Pk","W/m3",1);
    TurbulenceTransportResult result;

    for(std::size_t iter=1;iter<=max_iterations;++iter) {
        auto oldk=k;
        auto olde=epsilon;
        compute_kepsilon_production_field(mesh,strain_rate,k,epsilon,P,controls);

        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(
            n,"Sk","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> se(
            n,"Se","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(
            n,"Spk","kg/m3/s",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spe(
            n,"Spe","kg/m3/s",1);

        for(std::size_t i=0;i<n;++i) {
            const double ki=std::max(k(i),controls.k_min);
            const double ei=std::max(epsilon(i),controls.epsilon_min);
            const double nut=controls.C_mu*ki*ki/ei;
            sk(i)=P(i);
            spk(i)=-controls.density*ei/ki;
            se(i)=controls.C1*P(i)*ei/ki;
            spe(i)=-controls.density*controls.C2*ei/ki;
        }

        std::vector<double> gamma_k(n), gamma_e(n);
        for(std::size_t i=0;i<n;++i) {
            const double ki=std::max(k(i),controls.k_min);
            const double ei=std::max(epsilon(i),controls.epsilon_min);
            const double nut=controls.C_mu*ki*ki/ei;
            gamma_k[i]=controls.density*
                (controls.molecular_viscosity+nut/controls.sigma_k);
            gamma_e[i]=controls.density*
                (controls.molecular_viscosity+nut/controls.sigma_epsilon);
        }

        auto eqk=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
            nullptr,nullptr,&gamma_k);
        auto eqe=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,se,spe,epsilon_bcs,true,
            nullptr,nullptr,&gamma_e);

        ScalarSolveControls sc;
        sc.max_iterations=2000;
        sc.tolerance=tolerance;
        sc.relaxation=0.7;
        cfdx::core::Vector k_solution(n,0.0);
        cfdx::core::Vector epsilon_solution(n,0.0);
        for(std::size_t i=0;i<n;++i) {
            k_solution(i)=k(i);
            epsilon_solution(i)=epsilon(i);
        }
        auto rk=solve_scalar_equation(eqk,k_solution,sc);
        auto re=solve_scalar_equation(eqe,epsilon_solution,sc);
        for(std::size_t i=0;i<n;++i) {
            k(i)=k_solution(i);
            epsilon(i)=epsilon_solution(i);
        }
        enforce_turbulence_bounds(k,epsilon,controls);

        result.k_residual=scalar_equation_residual_inf(eqk,k_solution);
        result.second_residual=scalar_equation_residual_inf(eqe,epsilon_solution);
        result.iterations=iter;

        double dk=0.0,de=0.0;
        for(std::size_t i=0;i<n;++i) {
            dk=std::max(dk,std::abs(k(i)-oldk(i)));
            de=std::max(de,std::abs(epsilon(i)-olde(i)));
        }
        if(rk.status==cfdx::core::SolverStatus::CONVERGED &&
           re.status==cfdx::core::SolverStatus::CONVERGED &&
           std::max(dk,de)<=tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}


inline TurbulenceTransportResult solve_komega_transport(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const TurbulenceTransportControls& controls,const ScalarBoundaryConditions& k_bcs={},
    const ScalarBoundaryConditions& omega_bcs={},std::size_t max_iterations=100,double tolerance=1e-8)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::KOMEGA) throw std::invalid_argument("solve_komega_transport requires KOMEGA model");
    const std::size_t n=mesh.n_cells();
    if(k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("k-omega field size mismatch");
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old_k=k,old_w=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1),spw(n,"Spw","kg/m3/s",1);
        std::vector<double> gk(n),gw(n);
        for(std::size_t i=0;i<n;++i){
            const double ki=std::max(k(i),controls.k_min), wi=std::max(omega(i),controls.omega_min);
            const double nut=controls.a1*ki/wi, P=2.0*nut*strain_rate(i)*strain_rate(i);
            sk(i)=P; spk(i)=-controls.density*controls.beta_star*wi;
            sw(i)=controls.gamma1*P/std::max(nut,1e-20); spw(i)=-controls.density*controls.beta1*wi;
            gk[i]=controls.density*(controls.molecular_viscosity+controls.sigma_k*nut);
            gw[i]=controls.density*(controls.molecular_viscosity+controls.sigma_epsilon*nut);
        }
        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,&gk);
        auto eqw=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,nullptr,nullptr,&gw);
        ScalarSolveControls sc{2000,tolerance,0.7}; cfdx::core::Vector ks(n,0.0),ws(n,0.0);
        for(std::size_t i=0;i<n;++i){ks(i)=k(i);ws(i)=omega(i);}
        const auto rk=solve_scalar_equation(eqk,ks,sc), rw=solve_scalar_equation(eqw,ws,sc);
        for(std::size_t i=0;i<n;++i){k(i)=ks(i);omega(i)=ws(i);}
        enforce_turbulence_bounds(k,omega,controls);
        double dk=0,dw=0; for(std::size_t i=0;i<n;++i){dk=std::max(dk,std::abs(k(i)-old_k(i)));dw=std::max(dw,std::abs(omega(i)-old_w(i)));}
        result.k_residual=scalar_equation_residual_inf(eqk,ks); result.second_residual=scalar_equation_residual_inf(eqw,ws); result.iterations=iter;
        if(rk.status==cfdx::core::SolverStatus::CONVERGED&&rw.status==cfdx::core::SolverStatus::CONVERGED&&std::max(dk,dw)<=tolerance){result.converged=true;break;}
    }
    return result;
}

inline TurbulenceTransportResult solve_spalart_allmaras_transport(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& nu_tilde,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& vorticity,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,const ScalarBoundaryConditions& bcs={},
    std::size_t max_iterations=100,double tolerance=1e-8)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::SPALART_ALLMARAS) throw std::invalid_argument("solve_spalart_allmaras_transport requires SPALART_ALLMARAS model");
    const std::size_t n=mesh.n_cells();
    if(nu_tilde.size()!=n||vorticity.size()!=n||wall_distance.size()!=n||mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("Spalart-Allmaras field size mismatch");
    SpalartAllmarasModel sa; sa.cb1=controls.sa_cb1;sa.cb2=controls.sa_cb2;sa.sigma=controls.sa_sigma;sa.kappa=controls.sa_kappa;sa.cw2=controls.sa_cw2;sa.cw3=controls.sa_cw3;sa.cv1=controls.sa_cv1;sa.ct3=controls.sa_ct3;sa.ct4=controls.sa_ct4;
    sa.cw1=sa.cb1/(sa.kappa*sa.kappa)+(1.0+sa.cb2)/sa.sigma;
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old=nu_tilde;
        cfdx::core::Field<double,cfdx::core::Location::CELL> su(n,"Snu","m2/s3",1),spu(n,"Spnu","1/s",1);
        std::vector<double> gamma(n);
        const auto grad_nu_tilde = cfdx::core::compute_gradient_gauss(nu_tilde,mesh);
        const double* gx = grad_nu_tilde.component_data(0);
        const double* gy = grad_nu_tilde.component_data(1);
        const double* gz = grad_nu_tilde.component_data(2);
        for(std::size_t i=0;i<n;++i){
            const double wt=std::max(nu_tilde(i),0.0), nu=controls.molecular_viscosity, d=wall_distance(i), vort=std::max(vorticity(i),1e-20);
            if(!(d>0)||!std::isfinite(d)||!std::isfinite(vort)) throw std::invalid_argument("invalid SA wall/vorticity field");
            const double chi=wt/nu, st=std::max(vort+wt*sa.fv2(chi)/(sa.kappa*sa.kappa*d*d),1e-20);
            const double r=wt/(st*sa.kappa*sa.kappa*d*d), fw=sa.destruction_coefficient(r), ft2=sa.ft2(chi);
            const double prod=sa.cb1*(1-ft2)*st;
            const double destr=std::max(sa.cw1*fw-sa.cb1*ft2/(sa.kappa*sa.kappa),0.0)/(d*d);
            const double grad2=gx[i]*gx[i]+gy[i]*gy[i]+gz[i]*gz[i];
            su(i)=controls.density*(prod*wt + sa.cb2/sa.sigma*grad2); spu(i)=-controls.density*destr;
            gamma[i]=controls.density*(nu+wt)/sa.sigma;
        }
        auto eq=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,su,spu,bcs,true,nullptr,nullptr,&gamma);
        ScalarSolveControls sc{2000,tolerance,0.7}; cfdx::core::Vector sol(n,0.0); for(std::size_t i=0;i<n;++i)sol(i)=nu_tilde(i);
        const auto rr=solve_scalar_equation(eq,sol,sc); for(std::size_t i=0;i<n;++i)nu_tilde(i)=std::max(sol(i),0.0);
        double du=0;for(std::size_t i=0;i<n;++i)du=std::max(du,std::abs(nu_tilde(i)-old(i)));
        result.k_residual=scalar_equation_residual_inf(eq,sol);result.second_residual=result.k_residual;result.iterations=iter;
        if(rr.status==cfdx::core::SolverStatus::CONVERGED&&du<=tolerance){result.converged=true;break;}
    }
    return result;
}

inline ScalarBoundaryCondition wall_k_epsilon(double k_value, double epsilon_value)
{
    if(k_value<0.0 || epsilon_value<=0.0)
        throw std::invalid_argument("invalid wall turbulence values");
    return {ScalarBoundaryType::FIXED_VALUE,k_value,0.0};
}

inline double wall_epsilon_from_k(double k, double y, double Cmu=0.09)
{
    if(k<0.0 || y<=0.0 || Cmu<=0.0)
        throw std::invalid_argument("invalid wall epsilon input");
    return std::pow(Cmu,0.75)*std::pow(std::max(k,0.0),1.5)/y;
}

inline double wall_omega_from_k(double k, double y, double betaStar=0.09)
{
    if(k<0.0 || y<=0.0 || betaStar<=0.0)
        throw std::invalid_argument("invalid wall omega input");
    return std::sqrt(std::max(k,0.0))/(std::sqrt(betaStar)*y);
}

} // namespace cfdx::physics
