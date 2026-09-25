#pragma once

#include "cfdx/physics/turbulence_solver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

inline std::pair<double,double> compute_sst_blending(
    double k, double omega, double wall_distance,
    double molecular_viscosity, double beta_star,
    double grad_k_dot_grad_omega, double density, double sigma_w2)
{
    if(!std::isfinite(grad_k_dot_grad_omega) || !(density>0.0) || !(sigma_w2>0.0))
        throw std::invalid_argument("compute_sst_blending: invalid cross-diffusion inputs");
    const double ki=std::max(k,0.0), wi=std::max(omega,1e-20), y=std::max(wall_distance,1e-12);
    const double cd_kw=std::max(2.0*density*sigma_w2/wi*grad_k_dot_grad_omega,1e-10);
    const double arg1=std::min(
        std::max(std::sqrt(ki)/(beta_star*wi*y),500.0*molecular_viscosity/(y*y*wi)),
        4.0*density*sigma_w2*ki/(cd_kw*y*y));
    const double arg2=std::max(2.0*std::sqrt(ki)/(beta_star*wi*y),
                                500.0*molecular_viscosity/(y*y*wi));
    return {std::clamp(std::tanh(std::pow(arg1,4.0)),0.0,1.0),
            std::clamp(std::tanh(arg2*arg2),0.0,1.0)};
}

inline std::pair<double,double> compute_sst_blending(
    double k, double omega, double wall_distance,
    double molecular_viscosity, double beta_star)
{
    return compute_sst_blending(k,omega,wall_distance,molecular_viscosity,beta_star,
                                0.0,1.0,0.856);
}

struct SSTCellSources {
    double sk=0.0, spk=0.0, sw=0.0, spw=0.0;
    double gamma_k=0.0, gamma_w=0.0, nut=0.0;
};

inline SSTCellSources sst_cell_sources(
    double k, double omega, double strain, double F1, double F2,
    double grad_k_dot_grad_omega, const TurbulenceTransportControls& c)
{
    const double ki=std::max(k,c.k_min), wi=std::max(omega,c.omega_min);
    const double S=std::max(strain,0.0), f1=std::clamp(F1,0.0,1.0), f2=std::clamp(F2,0.0,1.0);
    const auto blend=[f1](double a,double b){return f1*a+(1.0-f1)*b;};
    SSTCellSources s;
    s.nut=c.a1*ki/std::max(c.a1*wi,S*f2);
    const double Pk=std::min(c.density*s.nut*S*S,
                             c.sst_production_limiter*c.beta_star*c.density*ki*wi);
    s.sk=Pk;
    s.spk=-c.density*c.beta_star*wi;
    const double cross=2.0*(1.0-f1)*c.density*c.sst_sigma_w2/wi*grad_k_dot_grad_omega;
    s.sw=blend(c.gamma1,c.gamma2)*c.density*S*S+std::max(cross,0.0);
    s.spw=-c.density*blend(c.beta1,c.beta2)*wi+std::min(cross,0.0)/wi;
    s.gamma_k=c.density*(c.molecular_viscosity+
        blend(c.sst_sigma_k1,c.sst_sigma_k2)*s.nut);
    s.gamma_w=c.density*(c.molecular_viscosity+
        blend(c.sst_sigma_w1,c.sst_sigma_w2)*s.nut);
    return s;
}

inline TurbulenceTransportResult solve_sst_transport(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& F1,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& F2,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& omega_bcs = {},
    std::size_t max_iterations = 100,
    double tolerance = 1e-8)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::SST)
        throw std::invalid_argument("solve_sst_transport requires SST model");
    const std::size_t n=mesh.n_cells();
    if(F1.size()!=n || F2.size()!=n || strain_rate.size()!=n)
        throw std::invalid_argument("SST blending field size mismatch");

    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter) {
        const auto cross = cell_gradient_dot(k,omega,mesh);
        auto oldk=k;
        auto oldw=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1),spw(n,"Spw","kg/m3/s",1);
        std::vector<double> gamma_k(n),gamma_w(n);
        for(std::size_t i=0;i<n;++i) {
            const auto q=sst_cell_sources(k(i),omega(i),strain_rate(i),F1(i),F2(i),cross[i],controls);
            sk(i)=q.sk; spk(i)=q.spk; sw(i)=q.sw; spw(i)=q.spw;
            gamma_k[i]=q.gamma_k; gamma_w[i]=q.gamma_w;
        }

        auto eqk=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
            nullptr,nullptr,nullptr,&gamma_k);
        auto eqw=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
            nullptr,nullptr,nullptr,&gamma_w);

        ScalarSolveControls sc{2000,tolerance,0.7};
        cfdx::core::Vector k_solution(n, 0.0);
        cfdx::core::Vector omega_solution(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            k_solution(i) = k(i);
            omega_solution(i) = omega(i);
        }
        const auto rk = solve_scalar_equation(eqk, k_solution, sc);
        const auto rw = solve_scalar_equation(eqw, omega_solution, sc);
        for (std::size_t i = 0; i < n; ++i) {
            k(i) = k_solution(i);
            omega(i) = omega_solution(i);
        }
        enforce_turbulence_bounds(k,omega,controls);

        double dk=0.0,dw=0.0;
        for(std::size_t i=0;i<n;++i) {
            dk=std::max(dk,std::abs(k(i)-oldk(i)));
            dw=std::max(dw,std::abs(omega(i)-oldw(i)));
        }
        result.k_residual=scalar_equation_residual_inf(eqk,k_solution);
        result.second_residual=scalar_equation_residual_inf(eqw,omega_solution);
        result.iterations=iter;
        if(rk.status==cfdx::core::SolverStatus::CONVERGED &&
           rw.status==cfdx::core::SolverStatus::CONVERGED &&
           std::max(dk,dw)<=tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

inline TurbulenceTransportResult solve_sst_transport_dynamic_blending(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,const ScalarBoundaryConditions& k_bcs={},
    const ScalarBoundaryConditions& omega_bcs={},std::size_t max_iterations=100,double tolerance=1e-8)
{
    validate_turbulence_controls(controls);
    const std::size_t n=mesh.n_cells();
    if(wall_distance.size()!=n||k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("SST dynamic-blending field size mismatch");
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old_k=k,old_w=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> F1(n,"F1","1",1),F2(n,"F2","1",1);
        for(std::size_t i=0;i<n;++i){
            const auto cross=cell_gradient_dot(k,omega,mesh);
            const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
            F1(i)=b.first; F2(i)=b.second;
        }
        const auto inner=solve_sst_transport(mesh,geometry,mass_flux,k,omega,strain_rate,F1,F2,controls,k_bcs,omega_bcs,1,tolerance);
        double dk=0.0,dw=0.0;
        for(std::size_t i=0;i<n;++i){dk=std::max(dk,std::abs(k(i)-old_k(i)));dw=std::max(dw,std::abs(omega(i)-old_w(i)));}
        result.k_residual=inner.k_residual; result.second_residual=inner.second_residual; result.iterations=iter;
        if(inner.converged&&std::max(dk,dw)<=tolerance){result.converged=true;break;}
    }
    return result;
}

} // namespace cfdx::physics
