#pragma once

#include "cfdx/physics/turbulence_solver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

inline std::pair<double,double> compute_sst_blending(
    double k, double omega, double wall_distance,
    double molecular_viscosity, double beta_star)
{
    if (!std::isfinite(k) || !std::isfinite(omega) ||
        !std::isfinite(wall_distance) || !std::isfinite(molecular_viscosity) ||
        !std::isfinite(beta_star) || beta_star <= 0.0 ||
        molecular_viscosity <= 0.0 || wall_distance <= 0.0)
        throw std::invalid_argument("compute_sst_blending: invalid inputs");

    const double ki = std::max(k, 0.0);
    const double wi = std::max(omega, 1e-20);
    const double y = std::max(wall_distance, 1e-12);
    const double arg1 = std::min(
        std::max(std::sqrt(ki) / (beta_star * wi * y),
                 500.0 * molecular_viscosity / (y*y*wi)),
        1.0e10);
    const double arg2 = std::max(
        2.0 * std::sqrt(ki) / (beta_star * wi * y),
        500.0 * molecular_viscosity / (y*y*wi));
    return {
        std::clamp(std::tanh(std::pow(arg1, 4.0)), 0.0, 1.0),
        std::clamp(std::tanh(arg2 * arg2), 0.0, 1.0)
    };
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
        auto oldk=k;
        auto oldw=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> P(
            n,"Pk","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(
            n,"Sk","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> sw(
            n,"Sw","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(
            n,"Spk","kg/m3/s",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spw(
            n,"Spw","kg/m3/s",1);
        std::vector<double> gamma_k(n),gamma_w(n);

        for(std::size_t i=0;i<n;++i) {
            const double ki=std::max(k(i),controls.k_min);
            const double wi=std::max(omega(i),controls.omega_min);
            const double f1=std::clamp(F1(i),0.0,1.0);
            const double f2=std::clamp(F2(i),0.0,1.0);
            const double beta=f1*controls.beta1+(1.0-f1)*controls.beta2;
            const double gamma=f1*(5.0/9.0)+(1.0-f1)*0.44;
            const double nut=controls.a1*ki/
                std::max(controls.a1*wi,strain_rate(i)*f2);
            P(i)=2.0*nut*strain_rate(i)*strain_rate(i);
            gamma_k[i]=controls.density*
                (controls.molecular_viscosity+controls.sigma_k*nut);
            gamma_w[i]=controls.density*
                (controls.molecular_viscosity+controls.sigma_epsilon*nut);
            sk(i)=P(i);
            spk(i)=-controls.density*controls.beta_star*wi;
            sw(i)=gamma*P(i)/std::max(nut,1e-20);
            spw(i)=-controls.density*beta*wi;
        }

        auto eqk=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
            nullptr,nullptr,&gamma_k);
        auto eqw=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
            nullptr,nullptr,&gamma_w);

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
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& omega_bcs = {},
    std::size_t max_iterations = 100,
    double tolerance = 1e-8)
{
    if(wall_distance.size()!=mesh.n_cells())
        throw std::invalid_argument("SST wall-distance field size mismatch");
    cfdx::core::Field<double,cfdx::core::Location::CELL> F1(mesh.n_cells(),"F1","1",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> F2(mesh.n_cells(),"F2","1",1);
    for(std::size_t i=0;i<mesh.n_cells();++i) {
        const auto blend=compute_sst_blending(k(i),omega(i),wall_distance(i),
                                              controls.molecular_viscosity,
                                              controls.beta_star);
        F1(i)=blend.first; F2(i)=blend.second;
    }
    return solve_sst_transport(mesh,geometry,mass_flux,k,omega,strain_rate,
                               F1,F2,controls,k_bcs,omega_bcs,max_iterations,tolerance);
}

} // namespace cfdx::physics
