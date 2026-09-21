#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/turbulence_transport.h"
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
            const double gamma_k=controls.density*
                (controls.molecular_viscosity+nut/controls.sigma_k);
            const double gamma_e=controls.density*
                (controls.molecular_viscosity+nut/controls.sigma_epsilon);
            (void)gamma_k; (void)gamma_e;
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
        sc.relaxation=1.0;
        cfdx::core::Vector k_solution(n,0.0);
        cfdx::core::Vector epsilon_solution(n,0.0);
        for(std::size_t i=0;i<n;++i) {
            k_solution(i)=k(i);
            epsilon_solution(i)=epsilon(i);
        }
        auto rk=solve_scalar_equation(eqk,k_solution,sc);
        auto re=solve_scalar_equation(eqe,epsilon_solution,sc);
        for(std::size_t i=0;i<n;++i) {
            k(i)=k(i)+0.7*(k_solution(i)-k(i));
            epsilon(i)=epsilon(i)+0.7*(epsilon_solution(i)-epsilon(i));
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
