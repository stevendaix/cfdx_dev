#pragma once

#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/scalar_transport_gradient.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

inline std::vector<double> sst_boundary_gradient_dot(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
{
    const auto grad_k = scalar_gradient_with_transport_boundaries(k, mesh, geometry, k_bcs);
    const auto grad_omega = scalar_gradient_with_transport_boundaries(
        omega, mesh, geometry, omega_bcs, omega_face_values);
    std::vector<double> result(mesh.n_cells());
    for (std::size_t cell = 0; cell < result.size(); ++cell) {
        result[cell] = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            result[cell] += grad_k(cell, component) * grad_omega(cell, component);
        if (!std::isfinite(result[cell]))
            throw std::invalid_argument("SST boundary gradient contraction is non-finite");
    }
    return result;
}

// Residual in solution units. Scaling each row by sum_j |a_ij| makes this
// comparable with the existing absolute iterate-change tolerance. Unlike a
// RHS-relative residual, it remains meaningful for a homogeneous equation.
// The unscaled physical residual is still reported in the transport result.
inline double sst_equation_solution_defect(
    const ScalarEquation& equation, const cfdx::core::Vector& solution)
{
    if (solution.size() != equation.rhs.size())
        throw std::invalid_argument("SST residual: field size mismatch");
    double maximum = 0.0;
    for (std::size_t cell = 0; cell < solution.size(); ++cell) {
        double residual = -equation.rhs(cell), scale = 0.0;
        for (auto entry = equation.matrix.row_offsets_data()[cell];
             entry < equation.matrix.row_offsets_data()[cell + 1]; ++entry) {
            const double coefficient = equation.matrix.values_data()[entry];
            residual += coefficient * solution(equation.matrix.columns_data()[entry]);
            scale += std::abs(coefficient);
        }
        if (!std::isfinite(residual) || !std::isfinite(scale) || !(scale > 0.0))
            return std::numeric_limits<double>::infinity();
        maximum = std::max(maximum, std::abs(residual) / scale);
    }
    return maximum;
}

// Translate the absolute solution-unit gate into a sufficient relative
// Krylov tolerance. Otherwise a warm start can satisfy ||r||2/||b||2 while
// failing the row-scaled SST gate, so every solve returns the same rejected
// iterate forever. This only tightens the requested linear tolerance; the
// external tolerance, iteration budgets and relaxation are unchanged.
inline double sst_linear_tolerance_for_solution_defect(
    const ScalarEquation& equation, double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance<=0.0 || equation.rhs.size()==0)
        throw std::invalid_argument("SST linear tolerance: invalid requested tolerance or empty equation");
    double rhs_norm=0.0;
    double minimum_row_norm=std::numeric_limits<double>::infinity();
    for (std::size_t cell=0; cell<equation.rhs.size(); ++cell) {
        if (!std::isfinite(equation.rhs(cell)))
            throw std::invalid_argument("SST linear tolerance: non-finite right-hand side");
        rhs_norm=std::hypot(rhs_norm,equation.rhs(cell));
        double row_norm=0.0;
        for (auto entry=equation.matrix.row_offsets_data()[cell];
             entry<equation.matrix.row_offsets_data()[cell+1]; ++entry)
            row_norm+=std::abs(equation.matrix.values_data()[entry]);
        if (!std::isfinite(row_norm) || !(row_norm>0.0))
            throw std::invalid_argument("SST linear tolerance: zero or non-finite matrix row norm");
        minimum_row_norm=std::min(minimum_row_norm,row_norm);
    }
    if (!std::isfinite(rhs_norm))
        throw std::invalid_argument("SST linear tolerance: right-hand side norm overflow");
    // Native CG/BiCGStab use max(||b||2,1e-15). GMRES/dense fallbacks use
    // smaller RHS scales, so this is sufficient for those paths as well.
    const double ratio=std::min(1.0,minimum_row_norm/std::max(rhs_norm,1e-15));
    const double translated=tolerance*ratio;
    if (!std::isfinite(translated) || !(translated>0.0))
        throw std::domain_error("SST linear tolerance: requested absolute defect cannot be represented");
    return translated;
}

inline std::pair<double,double> compute_sst_blending(
    double k, double omega, double wall_distance,
    double molecular_viscosity, double beta_star,
    double grad_k_dot_grad_omega, double density, double sigma_w2)
{
    if(!std::isfinite(k) || !std::isfinite(omega) || !std::isfinite(wall_distance) ||
       !std::isfinite(molecular_viscosity) || !std::isfinite(beta_star) ||
       !std::isfinite(grad_k_dot_grad_omega) ||
       !(k >= 0.0) || !(omega > 0.0) || !(wall_distance > 0.0) ||
       !(molecular_viscosity >= 0.0) || !(beta_star > 0.0) ||
       !(density > 0.0) || !(sigma_w2 > 0.0))
        throw std::invalid_argument("compute_sst_blending: invalid inputs");
    const double ki=k, wi=omega, y=wall_distance;
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
    // SST-2003 uses the limited k-production in the omega production term.
    // If nut is zero, the limited production is also zero; avoid a 0/0 path.
    const double omega_production = s.nut > 0.0 ? Pk / s.nut : 0.0;
    s.sw=blend(c.gamma1,c.gamma2)*omega_production+std::max(cross,0.0);
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
    double tolerance = 1e-8,
    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::SST)
        throw std::invalid_argument("solve_sst_transport requires SST model");
    const std::size_t n=mesh.n_cells();
    if(F1.size()!=n || F2.size()!=n || strain_rate.size()!=n ||
       k.size()!=n || omega.size()!=n || mass_flux.size()!=mesh.n_faces() ||
       k.dimension()!=1 || omega.dimension()!=1 || strain_rate.dimension()!=1 ||
       F1.dimension()!=1 || F2.dimension()!=1 || mass_flux.dimension()!=1)
        throw std::invalid_argument("SST blending field size mismatch");
    if (max_iterations == 0 || !std::isfinite(tolerance) || tolerance <= 0.0)
        throw std::invalid_argument("SST transport requires valid iteration controls");
    if (omega_face_values)
        validate_transport_boundary_face_values(mesh, geometry, *omega_face_values, true);

    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter) {
        const auto cross = sst_boundary_gradient_dot(
            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
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
            nullptr,nullptr,nullptr,&gamma_k, ConvectionScheme::UPWIND, &k);
        auto eqw=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
            omega_face_values,nullptr,nullptr,&gamma_w, ConvectionScheme::UPWIND, &omega);

        const ScalarSolveControls k_sc{
            2000,sst_linear_tolerance_for_solution_defect(eqk,tolerance),0.7};
        const ScalarSolveControls omega_sc{
            2000,sst_linear_tolerance_for_solution_defect(eqw,tolerance),0.7};
        cfdx::core::Vector k_solution(n, 0.0);
        cfdx::core::Vector omega_solution(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            k_solution(i) = k(i);
            omega_solution(i) = omega(i);
        }
        const auto rk = solve_scalar_equation(eqk, k_solution, k_sc);
        const auto rw = solve_scalar_equation(eqw, omega_solution, omega_sc);
        for (std::size_t i = 0; i < n; ++i) {
            k(i) = k_solution(i);
            omega(i) = omega_solution(i);
        }
        record_turbulence_bounds(result, enforce_turbulence_bounds(k,omega,controls));

        double dk=0.0,dw=0.0;
        for(std::size_t i=0;i<n;++i) {
            dk=std::max(dk,std::abs(k(i)-oldk(i)));
            dw=std::max(dw,std::abs(omega(i)-oldw(i)));
            // Residuals must describe the state returned to the caller, after
            // both linear relaxation and turbulence positivity enforcement.
            k_solution(i)=k(i);
            omega_solution(i)=omega(i);
        }
        result.k_residual=scalar_equation_residual_inf(eqk,k_solution);
        result.second_residual=scalar_equation_residual_inf(eqw,omega_solution);
        result.iterations=iter;
        if (result.nonfinite_bound_events==0 && rk.status==cfdx::core::SolverStatus::CONVERGED &&
           rw.status==cfdx::core::SolverStatus::CONVERGED &&
           std::isfinite(result.k_residual) && std::isfinite(result.second_residual) &&
           std::max(dk,dw)<=tolerance &&
           sst_equation_solution_defect(eqk,k_solution)<=tolerance &&
           sst_equation_solution_defect(eqw,omega_solution)<=tolerance) {
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
    const ScalarBoundaryConditions& omega_bcs={},std::size_t max_iterations=100,double tolerance=1e-8,
    const ScalarBoundaryFaceValues* omega_face_values=nullptr)
{
    validate_turbulence_controls(controls);
    const std::size_t n=mesh.n_cells();
    if(n==0||wall_distance.size()!=n||k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces()||
       wall_distance.dimension()!=1||k.dimension()!=1||omega.dimension()!=1||strain_rate.dimension()!=1||mass_flux.dimension()!=1)
        throw std::invalid_argument("SST dynamic-blending field size mismatch");
    if (controls.model!=TurbulenceModel::SST || max_iterations==0 ||
        !std::isfinite(tolerance) || tolerance<=0.0)
        throw std::invalid_argument("SST dynamic blending requires SST and valid iteration controls");
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old_k=k,old_w=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> F1(n,"F1","1",1),F2(n,"F2","1",1);
        // Compute the full-cell gradient contraction once per outer iteration.
        // Recomputing this n-cell vector inside the n-cell loop made this O(n^2).
        const auto cross=sst_boundary_gradient_dot(
            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
        for(std::size_t i=0;i<n;++i){
            const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
            F1(i)=b.first; F2(i)=b.second;
        }
        const auto inner=solve_sst_transport(
            mesh,geometry,mass_flux,k,omega,strain_rate,F1,F2,controls,
            k_bcs,omega_bcs,1,tolerance,omega_face_values);
        double dk=0.0,dw=0.0;
        for(std::size_t i=0;i<n;++i){dk=std::max(dk,std::abs(k(i)-old_k(i)));dw=std::max(dw,std::abs(omega(i)-old_w(i)));}
        accumulate_turbulence_diagnostics(result, inner);
        result.k_residual = inner.k_residual;
        result.second_residual = inner.second_residual;
        result.iterations = iter;
        if (inner.converged && result.nonfinite_bound_events == 0 &&
            std::max(dk, dw) <= tolerance) {
            result.converged = true;
            break;
        }
    }
    return result;
}

} // namespace cfdx::physics
